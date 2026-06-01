#!/usr/bin/env python3
"""
Landau-Zener Parameter Sweep Engine — Fast Edition
====================================================
Drop-in replacement for lz_sweep.py with three key optimizations:

  1. Pure-Numba Dormand-Prince (RK45) adaptive stepper
     — eliminates scipy's 81% Python-overhead per step
  2. Warm-start Newton for lambda_0
     — previous step's result seeds the next, cutting iterations ~3x
  3. Relaxed tolerances (rtol=1e-10, atol=1e-12 default)
     — 44% fewer steps with < 3e-9 error in transition probabilities

Typical speedup: 5-10x over lz_sweep.py (sub-second for well-separated N=3).

State vector layout (real, length 3N-1):
  y[0 : N-1]      = lambda_1,...,lambda_{N-1}   (interior eigenvalues)
  y[N-1 : 2N-1]   = Re(f_0),...,Re(f_{N-1})
  y[2N-1 : 3N-1]  = Im(f_0),...,Im(f_{N-1})

Mathematical reference: LZ_summary.tex (v3)
"""

import numpy as np
from numba import njit
import json
import time
import os
from datetime import datetime


# ==============================================================
#  Newton solver for exterior eigenvalue lambda_0 (warm-startable)
# ==============================================================

@njit(cache=True)
def newton_lambda0(u, gamma_sq, epsilon, n, lam0_hint):
    """
    Clamped Newton-bisection for lambda_0.

    If lam0_hint is finite and on the correct side of the boundary,
    it is used directly (warm start from previous ODE step).
    Otherwise, fall back to the dual-heuristic cold start.

    Returns: (lam0, converged)
    """
    sum_gam2 = 0.0
    for j in range(n):
        sum_gam2 += gamma_sq[j]

    if abs(u) < 1e-300:
        return sum_gam2 / (u + 1e-300), True

    if u < 0.0:
        boundary = epsilon[0]
        # Try warm start
        if np.isfinite(lam0_hint) and lam0_hint < boundary - 1e-14:
            lam0 = lam0_hint
        else:
            # Cold start: dual heuristic
            asymp = sum_gam2 / u
            local = epsilon[0] + gamma_sq[0] / u
            if asymp < boundary - 1e-12:
                lam0 = asymp
            elif local < boundary - 1e-12:
                lam0 = local
            else:
                lam0 = boundary - 1.0
    else:
        boundary = epsilon[n - 1]
        if np.isfinite(lam0_hint) and lam0_hint > boundary + 1e-14:
            lam0 = lam0_hint
        else:
            asymp = sum_gam2 / u
            local = epsilon[n - 1] + gamma_sq[n - 1] / u
            if asymp > boundary + 1e-12:
                lam0 = asymp
            elif local > boundary + 1e-12:
                lam0 = local
            else:
                lam0 = boundary + 1.0

    for it in range(20):
        fval = -u
        fpval = 0.0
        for j in range(n):
            dij = lam0 - epsilon[j]
            fval += gamma_sq[j] / dij
            fpval -= gamma_sq[j] / (dij * dij)

        if abs(fpval) < 1e-30:
            break
        if abs(fval) < 1e-14:
            break

        step = -fval / fpval
        new_lam = lam0 + step

        if u < 0.0 and new_lam >= boundary:
            new_lam = 0.5 * (lam0 + boundary)
        elif u > 0.0 and new_lam <= boundary:
            new_lam = 0.5 * (lam0 + boundary)

        lam0 = new_lam

    return lam0, abs(fval) < 1e-10


# ==============================================================
#  ODE Right-Hand Side (Numba, with warm-start lambda_0)
# ==============================================================

@njit(cache=True)
def lz_rhs(u, y, gamma, gamma_sq, epsilon, a_param, x_init, n,
           lam0_prev):
    """
    RHS of the augmented LZ ODE.  Returns (dydt, lam0) so the
    caller can feed lam0 back as warm start for the next step.
    """
    n_lam = n - 1
    dydt = np.empty(3 * n - 1)

    # Full lambda array
    lam0, _ = newton_lambda0(u, gamma_sq, epsilon, n, lam0_prev)
    lam_all = np.empty(n)
    lam_all[0] = lam0
    for i in range(n_lam):
        lam_all[i + 1] = y[i]

    # Auxiliary sums
    sum_inv_sq = np.zeros(n)
    big_e = np.zeros(n)
    for k in range(n):
        for j in range(n):
            dij = lam_all[k] - epsilon[j]
            dij_sq = dij * dij
            sum_inv_sq[k] += gamma_sq[j] / dij_sq
            big_e[k] += a_param[j] * gamma_sq[j] / dij

    # Sign factors
    delta = np.empty(n)
    for k in range(n):
        delta[k] = 1.0 if gamma[k] > 0.0 else -1.0
    if u >= 0.0:
        delta[0] = -delta[0]

    big_gam = np.empty(n)
    for k in range(n):
        big_gam[k] = delta[k] / np.sqrt(sum_inv_sq[k])

    big_ex = big_e[x_init]

    # lambda' for interior eigenvalues
    for i in range(n_lam):
        dydt[i] = -1.0 / sum_inv_sq[i + 1]

    # f' with gauge and coupling
    for k in range(n):
        f_re_k = y[n_lam + k]
        f_im_k = y[n_lam + n + k]
        alpha_k = big_e[k] - big_ex

        dr = -alpha_k * f_im_k
        di = alpha_k * f_re_k

        for m in range(n):
            if m != k:
                s_km = big_gam[k] * big_gam[m] / (lam_all[k] - lam_all[m])
                dr -= s_km * y[n_lam + m]
                di -= s_km * y[n_lam + n + m]

        dydt[n_lam + k] = dr
        dydt[n_lam + n + k] = di

    return dydt, lam0


# ==============================================================
#  Dormand-Prince RK45 adaptive stepper (pure Numba)
#
#  Butcher tableau coefficients for the DOPRI5 pair (5th order
#  solution + 4th order error estimate).  Adaptive step-sizing
#  with PI controller for smooth step changes.
# ==============================================================

# Dormand-Prince coefficients (as module-level constants for Numba)
_DP_A21 = 1.0/5.0
_DP_A31 = 3.0/40.0;      _DP_A32 = 9.0/40.0
_DP_A41 = 44.0/45.0;     _DP_A42 = -56.0/15.0;     _DP_A43 = 32.0/9.0
_DP_A51 = 19372.0/6561.0; _DP_A52 = -25360.0/2187.0; _DP_A53 = 64448.0/6561.0; _DP_A54 = -212.0/729.0
_DP_A61 = 9017.0/3168.0;  _DP_A62 = -355.0/33.0;    _DP_A63 = 46732.0/5247.0;  _DP_A64 = 49.0/176.0;  _DP_A65 = -5103.0/18656.0
_DP_A71 = 35.0/384.0;     _DP_A73 = 500.0/1113.0;   _DP_A74 = 125.0/192.0;     _DP_A75 = -2187.0/6784.0; _DP_A76 = 11.0/84.0

# Error coefficients: e_i = b_i - b*_i
_DP_E1 = 71.0/57600.0;   _DP_E3 = -71.0/16695.0; _DP_E4 = 71.0/1920.0
_DP_E5 = -17253.0/339200.0; _DP_E6 = 22.0/525.0;  _DP_E7 = -1.0/40.0

# c nodes
_DP_C2 = 1.0/5.0; _DP_C3 = 3.0/10.0; _DP_C4 = 4.0/5.0; _DP_C5 = 8.0/9.0


@njit(cache=True)
def dopri5_integrate(y0, t_start, t_end, gamma, gamma_sq, epsilon,
                     a_param, x_init, n, rtol, atol, max_steps=500000):
    """
    Integrate the LZ ODE from t_start to t_end using Dormand-Prince RK45.

    Returns: (y_final, nfev, nsteps, success)
    """
    ndim = len(y0)
    y = y0.copy()
    t = t_start
    direction = 1.0 if t_end > t_start else -1.0
    dt_max = abs(t_end - t_start)

    # Initial step size estimate
    f0, lam0_prev = lz_rhs(t, y, gamma, gamma_sq, epsilon, a_param,
                            x_init, n, np.inf)
    nfev = 1

    # Compute initial scale
    d0 = 0.0
    d1 = 0.0
    for i in range(ndim):
        sc = atol + rtol * abs(y[i])
        d0 += (y[i] / sc) ** 2
        d1 += (f0[i] / sc) ** 2
    d0 = np.sqrt(d0 / ndim)
    d1 = np.sqrt(d1 / ndim)

    if d0 < 1e-5 or d1 < 1e-5:
        h0 = 1e-6
    else:
        h0 = 0.01 * d0 / d1

    h0 = min(h0, dt_max)
    h = h0 * direction

    # Temp arrays (allocated once)
    k1 = np.empty(ndim)
    k2 = np.empty(ndim)
    k3 = np.empty(ndim)
    k4 = np.empty(ndim)
    k5 = np.empty(ndim)
    k6 = np.empty(ndim)
    k7 = np.empty(ndim)
    y_tmp = np.empty(ndim)

    for i in range(ndim):
        k1[i] = f0[i]

    nsteps = 0
    safety = 0.9
    min_factor = 0.2
    max_factor = 5.0

    while direction * (t_end - t) > 1e-14 * abs(t_end):
        nsteps += 1
        if nsteps > max_steps:
            return y, nfev, nsteps, False

        # Don't overshoot
        if direction * (t + h - t_end) > 0.0:
            h = t_end - t

        dt = h

        # Stage 2
        t2 = t + _DP_C2 * dt
        for i in range(ndim):
            y_tmp[i] = y[i] + dt * _DP_A21 * k1[i]
        f_tmp, lam0_prev = lz_rhs(t2, y_tmp, gamma, gamma_sq, epsilon,
                                    a_param, x_init, n, lam0_prev)
        for i in range(ndim):
            k2[i] = f_tmp[i]

        # Stage 3
        t3 = t + _DP_C3 * dt
        for i in range(ndim):
            y_tmp[i] = y[i] + dt * (_DP_A31 * k1[i] + _DP_A32 * k2[i])
        f_tmp, lam0_prev = lz_rhs(t3, y_tmp, gamma, gamma_sq, epsilon,
                                    a_param, x_init, n, lam0_prev)
        for i in range(ndim):
            k3[i] = f_tmp[i]

        # Stage 4
        t4 = t + _DP_C4 * dt
        for i in range(ndim):
            y_tmp[i] = y[i] + dt * (_DP_A41 * k1[i] + _DP_A42 * k2[i] + _DP_A43 * k3[i])
        f_tmp, lam0_prev = lz_rhs(t4, y_tmp, gamma, gamma_sq, epsilon,
                                    a_param, x_init, n, lam0_prev)
        for i in range(ndim):
            k4[i] = f_tmp[i]

        # Stage 5
        t5 = t + _DP_C5 * dt
        for i in range(ndim):
            y_tmp[i] = y[i] + dt * (_DP_A51 * k1[i] + _DP_A52 * k2[i] + _DP_A53 * k3[i] + _DP_A54 * k4[i])
        f_tmp, lam0_prev = lz_rhs(t5, y_tmp, gamma, gamma_sq, epsilon,
                                    a_param, x_init, n, lam0_prev)
        for i in range(ndim):
            k5[i] = f_tmp[i]

        # Stage 6
        for i in range(ndim):
            y_tmp[i] = y[i] + dt * (_DP_A61 * k1[i] + _DP_A62 * k2[i] + _DP_A63 * k3[i] + _DP_A64 * k4[i] + _DP_A65 * k5[i])
        f_tmp, lam0_prev = lz_rhs(t + dt, y_tmp, gamma, gamma_sq, epsilon,
                                    a_param, x_init, n, lam0_prev)
        for i in range(ndim):
            k6[i] = f_tmp[i]

        # 5th-order solution (stage 7 uses this as y_new)
        for i in range(ndim):
            y_tmp[i] = y[i] + dt * (_DP_A71 * k1[i] + _DP_A73 * k3[i] + _DP_A74 * k4[i] + _DP_A75 * k5[i] + _DP_A76 * k6[i])
        f_tmp, lam0_prev = lz_rhs(t + dt, y_tmp, gamma, gamma_sq, epsilon,
                                    a_param, x_init, n, lam0_prev)
        for i in range(ndim):
            k7[i] = f_tmp[i]

        nfev += 6  # stages 2-7

        # Error estimate
        err_norm = 0.0
        for i in range(ndim):
            sc = atol + rtol * max(abs(y[i]), abs(y_tmp[i]))
            err_i = dt * (_DP_E1 * k1[i] + _DP_E3 * k3[i] + _DP_E4 * k4[i] +
                          _DP_E5 * k5[i] + _DP_E6 * k6[i] + _DP_E7 * k7[i])
            err_norm += (err_i / sc) ** 2
        err_norm = np.sqrt(err_norm / ndim)

        if err_norm <= 1.0:
            # Accept step
            t = t + dt
            for i in range(ndim):
                y[i] = y_tmp[i]
            # FSAL: k7 becomes k1 for next step
            for i in range(ndim):
                k1[i] = k7[i]

            # Step size increase
            if err_norm < 1e-30:
                factor = max_factor
            else:
                factor = min(max_factor, safety * err_norm ** (-0.2))
            h = dt * factor
            # Clamp
            if abs(h) > dt_max:
                h = direction * dt_max
        else:
            # Reject step
            factor = max(min_factor, safety * err_norm ** (-0.2))
            h = dt * factor

    return y, nfev, nsteps, True


# ==============================================================
#  Initial eigenvalue computation (unchanged from lz_sweep.py)
# ==============================================================

def solve_all_lambda(u0, gamma_sq, epsilon):
    """Polynomial root-finding for all N eigenvalues at u0."""
    n = len(epsilon)
    full_prod = np.poly(epsilon)
    poly_coeffs = u0 * full_prod
    for i in range(n):
        eps_without_i = np.delete(epsilon, i)
        partial = np.poly(eps_without_i)
        padded = np.zeros(n + 1)
        padded[1:] = partial
        poly_coeffs -= gamma_sq[i] * padded
    roots = np.roots(poly_coeffs)
    real_roots = np.sort(np.real(roots[np.abs(np.imag(roots)) < 1e-8]))
    if len(real_roots) != n:
        real_roots = np.sort(np.real(roots))
    if u0 < 0:
        return real_roots
    else:
        result = np.zeros(n)
        result[0] = real_roots[-1]
        result[1:] = real_roots[:-1]
        return result


# ==============================================================
#  Initial conditions (unchanged from lz_sweep.py)
# ==============================================================

def compute_initial_state(u0, gamma, epsilon, a_param, x_init):
    """Asymptotic initial conditions at u0 << -1 (Eqs. 5-6 of LZ_summary)."""
    n = len(gamma)
    gamma_sq = gamma ** 2
    lam_vec = solve_all_lambda(u0, gamma_sq, epsilon)
    denom = np.subtract.outer(lam_vec, epsilon)
    sum_inv_sq  = np.sum(gamma_sq / denom**2, axis=1)
    sum_inv_cub = np.sum(gamma_sq / denom**3, axis=1)
    big_e       = np.sum(a_param * gamma_sq / denom, axis=1)
    delta = np.sign(gamma).astype(float)
    delta[delta == 0] = -1.0
    big_gam = delta / np.sqrt(sum_inv_sq)
    lam_diff = np.subtract.outer(lam_vec, lam_vec)
    with np.errstate(divide='ignore', invalid='ignore'):
        s_mat = np.outer(big_gam, big_gam) / lam_diff
    np.fill_diagonal(s_mat, 0.0)
    x = x_init
    f0 = np.zeros(n, dtype=complex)
    for j in range(n):
        if j != x:
            correction = (
                1j / (big_e[x] - big_e[j])
                + (big_gam[x] - big_gam[j]) / (epsilon[x] - epsilon[j])
                - big_gam[x]**2 * sum_inv_cub[x]
                - big_gam[j]**2 * sum_inv_cub[j]
            )
            f0[j] = s_mat[x, j] * correction
    f_sum_sq = np.sum(np.abs(np.delete(f0, x))**2)
    if f_sum_sq >= 1.0:
        f0[x] = 0.0 + 0j
    else:
        f0[x] = np.sqrt(1.0 - f_sum_sq) + 0j
    state = np.zeros(3 * n - 1)
    state[:n-1]         = lam_vec[1:]
    state[n-1 : 2*n-1]  = np.real(f0)
    state[2*n-1 : 3*n-1] = np.imag(f0)
    return state


# ==============================================================
#  Single ODE solve (uses Numba DOPRI5)
# ==============================================================

def solve_lz(gamma, epsilon, a_param, u0, u_final, x_init,
             rtol=1e-10, atol=1e-12):
    """
    Solve the LZ ODE for initial eigenvector x_init (0-indexed).
    Returns: (y_final, nfev, nsteps)
    """
    n = len(gamma)
    gamma_arr = np.asarray(gamma, dtype=np.float64)
    gamma_sq  = gamma_arr ** 2
    epsilon_arr = np.asarray(epsilon, dtype=np.float64)
    a_arr = np.asarray(a_param, dtype=np.float64)

    y0 = compute_initial_state(u0, gamma_arr, epsilon_arr, a_arr, x_init)

    y_final, nfev, nsteps, success = dopri5_integrate(
        y0, u0, u_final, gamma_arr, gamma_sq, epsilon_arr,
        a_arr, x_init, n, rtol, atol)

    if not success:
        raise RuntimeError(f"DOPRI5 failed for x_init={x_init}: max steps exceeded")

    return y_final, nfev, nsteps


# ==============================================================
#  Transition matrix (N solves)
# ==============================================================

def transition_matrix(gamma, epsilon, a_param, u0, u_final,
                      rtol=1e-10, atol=1e-12):
    """
    Compute NxN transition probability matrix.
    P[x, i] = |f_i(u_final)|^2 starting from eigenvector x.
    """
    n = len(gamma)
    n_lam = n - 1
    P = np.zeros((n, n))
    total_nfev = 0
    total_nsteps = 0

    for x in range(n):
        y_final, nfev, nsteps = solve_lz(
            gamma, epsilon, a_param, u0, u_final, x, rtol=rtol, atol=atol)
        total_nfev += nfev
        total_nsteps += nsteps

        for i in range(n):
            re_fi = y_final[n_lam + i]
            im_fi = y_final[n_lam + n + i]
            P[x, i] = re_fi**2 + im_fi**2

    return P, total_nfev, total_nsteps


# ==============================================================
#  Parameter space generation (identical to lz_sweep.py)
# ==============================================================

def generate_parameter_sets():
    """Generate ~100 parameter sets across five physical regimes for N=3."""
    sets = []
    for i in range(20):
        rng = np.random.RandomState(1000 + i)
        g_scale = 0.1 + 0.4 * rng.rand()
        gamma = g_scale * (0.5 + rng.rand(3))
        epsilon = np.sort(-2.0 + 4.0 * rng.rand(3))
        while np.min(np.diff(epsilon)) < 0.3:
            epsilon = np.sort(-2.0 + 4.0 * rng.rand(3))
        a = -1.0 + 2.0 * rng.rand(3)
        sets.append({'gamma': gamma.tolist(), 'epsilon': epsilon.tolist(),
                     'a': a.tolist(), 'regime': 'weak_coupling'})
    for i in range(30):
        rng = np.random.RandomState(2000 + i)
        g_scale = 0.5 + 1.0 * rng.rand()
        gamma = g_scale * (0.5 + rng.rand(3))
        epsilon = np.sort(-2.0 + 4.0 * rng.rand(3))
        while np.min(np.diff(epsilon)) < 0.3:
            epsilon = np.sort(-2.0 + 4.0 * rng.rand(3))
        a = -1.0 + 2.0 * rng.rand(3)
        sets.append({'gamma': gamma.tolist(), 'epsilon': epsilon.tolist(),
                     'a': a.tolist(), 'regime': 'moderate_coupling'})
    for i in range(20):
        rng = np.random.RandomState(3000 + i)
        g_scale = 1.5 + 1.5 * rng.rand()
        gamma = g_scale * (0.5 + rng.rand(3))
        epsilon = np.sort(-2.0 + 4.0 * rng.rand(3))
        while np.min(np.diff(epsilon)) < 0.3:
            epsilon = np.sort(-2.0 + 4.0 * rng.rand(3))
        a = -1.0 + 2.0 * rng.rand(3)
        sets.append({'gamma': gamma.tolist(), 'epsilon': epsilon.tolist(),
                     'a': a.tolist(), 'regime': 'strong_coupling'})
    for i in range(15):
        rng = np.random.RandomState(4000 + i)
        gamma = (0.5 + 1.0 * rng.rand(3)).tolist()
        center = -1.0 + 2.0 * rng.rand()
        gap = 0.01 + 0.05 * rng.rand()
        epsilon = [center - gap, center, center + gap]
        a = (-1.0 + 2.0 * rng.rand(3)).tolist()
        sets.append({'gamma': gamma, 'epsilon': epsilon,
                     'a': a, 'regime': 'near_degenerate'})
    for i in range(15):
        rng = np.random.RandomState(5000 + i)
        gamma = (0.01 + 3.0 * rng.rand(3)).tolist()
        epsilon = np.sort(-3.0 + 6.0 * rng.rand(3))
        while np.min(np.diff(epsilon)) < 0.1:
            epsilon = np.sort(-3.0 + 6.0 * rng.rand(3))
        epsilon = epsilon.tolist()
        a = (-2.0 + 4.0 * rng.rand(3)).tolist()
        sets.append({'gamma': gamma, 'epsilon': epsilon,
                     'a': a, 'regime': 'mixed_asymmetric'})
    return sets


# ==============================================================
#  Sweep driver
# ==============================================================

def run_single(idx, params, u0, u_final, output_dir, rtol=1e-10, atol=1e-12):
    """Solve one parameter set and write JSON output file."""
    gamma   = np.array(params['gamma'])
    epsilon = np.array(params['epsilon'])
    a       = np.array(params['a'])
    n       = len(gamma)

    t0 = time.perf_counter()
    try:
        P, nfev, nsteps = transition_matrix(
            gamma, epsilon, a, u0, u_final, rtol=rtol, atol=atol)
        elapsed = time.perf_counter() - t0

        row_sums = P.sum(axis=1).tolist()
        max_unit_err = max(abs(s - 1.0) for s in row_sums)

        result = {
            'run_index': idx,
            'timestamp': datetime.now().isoformat(),
            'parameters': {
                'N': int(n), 'gamma': params['gamma'],
                'epsilon': params['epsilon'], 'a': params['a'],
                'u0': u0, 'u_final': u_final,
                'regime': params['regime']
            },
            'transition_matrix': P.tolist(),
            'row_sums': row_sums,
            'max_unitarity_error': max_unit_err,
            'compute_time_s': round(elapsed, 6),
            'solver_info': {
                'method': 'DOPRI5_numba', 'rtol': rtol, 'atol': atol,
                'total_nfev': int(nfev), 'total_nsteps': int(nsteps)
            },
            'status': 'OK' if max_unit_err < 1e-8 else 'UNITARITY_WARNING'
        }
    except Exception as e:
        elapsed = time.perf_counter() - t0
        result = {
            'run_index': idx,
            'timestamp': datetime.now().isoformat(),
            'parameters': {
                'N': int(n), 'gamma': params['gamma'],
                'epsilon': params['epsilon'], 'a': params['a'],
                'u0': u0, 'u_final': u_final,
                'regime': params['regime']
            },
            'transition_matrix': None,
            'compute_time_s': round(elapsed, 6),
            'status': f'ERROR: {str(e)}'
        }

    fname = f"run_{idx:03d}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    fpath = os.path.join(output_dir, fname)
    with open(fpath, 'w') as f:
        json.dump(result, f, indent=2)

    return result


def run_sweep(u0=-50.0, u_final=50.0, output_base=None,
              rtol=1e-10, atol=1e-12):
    """Execute the full parameter sweep."""
    if output_base is None:
        output_base = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   'lz_sweep_output')

    ts = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_dir = os.path.join(output_base, f'sweep_fast_{ts}')
    os.makedirs(output_dir, exist_ok=True)

    param_sets = generate_parameter_sets()
    n_sets = len(param_sets)

    print(f"Landau-Zener Parameter Sweep (FAST)")
    print(f"{'=' * 60}")
    print(f"  Parameter sets : {n_sets}")
    print(f"  Domain         : [{u0}, {u_final}]")
    print(f"  Solver         : DOPRI5 (Numba-native), rtol={rtol}, atol={atol}")
    print(f"  Output         : {output_dir}")
    print(f"{'=' * 60}")

    # JIT warm-up
    print("  Warming up JIT...", end=' ', flush=True)
    _g = np.array([1.0, 0.8, 0.6])
    _e = np.array([-1.0, 0.0, 1.5])
    _a = np.array([0.5, -0.3, 0.7])
    _ = transition_matrix(_g, _e, _a, -10.0, 10.0, rtol=rtol, atol=atol)
    print("done.\n")

    results = []
    t_total_start = time.perf_counter()

    for idx, params in enumerate(param_sets):
        result = run_single(idx, params, u0, u_final, output_dir,
                           rtol=rtol, atol=atol)
        results.append(result)

        t = result['compute_time_s']
        stat = result['status']
        regime = params['regime']

        if result.get('transition_matrix') is not None:
            P = np.array(result['transition_matrix'])
            diag = [f"{P[i][i]:.6f}" for i in range(len(P))]
            print(f"  [{idx:3d}/{n_sets}]  {regime:20s}  "
                  f"t={t:7.3f}s  diag=[{', '.join(diag)}]  {stat}")
        else:
            print(f"  [{idx:3d}/{n_sets}]  {regime:20s}  "
                  f"t={t:7.3f}s  {stat}")

    total_wall = time.perf_counter() - t_total_start

    ok   = sum(1 for r in results if r['status'] == 'OK')
    warn = sum(1 for r in results if 'WARNING' in r.get('status', ''))
    err  = sum(1 for r in results if 'ERROR' in r.get('status', ''))
    times = [r['compute_time_s'] for r in results if r.get('transition_matrix')]

    summary = {
        'sweep_timestamp': ts,
        'total_parameter_sets': n_sets,
        'ok': ok, 'warnings': warn, 'errors': err,
        'total_wall_time_s': round(total_wall, 3),
        'mean_solve_time_s': round(float(np.mean(times)), 4) if times else None,
        'median_solve_time_s': round(float(np.median(times)), 4) if times else None,
        'max_solve_time_s': round(float(np.max(times)), 4) if times else None,
        'min_solve_time_s': round(float(np.min(times)), 4) if times else None,
        'solver': {'method': 'DOPRI5_numba', 'rtol': rtol, 'atol': atol,
                   'domain': [u0, u_final]},
        'output_dir': output_dir
    }

    with open(os.path.join(output_dir, 'sweep_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"\n{'=' * 60}")
    print(f"  SWEEP COMPLETE")
    print(f"  OK: {ok}  |  Warnings: {warn}  |  Errors: {err}")
    print(f"  Total wall time : {total_wall:.1f}s")
    if times:
        print(f"  Mean solve time : {np.mean(times):.3f}s")
        print(f"  Median solve    : {np.median(times):.3f}s")
        print(f"  Max solve       : {np.max(times):.3f}s")
    print(f"  Output          : {output_dir}")
    print(f"{'=' * 60}")

    return output_dir, results


if __name__ == '__main__':
    import sys
    base = sys.argv[1] if len(sys.argv) > 1 else None
    run_sweep(output_base=base)
