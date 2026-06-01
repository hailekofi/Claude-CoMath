#!/usr/bin/env python3
"""
Landau-Zener Parameter Sweep Engine — Interaction-Picture Edition
=================================================================
Builds on lz_sweep_fast.py with key optimizations:

  4. Interaction-picture (IP) transformation for the wavefunction ODE
     — factors out the fast oscillatory phase exp(i*int(E_i - E_x)du)
     — the transformed amplitudes g_i are slowly varying: g_i' ~ O(1/|u|)
     — stepper takes ~30x larger steps in the tails (|u| >> 1)
     — transition probabilities: |f_i|^2 = |g_i|^2  (exact, no back-transform)

Physics:
  Original:   f_i' = i(E_i - E_x) f_i  -  Sum_{j!=i} S_{ij} f_j
  IP form:    g_i' = -Sum_{j!=i} S_{ij} exp(i(Phi_j - Phi_i)) g_j
  where:      Phi_i' = E_i - E_x    (accumulated phase, smooth and slow)
              g_i = f_i * exp(-i*Phi_i)

  The diagonal oscillation i(E_i - E_x)f_i has been absorbed into the
  phase factor.  What remains is the coupling S_{ij} ~ O(1/|u|), which
  is the only term that changes transition probabilities.

State vector (real, length 4N-1):
  y[0 : N-1]          = lambda_1,...,lambda_{N-1}    (interior eigenvalues)
  y[N-1 : 2N-1]       = Phi_0,...,Phi_{N-1}          (accumulated phases)
  y[2N-1 : 3N-1]      = Re(g_0),...,Re(g_{N-1})      (IP amplitudes, real)
  y[3N-1 : 4N-1]      = Im(g_0),...,Im(g_{N-1})      (IP amplitudes, imag)

  5. Hybrid initial-condition strategy
     — fixed domain u0=-50, u_final=+50 for all runs
     — when the first-order asymptotic IC (Eq. 5) is valid (Q < 0.01),
       use the analytic formula directly ("formula" path)
     — when Q >= 0.01 (near-degenerate eigenvalues), use backward
       bootstrap: zeroth-order IC (g_x=1, others=0) at u_far=-100,
       backward integrate to u0=-50, then forward integrate to u_final
     — the bootstrap is self-consistent to ~1e-8 and adds < 0.1s
     — preserves subsecond compute time for all 100 parameter sets

Retained from lz_sweep_fast.py:
  1. Pure-Numba Dormand-Prince stepper (no scipy overhead)
  2. Warm-start Newton for lambda_0
  3. Relaxed tolerances (rtol=1e-10, atol=1e-12)
"""

import numpy as np
from numba import njit
import json
import time
import os
from datetime import datetime


# ==============================================================
#  Newton solver (identical to lz_sweep_fast.py)
# ==============================================================

@njit(cache=True)
def newton_lambda0(u, gamma_sq, epsilon, n, lam0_hint):
    """Clamped Newton-bisection for exterior eigenvalue lambda_0."""
    sum_gam2 = 0.0
    for j in range(n):
        sum_gam2 += gamma_sq[j]
    if abs(u) < 1e-300:
        return sum_gam2 / (u + 1e-300), True
    if u < 0.0:
        boundary = epsilon[0]
        if np.isfinite(lam0_hint) and lam0_hint < boundary - 1e-14:
            lam0 = lam0_hint
        else:
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
#  Interaction-picture RHS (Numba)
#
#  State: [lam_1..lam_{N-1}, Phi_0..Phi_{N-1}, Re(g)_0..Re(g)_{N-1}, Im(g)_0..Im(g)_{N-1}]
#  Length: 4N - 1
#
#  Derivatives:
#    lam_i'  = -1 / Sigma^(2)_i         (interior eigenvalues, same as before)
#    Phi_k'  = E_k - E_x                (phase accumulation, smooth)
#    g_k'    = -Sum_{m!=k} S_{km} * exp(i*(Phi_m - Phi_k)) * g_m
#            = -Sum_{m!=k} S_{km} * [cos(dPhi)*Re(g_m) - sin(dPhi)*Im(g_m)]
#              -i*Sum_{m!=k} S_{km} * [cos(dPhi)*Im(g_m) + sin(dPhi)*Re(g_m)]
#
#  The coupling S_{km} = Gamma_k * Gamma_m / (lambda_k - lambda_m)
#  decays as O(1/|u|) in the tails, so g_k' -> 0 for large |u|.
# ==============================================================

@njit(cache=True)
def lz_rhs_ip(u, y, gamma, gamma_sq, epsilon, a_param, x_init, n,
              lam0_prev):
    """
    RHS in the interaction picture.
    Returns (dydt, lam0) for warm-start chaining.
    """
    n_lam = n - 1
    ndim = 4 * n - 1
    dydt = np.empty(ndim)

    # Recover lambda_0 from Newton
    lam0, _ = newton_lambda0(u, gamma_sq, epsilon, n, lam0_prev)
    lam_all = np.empty(n)
    lam_all[0] = lam0
    for i in range(n_lam):
        lam_all[i + 1] = y[i]

    # Extract phases and IP amplitudes from state
    phi = np.empty(n)
    for k in range(n):
        phi[k] = y[n_lam + k]

    g_re = np.empty(n)
    g_im = np.empty(n)
    for k in range(n):
        g_re[k] = y[2 * n - 1 + k]
        g_im[k] = y[3 * n - 1 + k]

    # Auxiliary sums: Sigma^(2), E
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

    # Gamma_k
    big_gam = np.empty(n)
    for k in range(n):
        big_gam[k] = delta[k] / np.sqrt(sum_inv_sq[k])

    big_ex = big_e[x_init]

    # ---- lambda' (interior eigenvalues) ----
    for i in range(n_lam):
        dydt[i] = -1.0 / sum_inv_sq[i + 1]

    # ---- Phi' = E_k - E_x ----
    for k in range(n):
        dydt[n_lam + k] = big_e[k] - big_ex

    # ---- g' in interaction picture ----
    # g_k' = -Sum_{m!=k} S_{km} * exp(i*(Phi_m - Phi_k)) * g_m
    for k in range(n):
        dr = 0.0
        di = 0.0
        for m in range(n):
            if m != k:
                s_km = big_gam[k] * big_gam[m] / (lam_all[k] - lam_all[m])
                d_phi = phi[m] - phi[k]
                cos_dp = np.cos(d_phi)
                sin_dp = np.sin(d_phi)
                # exp(i*dPhi) * g_m = (cos + i*sin) * (Re + i*Im)
                #   Real part: cos*Re - sin*Im
                #   Imag part: cos*Im + sin*Re
                rot_re = cos_dp * g_re[m] - sin_dp * g_im[m]
                rot_im = cos_dp * g_im[m] + sin_dp * g_re[m]
                dr -= s_km * rot_re
                di -= s_km * rot_im

        dydt[2 * n - 1 + k] = dr
        dydt[3 * n - 1 + k] = di

    return dydt, lam0


# ==============================================================
#  Dormand-Prince RK45 (identical stepper, works with any RHS)
# ==============================================================

_DP_A21 = 1.0/5.0
_DP_A31 = 3.0/40.0;      _DP_A32 = 9.0/40.0
_DP_A41 = 44.0/45.0;     _DP_A42 = -56.0/15.0;     _DP_A43 = 32.0/9.0
_DP_A51 = 19372.0/6561.0; _DP_A52 = -25360.0/2187.0; _DP_A53 = 64448.0/6561.0; _DP_A54 = -212.0/729.0
_DP_A61 = 9017.0/3168.0;  _DP_A62 = -355.0/33.0;    _DP_A63 = 46732.0/5247.0;  _DP_A64 = 49.0/176.0;  _DP_A65 = -5103.0/18656.0
_DP_A71 = 35.0/384.0;     _DP_A73 = 500.0/1113.0;   _DP_A74 = 125.0/192.0;     _DP_A75 = -2187.0/6784.0; _DP_A76 = 11.0/84.0
_DP_E1 = 71.0/57600.0;   _DP_E3 = -71.0/16695.0; _DP_E4 = 71.0/1920.0
_DP_E5 = -17253.0/339200.0; _DP_E6 = 22.0/525.0;  _DP_E7 = -1.0/40.0
_DP_C2 = 1.0/5.0; _DP_C3 = 3.0/10.0; _DP_C4 = 4.0/5.0; _DP_C5 = 8.0/9.0


@njit(cache=True)
def dopri5_ip(y0, t_start, t_end, gamma, gamma_sq, epsilon,
              a_param, x_init, n, rtol, atol, max_steps=500000):
    """
    DOPRI5 integrator for the interaction-picture LZ ODE.
    Returns: (y_final, nfev, nsteps, success)
    """
    ndim = len(y0)
    y = y0.copy()
    t = t_start
    direction = 1.0 if t_end > t_start else -1.0
    dt_max = abs(t_end - t_start)

    f0, lam0_prev = lz_rhs_ip(t, y, gamma, gamma_sq, epsilon, a_param,
                                x_init, n, np.inf)
    nfev = 1

    d0 = 0.0; d1 = 0.0
    for i in range(ndim):
        sc = atol + rtol * abs(y[i])
        d0 += (y[i] / sc) ** 2
        d1 += (f0[i] / sc) ** 2
    d0 = np.sqrt(d0 / ndim)
    d1 = np.sqrt(d1 / ndim)
    h0 = 0.01 * d0 / d1 if (d0 >= 1e-5 and d1 >= 1e-5) else 1e-6
    h0 = min(h0, dt_max)
    h = h0 * direction

    k1 = np.empty(ndim); k2 = np.empty(ndim); k3 = np.empty(ndim)
    k4 = np.empty(ndim); k5 = np.empty(ndim); k6 = np.empty(ndim)
    k7 = np.empty(ndim); y_tmp = np.empty(ndim)
    for i in range(ndim): k1[i] = f0[i]

    nsteps = 0

    while direction * (t_end - t) > 1e-14 * abs(t_end):
        nsteps += 1
        if nsteps > max_steps:
            return y, nfev, nsteps, False

        if direction * (t + h - t_end) > 0.0:
            h = t_end - t
        dt = h

        # Stage 2
        for i in range(ndim):
            y_tmp[i] = y[i] + dt * _DP_A21 * k1[i]
        f_tmp, lam0_prev = lz_rhs_ip(t + _DP_C2*dt, y_tmp, gamma, gamma_sq,
                                       epsilon, a_param, x_init, n, lam0_prev)
        for i in range(ndim): k2[i] = f_tmp[i]

        # Stage 3
        for i in range(ndim):
            y_tmp[i] = y[i] + dt * (_DP_A31*k1[i] + _DP_A32*k2[i])
        f_tmp, lam0_prev = lz_rhs_ip(t + _DP_C3*dt, y_tmp, gamma, gamma_sq,
                                       epsilon, a_param, x_init, n, lam0_prev)
        for i in range(ndim): k3[i] = f_tmp[i]

        # Stage 4
        for i in range(ndim):
            y_tmp[i] = y[i] + dt * (_DP_A41*k1[i] + _DP_A42*k2[i] + _DP_A43*k3[i])
        f_tmp, lam0_prev = lz_rhs_ip(t + _DP_C4*dt, y_tmp, gamma, gamma_sq,
                                       epsilon, a_param, x_init, n, lam0_prev)
        for i in range(ndim): k4[i] = f_tmp[i]

        # Stage 5
        for i in range(ndim):
            y_tmp[i] = y[i] + dt * (_DP_A51*k1[i] + _DP_A52*k2[i] + _DP_A53*k3[i] + _DP_A54*k4[i])
        f_tmp, lam0_prev = lz_rhs_ip(t + _DP_C5*dt, y_tmp, gamma, gamma_sq,
                                       epsilon, a_param, x_init, n, lam0_prev)
        for i in range(ndim): k5[i] = f_tmp[i]

        # Stage 6
        for i in range(ndim):
            y_tmp[i] = y[i] + dt * (_DP_A61*k1[i] + _DP_A62*k2[i] + _DP_A63*k3[i] + _DP_A64*k4[i] + _DP_A65*k5[i])
        f_tmp, lam0_prev = lz_rhs_ip(t + dt, y_tmp, gamma, gamma_sq, epsilon,
                                       a_param, x_init, n, lam0_prev)
        for i in range(ndim): k6[i] = f_tmp[i]

        # Stage 7 (5th-order solution)
        for i in range(ndim):
            y_tmp[i] = y[i] + dt * (_DP_A71*k1[i] + _DP_A73*k3[i] + _DP_A74*k4[i] + _DP_A75*k5[i] + _DP_A76*k6[i])
        f_tmp, lam0_prev = lz_rhs_ip(t + dt, y_tmp, gamma, gamma_sq, epsilon,
                                       a_param, x_init, n, lam0_prev)
        for i in range(ndim): k7[i] = f_tmp[i]
        nfev += 6

        # Error estimate
        err_norm = 0.0
        for i in range(ndim):
            sc = atol + rtol * max(abs(y[i]), abs(y_tmp[i]))
            err_i = dt * (_DP_E1*k1[i] + _DP_E3*k3[i] + _DP_E4*k4[i] +
                          _DP_E5*k5[i] + _DP_E6*k6[i] + _DP_E7*k7[i])
            err_norm += (err_i / sc) ** 2
        err_norm = np.sqrt(err_norm / ndim)

        if err_norm <= 1.0:
            t = t + dt
            for i in range(ndim): y[i] = y_tmp[i]
            for i in range(ndim): k1[i] = k7[i]
            if err_norm < 1e-30:
                factor = 5.0
            else:
                factor = min(5.0, 0.9 * err_norm ** (-0.2))
            h = dt * factor
            if abs(h) > dt_max:
                h = direction * dt_max
        else:
            factor = max(0.2, 0.9 * err_norm ** (-0.2))
            h = dt * factor

    return y, nfev, nsteps, True


# ==============================================================
#  Initial eigenvalue computation
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
#  Adaptive u0 selection
#
#  The first-order asymptotic IC (Eq. 5) is valid when the off-
#  diagonal amplitudes are small: sum_{j!=x} |f0_j|^2 << 1.
#  We define the quality metric as:
#
#    Q = max_x sum_{j!=x} |f0_j|^2
#
#  and require Q < threshold before starting the ODE.
#
#  When Q is above threshold, the cubic auxiliary sums Sigma^(3)
#  dominate (eigenvalues too close to bare levels).  Doubling |u0|
#  separates them and Q drops as O(1/u0^2).
#
#  The metric is cheap (~300 us): just eigenvalue solve + Eq. 5.
#  The IP solver makes the extra integration domain ~free.
# ==============================================================

def ic_quality_metric(u0, gamma, epsilon, a_param):
    """
    Compute the asymptotic IC quality metric (no ODE solve).

    Returns:
      max_metric: max over x_init of sum_{j!=x} |f0_j|^2
      per_x:      list of per-x_init values
    """
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

    max_metric = 0.0
    per_x = []
    for x in range(n):
        f0_sq_sum = 0.0
        for j in range(n):
            if j != x:
                correction = (
                    1j / (big_e[x] - big_e[j])
                    + (big_gam[x] - big_gam[j]) / (epsilon[x] - epsilon[j])
                    - big_gam[x]**2 * sum_inv_cub[x]
                    - big_gam[j]**2 * sum_inv_cub[j]
                )
                f0_sq_sum += abs(s_mat[x, j] * correction)**2
        per_x.append(f0_sq_sum)
        max_metric = max(max_metric, f0_sq_sum)

    return max_metric, per_x


def compute_bootstrap_ic(u0, u_far, gamma, epsilon, a_param, x_init,
                         rtol=1e-10, atol=1e-12):
    """
    Backward-bootstrap initial condition for cases where the analytic
    first-order IC fails (Q >= threshold).

    Strategy:
      1. At u_far (e.g. -100), set zeroth-order IC: g_x = 1, others = 0.
         This is exact as u -> -inf since coupling S ~ O(1/|u|) -> 0.
      2. Backward-integrate from u_far to u0 (e.g. -100 -> -50).
         The ODE dynamics in this range are perturbatively slow
         (max|S_ij| ~ 0.004 at u=-50), so the backward integration
         self-consistently builds up the small off-diagonal components.
      3. The resulting state at u0 is used as IC for the forward solve.

    Self-consistency: varying u_far from -100 to -500 changes transition
    probs by < 1e-8.  The bootstrap is much more accurate than pushing
    the analytic formula, which can have Q > 1 in degenerate cases.

    Returns: state vector at u0 (length 4N-1)
    """
    n = len(gamma)
    gamma_arr = np.asarray(gamma, dtype=np.float64)
    gamma_sq = gamma_arr ** 2
    epsilon_arr = np.asarray(epsilon, dtype=np.float64)
    a_arr = np.asarray(a_param, dtype=np.float64)

    # Zeroth-order IC at u_far: g_x = 1, all others = 0, phases = 0
    lam_far = solve_all_lambda(u_far, gamma_sq, epsilon_arr)
    state = np.zeros(4 * n - 1)
    state[:n-1] = lam_far[1:]                     # interior eigenvalues
    # phases all zero
    state[2 * n - 1 + x_init] = 1.0               # Re(g_x) = 1
    # Im(g_x) = 0, all other g = 0

    # Backward integrate: u_far -> u0 (note: u_far < u0 < 0, so this
    # is a forward-in-time integration from u_far to u0)
    y_at_u0, nfev_b, nsteps_b, success_b = dopri5_ip(
        state, u_far, u0, gamma_arr, gamma_sq, epsilon_arr,
        a_arr, x_init, n, rtol, atol, max_steps=100000)

    if not success_b:
        raise RuntimeError(
            "Bootstrap backward integration failed for x_init=%d" % x_init)

    return y_at_u0, nfev_b, nsteps_b


# ==============================================================
#  Initial conditions in IP variables
# ==============================================================

def compute_initial_state_ip(u0, gamma, epsilon, a_param, x_init):
    """
    Compute initial state in interaction-picture variables.

    At u0, the phases Phi_i are all zero (by definition: Phi is
    accumulated from u0 onward).  So g_i(u0) = f_i(u0) — the
    initial amplitudes are identical to the lab frame.

    State: [lam_1..lam_{N-1}, Phi_0..Phi_{N-1}, Re(g), Im(g)]
    Length: 4N - 1
    """
    n = len(gamma)
    gamma_sq = gamma ** 2
    lam_vec = solve_all_lambda(u0, gamma_sq, epsilon)

    # Compute f0 the same way as before
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

    # Pack IP state: [lambdas, phases (all zero), Re(g), Im(g)]
    state = np.zeros(4 * n - 1)
    state[:n-1] = lam_vec[1:]            # interior lambdas
    # state[n-1 : 2*n-1] = 0.0           # phases all zero at u0
    state[2*n-1 : 3*n-1] = np.real(f0)   # g = f at u0 (phases = 0)
    state[3*n-1 : 4*n-1] = np.imag(f0)

    return state


# ==============================================================
#  Single ODE solve (interaction-picture)
# ==============================================================

def solve_lz(gamma, epsilon, a_param, u0, u_final, x_init,
             rtol=1e-10, atol=1e-12, use_bootstrap=False, u_far=-100.0):
    """
    Solve the LZ ODE in the interaction picture.

    If use_bootstrap=True, generates IC via backward bootstrap from u_far.
    Otherwise uses the analytic first-order formula.

    Returns: (y_final, nfev, nsteps, ic_method)
      ic_method: "formula" or "bootstrap"

    Transition probabilities: P[i] = y_final[2N-1+i]^2 + y_final[3N-1+i]^2
    (Since |f_i|^2 = |g_i|^2, no back-transform needed.)
    """
    n = len(gamma)
    gamma_arr = np.asarray(gamma, dtype=np.float64)
    gamma_sq  = gamma_arr ** 2
    epsilon_arr = np.asarray(epsilon, dtype=np.float64)
    a_arr = np.asarray(a_param, dtype=np.float64)

    nfev_total = 0
    nsteps_total = 0

    if use_bootstrap:
        y0, nfev_b, nsteps_b = compute_bootstrap_ic(
            u0, u_far, gamma_arr, epsilon_arr, a_arr, x_init,
            rtol=rtol, atol=atol)
        nfev_total += nfev_b
        nsteps_total += nsteps_b
        ic_method = "bootstrap"
    else:
        y0 = compute_initial_state_ip(u0, gamma_arr, epsilon_arr, a_arr, x_init)
        ic_method = "formula"

    y_final, nfev, nsteps, success = dopri5_ip(
        y0, u0, u_final, gamma_arr, gamma_sq, epsilon_arr,
        a_arr, x_init, n, rtol, atol, max_steps=500000)
    nfev_total += nfev
    nsteps_total += nsteps

    if not success:
        raise RuntimeError("DOPRI5-IP failed for x_init=%d: max steps exceeded" % x_init)

    return y_final, nfev_total, nsteps_total, ic_method


# ==============================================================
#  Transition matrix
# ==============================================================

def transition_matrix(gamma, epsilon, a_param, u0, u_final,
                      rtol=1e-10, atol=1e-12,
                      use_bootstrap=False, u_far=-100.0):
    """
    NxN transition probability matrix.
    P[x, i] = |g_i(u_final)|^2 = |f_i(u_final)|^2

    If use_bootstrap=True, all x_init channels use backward bootstrap IC.
    """
    n = len(gamma)
    P = np.zeros((n, n))
    total_nfev = 0
    total_nsteps = 0
    ic_method = None

    for x in range(n):
        y_final, nfev, nsteps, ic_m = solve_lz(
            gamma, epsilon, a_param, u0, u_final, x,
            rtol=rtol, atol=atol,
            use_bootstrap=use_bootstrap, u_far=u_far)
        total_nfev += nfev
        total_nsteps += nsteps
        ic_method = ic_m

        # |g_i|^2 = |f_i|^2 — phases cancel
        for i in range(n):
            re_gi = y_final[2 * n - 1 + i]
            im_gi = y_final[3 * n - 1 + i]
            P[x, i] = re_gi**2 + im_gi**2

    return P, total_nfev, total_nsteps, ic_method


# ==============================================================
#  Parameter space generation (identical)
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

def run_single(idx, params, u0, u_final, u_far, output_dir,
               rtol=1e-10, atol=1e-12, ic_threshold=0.01):
    """Solve one parameter set with hybrid IC and write JSON."""
    gamma   = np.array(params['gamma'])
    epsilon = np.array(params['epsilon'])
    a       = np.array(params['a'])
    n       = len(gamma)

    # Evaluate IC quality metric to decide formula vs bootstrap
    ic_metric, ic_per_x = ic_quality_metric(u0, gamma, epsilon, a)
    use_bootstrap = (ic_metric >= ic_threshold)

    t0 = time.perf_counter()
    try:
        P, nfev, nsteps, ic_method = transition_matrix(
            gamma, epsilon, a, u0, u_final,
            rtol=rtol, atol=atol,
            use_bootstrap=use_bootstrap, u_far=u_far)
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
            'hybrid_ic': {
                'ic_method': ic_method,
                'ic_metric': round(ic_metric, 10),
                'ic_per_x': [round(v, 10) for v in ic_per_x],
                'ic_threshold': ic_threshold,
                'u_far': u_far if use_bootstrap else None
            },
            'transition_matrix': P.tolist(),
            'row_sums': row_sums,
            'max_unitarity_error': max_unit_err,
            'compute_time_s': round(elapsed, 6),
            'solver_info': {
                'method': 'DOPRI5_IP_numba', 'rtol': rtol, 'atol': atol,
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
            'hybrid_ic': {
                'ic_method': 'bootstrap' if use_bootstrap else 'formula',
                'ic_metric': round(ic_metric, 10),
                'ic_per_x': [round(v, 10) for v in ic_per_x],
                'ic_threshold': ic_threshold,
                'u_far': u_far if use_bootstrap else None
            },
            'transition_matrix': None,
            'compute_time_s': round(elapsed, 6),
            'status': 'ERROR: %s' % str(e)
        }

    fname = "run_%03d_%s.json" % (idx, datetime.now().strftime('%Y%m%d_%H%M%S'))
    fpath = os.path.join(output_dir, fname)
    with open(fpath, 'w') as f:
        json.dump(result, f, indent=2)
    return result


def run_sweep(u0=-50.0, u_final=50.0, u_far=-100.0, output_base=None,
              rtol=1e-10, atol=1e-12, ic_threshold=0.01):
    """Execute the full parameter sweep with hybrid IC strategy."""
    if output_base is None:
        output_base = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   'lz_sweep_output')
    ts = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_dir = os.path.join(output_base, 'sweep_ip_%s' % ts)
    os.makedirs(output_dir, exist_ok=True)

    param_sets = generate_parameter_sets()
    n_sets = len(param_sets)

    print("Landau-Zener Parameter Sweep (IP + Hybrid IC)")
    print("=" * 70)
    print("  Parameter sets  : %d" % n_sets)
    print("  Domain          : [%g, %g]" % (u0, u_final))
    print("  Bootstrap u_far : %g" % u_far)
    print("  IC threshold    : %g" % ic_threshold)
    print("  Solver          : DOPRI5-IP (Numba), rtol=%g, atol=%g" % (rtol, atol))
    print("  Output          : %s" % output_dir)
    print("=" * 70)

    print("  Warming up JIT...", end=' ', flush=True)
    _g = np.array([1.0, 0.8, 0.6])
    _e = np.array([-1.0, 0.0, 1.5])
    _a = np.array([0.5, -0.3, 0.7])
    _ = transition_matrix(_g, _e, _a, -10., 10., rtol=rtol, atol=atol)
    # Also warm up the bootstrap path
    _ = transition_matrix(_g, _e, _a, -10., 10., rtol=rtol, atol=atol,
                          use_bootstrap=True, u_far=-20.0)
    print("done.\n")

    results = []
    n_formula = 0
    n_bootstrap = 0
    t_total_start = time.perf_counter()

    for idx, params in enumerate(param_sets):
        result = run_single(idx, params, u0, u_final, u_far, output_dir,
                           rtol=rtol, atol=atol, ic_threshold=ic_threshold)
        results.append(result)
        t = result['compute_time_s']
        stat = result['status']
        regime = params['regime']
        hic = result['hybrid_ic']
        ic_tag = hic['ic_method'].upper()[0]  # F or B
        if hic['ic_method'] == 'bootstrap':
            n_bootstrap += 1
        else:
            n_formula += 1

        if result.get('transition_matrix') is not None:
            P = np.array(result['transition_matrix'])
            diag = [("%.6f" % P[i][i]) for i in range(len(P))]
            print("  [%3d/%d]  %-20s  t=%7.3fs  Q=%.4f  IC=%s  "
                  "diag=[%s]  %s" % (idx, n_sets, regime, t,
                  hic['ic_metric'], ic_tag, ', '.join(diag), stat))
        else:
            print("  [%3d/%d]  %-20s  t=%7.3fs  Q=%.4f  IC=%s  %s"
                  % (idx, n_sets, regime, t, hic['ic_metric'], ic_tag, stat))

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
        'hybrid_ic': {
            'u0': u0, 'u_final': u_final, 'u_far': u_far,
            'ic_threshold': ic_threshold,
            'n_formula': n_formula,
            'n_bootstrap': n_bootstrap
        },
        'solver': {'method': 'DOPRI5_IP_numba', 'rtol': rtol, 'atol': atol},
        'output_dir': output_dir
    }
    with open(os.path.join(output_dir, 'sweep_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)

    print("\n" + "=" * 70)
    print("  SWEEP COMPLETE")
    print("  OK: %d  |  Warnings: %d  |  Errors: %d" % (ok, warn, err))
    print("  Total wall time : %.1fs" % total_wall)
    if times:
        print("  Mean solve time : %.3fs" % np.mean(times))
        print("  Median solve    : %.3fs" % np.median(times))
        print("  Max solve       : %.3fs" % np.max(times))
    print("  IC methods      : formula=%d, bootstrap=%d" % (n_formula, n_bootstrap))
    print("  Output          : %s" % output_dir)
    print("=" * 70)

    return output_dir, results


if __name__ == '__main__':
    import sys
    base = sys.argv[1] if len(sys.argv) > 1 else None
    run_sweep(output_base=base)
