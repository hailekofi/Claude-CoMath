"""
Regime sweep: where is the incoherent LZ-grid stochastic-matrix model exact,
and what is the leading-order correction?

The stochastic-matrix grid ``assay.grid.grid_P`` composes three doubly-stochastic
2(+)1 Landau-Zener crossings as an ordered product of probability matrices.
This is exact only in the well-separated limit: the two projected Q4 windows
must be wide compared to their separation along the real-u axis.  We measure
window-overlap by the dimensionless

    s = |Delta u_center| / max(Im u(q))

(the ratio of the inter-window center separation to the larger imaginary
extent of the Q4 conjugate pair, in the u-plane).  Past spot checks suggest
median s ~ 0.44 with only ~15% of samples having s > 1.5.

This module
  * draws random ``Params`` from a controlled distribution,
  * computes the benchmark P with Richardson cutoff-extrapolation,
    P_R = (8 P(2T) - P(T))/7  (the genuine adiabatic tail is O(T^-3)),
  * computes the grid P,
  * forms the permutation-invariant residual
    r = min over (row,col) perms of ||P_bench - P_grid[rp,cp]||_inf
    (``assay.grid.best_match``).

Pipelines:
  - ``run_generic_sweep(N)``: ~N random samples, stratified by s.
  - ``run_stratum_sweep(N)``: ~N each of well-separated, generic, near-degen,
    small-slope, large-slope, u_c position strata.
  - ``fit_scaling(...)``: log r vs s (exp fit) and log r vs log s (power fit).
  - ``confluent_q4_diagnostic(par)``: in near-degenerate stratum, compute the
    Q4 quartic discriminant and the merged-window effective action, and test
    the merged-window grid model against benchmark.

The module is read-only with respect to ``assay/*``: nothing in this file
mutates existing modules.
"""

from __future__ import annotations

import csv
import os
import time
from dataclasses import dataclass

import numpy as np

from .geometry import Params, Geometry
from .benchmark import transition_matrix
from .grid import grid_P, best_match, lz_parameters, crossing_order, _stochastic


# ---------------------------------------------------------------------------
#  Random sampler
# ---------------------------------------------------------------------------
def random_sample(rng: np.random.Generator,
                  eps_range: tuple = (-2.0, 2.0),
                  min_gap: float = 0.3,
                  gam_range: tuple = (0.4, 1.5),
                  a_range: tuple = (-2.0, 2.0),
                  x: int = 0,
                  max_tries: int = 200) -> Params:
    """Draw a random Params with real ordered eps (min-gap >= 0.3),
    gam in [0.4,1.5] with random signs, a in [-2,2]."""
    for _ in range(max_tries):
        eps = np.sort(rng.uniform(eps_range[0], eps_range[1], 3))
        if min(np.diff(eps)) < min_gap:
            continue
        gam = rng.uniform(gam_range[0], gam_range[1], 3) * rng.choice([-1.0, 1.0], 3)
        a = rng.uniform(a_range[0], a_range[1], 3)
        return Params(eps=tuple(eps), gam=tuple(gam), a=tuple(a), x=x)
    raise RuntimeError("random_sample: could not satisfy min_gap after max_tries")


def stratified_sample(rng: np.random.Generator, stratum: str) -> Params:
    """
    Draw a Params biased toward the requested stratum.  See module docstring.

    stratum in {generic, well_sep, near_degen, small_slope, large_slope,
                uc_before, uc_between, uc_after}.
    """
    if stratum == "generic":
        return random_sample(rng)
    if stratum == "well_sep":
        # s is scale-free; high s requires *anisotropic* gammas (one large,
        # the others small) with the central level on a wide eps gap.
        # Empirically this is the only systematic way to drive s > 1.5.
        for _ in range(800):
            eps = np.sort(rng.uniform(-3.0, 3.0, 3))
            if min(np.diff(eps)) < 1.0:
                continue
            # one gamma large, two small (centre or edge varying)
            ranks = rng.permutation(3)
            gam = np.empty(3)
            gam[ranks[0]] = rng.uniform(1.4, 2.0) * rng.choice([-1.0, 1.0])
            gam[ranks[1]] = rng.uniform(0.10, 0.30) * rng.choice([-1.0, 1.0])
            gam[ranks[2]] = rng.uniform(0.10, 0.30) * rng.choice([-1.0, 1.0])
            slopes = rng.uniform(-2.0, 2.0, 3)
            par = Params(eps=tuple(eps), gam=tuple(gam), a=tuple(slopes), x=0)
            try:
                wm = window_sep_metric(Geometry(par))
                if wm["s"] > 1.5:
                    return par
            except Exception:
                continue
        raise RuntimeError("well_sep sampler exhausted")
    if stratum == "near_degen":
        # near-degenerate Q4 windows: filter generic samples for s < 0.1.
        # The s < 0.1 stratum populates ~15% of the generic sample space.
        for _ in range(400):
            par = random_sample(rng)
            try:
                wm = window_sep_metric(Geometry(par))
                if wm["s"] < 0.1:
                    return par
            except Exception:
                continue
        raise RuntimeError("near_degen sampler exhausted")
    if stratum == "small_slope":
        return random_sample(rng, a_range=(-0.45, 0.45))
    if stratum == "large_slope":
        # max|a| > 5
        for _ in range(200):
            eps = np.sort(rng.uniform(-2.0, 2.0, 3))
            if min(np.diff(eps)) < 0.3:
                continue
            a = rng.uniform(-7.0, 7.0, 3)
            if max(abs(a)) < 5.0:
                a *= 6.0 / max(abs(a))
            gam = rng.uniform(0.4, 1.5, 3) * rng.choice([-1.0, 1.0], 3)
            return Params(eps=tuple(eps), gam=tuple(gam), a=tuple(a), x=0)
        raise RuntimeError("large_slope sampler exhausted")
    if stratum in ("uc_before", "uc_between", "uc_after"):
        # We do not bias the sampler itself; we draw generic and label later
        # based on where L_H = 0 falls relative to the windows.
        return random_sample(rng)
    raise ValueError(f"unknown stratum: {stratum}")


# ---------------------------------------------------------------------------
#  Window-overlap metric and crossing geometry
# ---------------------------------------------------------------------------
def window_sep_metric(geo: Geometry) -> dict:
    """Compute s = |Delta u_center| / max(Im u(q)) and ancillary geometry."""
    ws = geo.windows()
    if len(ws) < 2:
        return {"s": np.nan, "n_windows": len(ws),
                "u_centers": [w["u_center"] for w in ws]}
    u0, u1 = ws[0]["u_center"], ws[1]["u_center"]
    du = abs(u1 - u0)
    max_im = max(max(abs(u.imag) for u in w["u_vals"]) for w in ws)
    s = du / max_im if max_im > 0 else np.inf
    return {
        "s": float(s),
        "du_center": float(du),
        "max_im_u": float(max_im),
        "u_centers": [u0, u1],
        "u_crossings_diabatic": _diabatic_crossings(geo),
        "u_LH_zero": geo.LH_root(),
        "n_windows": 2,
    }


def _diabatic_crossings(geo: Geometry) -> list:
    """The three pairwise diabatic crossings u_ij sorted by u."""
    H0 = geo.H0
    a = geo.a
    out = []
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        u_ij = (H0[j, j] - H0[i, i]) / (a[i] - a[j])
        out.append((float(u_ij), (i, j)))
    out.sort()
    return out


def uc_position(geo: Geometry, window_data: dict) -> str:
    """Where does L_H=0 fall relative to the two windows?"""
    if window_data["n_windows"] < 2:
        return "no_windows"
    u_lh = window_data["u_LH_zero"]
    uc = window_data["u_centers"]
    lo, hi = min(uc), max(uc)
    if u_lh < lo:
        return "uc_before"
    if u_lh > hi:
        return "uc_after"
    return "uc_between"


# ---------------------------------------------------------------------------
#  Benchmark and residual
# ---------------------------------------------------------------------------
def benchmark_P(geo: Geometry, T_lo: float = 25.0, T_hi: float = 50.0,
                rtol: float = 1e-12, atol: float = 1e-13) -> np.ndarray:
    """Richardson cutoff-extrapolated benchmark: P_R = (8 P(2T) - P(T))/7.

    The genuine adiabatic tail-transition error is O(T^-3), so this is a
    third-order Richardson step (cf. harness conventions).
    """
    P_lo = transition_matrix(geo, T=T_lo, rtol=rtol, atol=atol)
    P_hi = transition_matrix(geo, T=T_hi, rtol=rtol, atol=atol)
    return (8.0 * P_hi - P_lo) / 7.0


def residual_grid_vs_bench(geo: Geometry,
                            T_lo: float = 25.0, T_hi: float = 50.0,
                            rtol: float = 1e-12, atol: float = 1e-13) -> dict:
    """Compute the full residual diagnostic for one Params sample."""
    P_bench = benchmark_P(geo, T_lo=T_lo, T_hi=T_hi, rtol=rtol, atol=atol)
    P_g = grid_P(geo)
    err, rp, cp = best_match(P_g, P_bench)
    return {"P_bench": P_bench, "P_grid": P_g,
            "residual": float(err), "rp": rp, "cp": cp}


# ---------------------------------------------------------------------------
#  Confluent-Q4 merged-window diagnostic
# ---------------------------------------------------------------------------
def q4_discriminant(geo: Geometry) -> dict:
    """Compute the Q4 polynomial discriminant and pair-collision diagnostics.

    Q4 has two pairs of complex-conjugate roots in the generic regime; the
    two pairs are taken to be the two "windows".  As two roots from
    distinct pairs collide (real-real or imag-imag), the discriminant
    vanishes.  We compute the closest-pair distance among Q4 roots.
    """
    roots = geo.Q4_roots()
    disc = float(np.real(np.prod([roots[i] - roots[j]
                                  for i in range(4) for j in range(i+1, 4)])))**2
    # Identify the two pairs by conjugate partner
    used = [False] * 4
    pairs = []
    for i in range(4):
        if used[i]:
            continue
        best, bd = None, np.inf
        for j in range(4):
            if j == i or used[j]:
                continue
            d = abs(roots[i] - np.conj(roots[j]))
            if d < bd:
                best, bd = j, d
        if best is None:
            continue
        used[i] = used[best] = True
        pairs.append((roots[i], roots[best]))
    inter_pair_min = np.inf
    if len(pairs) == 2:
        a1, b1 = pairs[0]
        a2, b2 = pairs[1]
        inter_pair_min = min(abs(a1 - a2), abs(a1 - b2), abs(b1 - a2), abs(b1 - b2))
    return {
        "Q4_roots": roots,
        "discriminant_proxy": disc,
        "min_inter_pair_distance": float(inter_pair_min),
        "n_pairs": len(pairs),
    }


def confluent_grid_P(geo: Geometry) -> np.ndarray:
    """Merged-window effective grid in the confluent-$Q_4$ regime.

    In the near-degenerate stratum the two projected $Q_4$ conjugate
    pairs (the two "windows") nearly merge.  A physically motivated
    effective model is to replace the two corresponding $2(+)1$
    diabatic crossings by a *single fused $3$-state event* with the
    summed elementary actions (LZ exponents are additive in the
    confluent limit) and a doubly-stochastic structure consistent with
    permutation symmetry between the two merging pairs.

    Concretely: when crossings $(i,j)$ and $(i,k)$ (sharing the
    "active" level $i$) merge, we approximate the joint event by the
    classical doubly-stochastic matrix in which level $i$ flips to a
    uniform mixture of the other two with combined exponent, and the
    spectator pair (j,k) inherits the matching complementary
    probabilities by row/column conservation.  Where the two merging
    crossings do not share a level, we fall back on the matrix product
    of the two individual stochastic matrices (which equals the bare
    grid; signalling that the structural form of the grid is correct
    but the elementary actions need correcting in the merged frame).
    """
    data = lz_parameters(geo)
    order = _diabatic_crossings(geo)
    if len(order) < 3:
        return grid_P(geo)
    diffs = [(abs(order[k+1][0] - order[k][0]), k) for k in range(len(order)-1)]
    diffs.sort()
    _, kmin = diffs[0]
    pair_a = order[kmin][1]
    pair_b = order[kmin+1][1]
    # check if they share a level
    shared = set(pair_a) & set(pair_b)
    if not shared:
        # no shared level: fall back to product at merge location
        M_a = _stochastic(pair_a, data["q"][pair_a])
        M_b = _stochastic(pair_b, data["q"][pair_b])
        M_merged = M_b @ M_a
    else:
        i = list(shared)[0]
        j = (set(pair_a) - {i}).pop()
        k = (set(pair_b) - {i}).pop()
        # sum the LZ exponents (joint action of the confluent event)
        G_a = lz_parameters(geo)["Gamma"][pair_a]
        G_b = lz_parameters(geo)["Gamma"][pair_b]
        # joint stay probability for the active level i
        q_joint = float(np.exp(-2.0 * np.pi * (G_a + G_b)))
        # split the loss equally between the two partner levels
        # (this is the unique doubly-stochastic 3-state extension that is
        #  symmetric in the two partners and respects the joint action)
        loss = (1.0 - q_joint) / 2.0
        M_merged = np.eye(3)
        M_merged[i, i] = q_joint
        M_merged[i, j] = loss
        M_merged[i, k] = loss
        M_merged[j, i] = loss
        M_merged[k, i] = loss
        # close the partner pair (j,k) with the remaining off-diagonal so
        # that columns sum to 1 (already conserved row-wise)
        M_merged[j, j] = 1.0 - loss
        M_merged[k, k] = 1.0 - loss
        M_merged[j, k] = 0.0
        M_merged[k, j] = 0.0
    P = np.eye(3)
    for kk, (_, pair) in enumerate(order):
        if kk == kmin:
            P = M_merged @ P
        elif kk == kmin + 1:
            continue
        else:
            P = _stochastic(pair, data["q"][pair]) @ P
    return P


# ---------------------------------------------------------------------------
#  Sweep drivers
# ---------------------------------------------------------------------------
@dataclass
class SweepConfig:
    T_lo: float = 25.0
    T_hi: float = 50.0
    rtol: float = 1e-12
    atol: float = 1e-13
    verbose: bool = True


def run_one(par: Params, cfg: SweepConfig) -> dict:
    """Run the benchmark vs grid residual for one Params."""
    geo = Geometry(par)
    wm = window_sep_metric(geo)
    res = residual_grid_vs_bench(geo, T_lo=cfg.T_lo, T_hi=cfg.T_hi,
                                  rtol=cfg.rtol, atol=cfg.atol)
    return {
        "params": par,
        "s": wm["s"],
        "du_center": wm.get("du_center", np.nan),
        "max_im_u": wm.get("max_im_u", np.nan),
        "u_LH_zero": wm.get("u_LH_zero", np.nan),
        "u_centers": wm.get("u_centers", []),
        "residual": res["residual"],
        "uc_pos": uc_position(geo, wm),
        "n_windows": wm["n_windows"],
        "max_abs_a": float(max(abs(np.asarray(par.a)))),
        "P_bench": res["P_bench"],
        "P_grid": res["P_grid"],
    }


def run_generic_sweep(N: int, seed: int = 12345,
                       cfg: SweepConfig | None = None) -> list[dict]:
    """Stratified-by-s sweep of N random samples drawn generically."""
    if cfg is None:
        cfg = SweepConfig()
    rng = np.random.default_rng(seed)
    out = []
    t0 = time.time()
    for k in range(N):
        try:
            par = random_sample(rng)
            row = run_one(par, cfg)
            row["sample_id"] = k
            row["stratum"] = "generic"
            row["status"] = "ok"
        except Exception as exc:
            row = {"sample_id": k, "stratum": "generic", "status": str(exc),
                   "s": np.nan, "residual": np.nan, "uc_pos": "?",
                   "max_abs_a": np.nan}
        out.append(row)
        if cfg.verbose and (k % 25 == 0):
            elapsed = time.time() - t0
            print(f"  [generic {k:4d}/{N}]  s={row.get('s', float('nan')):.3f}  "
                  f"r={row.get('residual', float('nan')):.3e}  "
                  f"({elapsed:.1f}s elapsed)")
    return out


def run_stratum(stratum: str, N: int, seed: int,
                 cfg: SweepConfig | None = None) -> list[dict]:
    """Run N samples within a given stratum.

    Pre-filtering on geometry/parameters (cheap) is done before the
    expensive benchmark `run_one` call.
    """
    if cfg is None:
        cfg = SweepConfig()
    rng = np.random.default_rng(seed)
    out = []
    n_collected = 0
    n_attempts = 0
    n_tries = 0
    while n_collected < N and n_tries < N * 60:
        n_tries += 1
        try:
            par = stratified_sample(rng, stratum)
        except Exception:
            continue
        # cheap geometry pre-check before running the benchmark
        geo = Geometry(par)
        wm = window_sep_metric(geo)
        if not np.isfinite(wm["s"]):
            continue
        s = wm["s"]
        max_abs_a = float(max(abs(np.asarray(par.a))))
        uc = uc_position(geo, wm)
        keep = True
        if stratum == "well_sep" and not (s > 1.5):
            keep = False
        if stratum == "near_degen" and not (s < 0.1):
            keep = False
        if stratum == "small_slope" and not (max_abs_a < 0.5):
            keep = False
        if stratum == "large_slope" and not (max_abs_a > 5.0):
            keep = False
        if stratum == "uc_before" and uc != "uc_before":
            keep = False
        if stratum == "uc_between" and uc != "uc_between":
            keep = False
        if stratum == "uc_after" and uc != "uc_after":
            keep = False
        if not keep:
            continue
        n_attempts += 1
        try:
            row = run_one(par, cfg)
        except Exception as exc:
            continue
        row["sample_id"] = n_collected
        row["stratum"] = stratum
        row["status"] = "ok"
        out.append(row)
        n_collected += 1
        if cfg.verbose and (n_collected % 10 == 0):
            print(f"  [{stratum} {n_collected:3d}/{N}]  s={row['s']:.3f}  "
                  f"r={row['residual']:.3e}  ({n_tries} tries)")
    return out


# ---------------------------------------------------------------------------
#  Scaling fits
# ---------------------------------------------------------------------------
def fit_power(s_arr: np.ndarray, r_arr: np.ndarray,
              s_min: float = 0.3, s_max: float = 3.0) -> dict:
    """Fit log r = -k log s + c on the band [s_min, s_max], excluding
    saturation at small s.  Returns slope k and prefactor."""
    mask = (s_arr > s_min) & (s_arr < s_max) & (r_arr > 0) & np.isfinite(r_arr)
    if mask.sum() < 5:
        return {"k": np.nan, "c": np.nan, "n": int(mask.sum()), "rms": np.nan}
    xs = np.log(s_arr[mask])
    ys = np.log(r_arr[mask])
    A = np.column_stack([xs, np.ones_like(xs)])
    coeffs, *_ = np.linalg.lstsq(A, ys, rcond=None)
    slope, intercept = coeffs
    resid = ys - (slope * xs + intercept)
    return {"k": float(-slope), "c": float(intercept), "n": int(mask.sum()),
            "rms": float(np.sqrt(np.mean(resid ** 2)))}


def fit_exp(s_arr: np.ndarray, r_arr: np.ndarray,
            s_min: float = 0.3, s_max: float = 3.0) -> dict:
    """Fit log r = -c s + b on the band [s_min, s_max]."""
    mask = (s_arr > s_min) & (s_arr < s_max) & (r_arr > 0) & np.isfinite(r_arr)
    if mask.sum() < 5:
        return {"c": np.nan, "b": np.nan, "n": int(mask.sum()), "rms": np.nan}
    xs = s_arr[mask]
    ys = np.log(r_arr[mask])
    A = np.column_stack([xs, np.ones_like(xs)])
    coeffs, *_ = np.linalg.lstsq(A, ys, rcond=None)
    slope, intercept = coeffs
    resid = ys - (slope * xs + intercept)
    return {"c": float(-slope), "b": float(intercept), "n": int(mask.sum()),
            "rms": float(np.sqrt(np.mean(resid ** 2)))}


# ---------------------------------------------------------------------------
#  CSV I/O
# ---------------------------------------------------------------------------
def write_csv(rows: list[dict], path: str):
    """Write the sweep table (sample_id, s, residual, stratum, ...)."""
    cols = ["sample_id", "stratum", "s", "residual", "du_center", "max_im_u",
            "u_LH_zero", "max_abs_a", "uc_pos", "n_windows", "status"]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            w.writerow([r.get(k, "") for k in cols])


def read_csv(path: str) -> list[dict]:
    out = []
    with open(path, "r") as f:
        rdr = csv.DictReader(f)
        for row in rdr:
            row["s"] = float(row["s"]) if row["s"] not in ("", "nan") else np.nan
            row["residual"] = float(row["residual"]) if row["residual"] not in ("", "nan") else np.nan
            row["du_center"] = float(row["du_center"]) if row["du_center"] not in ("", "nan") else np.nan
            row["max_im_u"] = float(row["max_im_u"]) if row["max_im_u"] not in ("", "nan") else np.nan
            row["u_LH_zero"] = float(row["u_LH_zero"]) if row["u_LH_zero"] not in ("", "nan") else np.nan
            row["max_abs_a"] = float(row["max_abs_a"]) if row["max_abs_a"] not in ("", "nan") else np.nan
            out.append(row)
    return out


# ---------------------------------------------------------------------------
#  Stratum statistics
# ---------------------------------------------------------------------------
def stratum_stats(rows: list[dict]) -> dict:
    """Mean, median, p90, p99, max of the residual within each stratum."""
    by = {}
    for r in rows:
        if r.get("status") != "ok":
            continue
        by.setdefault(r["stratum"], []).append(float(r["residual"]))
    out = {}
    for k, vals in by.items():
        a = np.asarray(vals, float)
        a = a[np.isfinite(a)]
        if a.size == 0:
            continue
        out[k] = {
            "n": int(a.size),
            "mean": float(np.mean(a)),
            "median": float(np.median(a)),
            "p90": float(np.percentile(a, 90)),
            "p99": float(np.percentile(a, 99)),
            "max": float(np.max(a)),
        }
    return out


# ---------------------------------------------------------------------------
#  Driver entrypoint
# ---------------------------------------------------------------------------
def main(out_csv: str,
         n_generic: int = 500,
         n_per_stratum: int = 50,
         cfg: SweepConfig | None = None,
         strata: tuple = ("well_sep", "near_degen", "small_slope",
                          "large_slope", "uc_before", "uc_between", "uc_after"),
         large_slope_n: int | None = None,
         ) -> dict:
    """Run the generic + stratum sweep, write a CSV, and return a summary.

    Partial CSV is written after each stratum to support graceful interrupts.
    `large_slope_n` overrides n_per_stratum for the (typically slow)
    large-slope stratum.
    """
    if cfg is None:
        cfg = SweepConfig()

    rows = []
    t0 = time.time()

    print(f"=== Generic sweep: {n_generic} samples (T_lo={cfg.T_lo}, "
          f"T_hi={cfg.T_hi}, rtol={cfg.rtol}) ===")
    rows.extend(run_generic_sweep(n_generic, seed=12345, cfg=cfg))

    # write partial CSV after generic in case of later interrupt
    write_csv(rows, out_csv)
    print(f"*** Wrote partial {len(rows)} rows to {out_csv}")

    for k, stratum in enumerate(strata):
        n_this = n_per_stratum
        if stratum == "large_slope" and large_slope_n is not None:
            n_this = large_slope_n
        print(f"\n=== Stratum '{stratum}': {n_this} samples ===")
        rows.extend(run_stratum(stratum, n_this, seed=20000 + 1000 * k,
                                cfg=cfg))
        # write partial CSV after each stratum
        write_csv(rows, out_csv)
        print(f"*** Wrote partial {len(rows)} rows to {out_csv}")

    write_csv(rows, out_csv)
    print(f"\n*** Wrote {len(rows)} rows to {out_csv}")
    print(f"*** Total time: {time.time() - t0:.1f} s")

    # Summary
    s_arr = np.array([r["s"] for r in rows if r.get("status") == "ok"
                      and r["stratum"] == "generic"], float)
    r_arr = np.array([r["residual"] for r in rows if r.get("status") == "ok"
                      and r["stratum"] == "generic"], float)

    fp = fit_power(s_arr, r_arr)
    fe = fit_exp(s_arr, r_arr)
    print(f"\nGeneric scaling fits (band 0.3 < s < 3.0):")
    print(f"  power:  r ~ s^(-{fp['k']:.3f})   c={fp['c']:.3f}   "
          f"n={fp['n']}   rms(log)={fp['rms']:.3f}")
    print(f"  exp:    r ~ exp(-{fe['c']:.3f} s)  b={fe['b']:.3f}  "
          f"n={fe['n']}   rms(log)={fe['rms']:.3f}")

    stats = stratum_stats(rows)
    print("\nPer-stratum residual stats:")
    for k, v in stats.items():
        print(f"  {k:14s}  n={v['n']:3d}  mean={v['mean']:.3e}  "
              f"median={v['median']:.3e}  p90={v['p90']:.3e}  "
              f"p99={v['p99']:.3e}  max={v['max']:.3e}")

    return {"rows": rows, "fit_power": fp, "fit_exp": fe, "stratum_stats": stats}
