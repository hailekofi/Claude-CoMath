"""
Validation matrix: stratified benchmark dataset and per-stratum residual API.

The Type-1, N=3 closed-form programme submits candidate formulas for the full
3x3 transition matrix P to a single, dense, stratified validation harness.
This module is that harness.

Pipeline
--------
1.  Sample-strata generator.   Eight strata stress different regimes:
      (S1) generic        -- random ordered eps, generic slopes
      (S2) small slope    -- max|a| in [0.1, 0.5]
      (S3) large slope    -- max|a| in [4, 8]
      (S4) near-degenerate Q4 (small conjugate-pair separation/width ratio)
      (S5) near-W4 wall   (one W4 root close to the real u-axis)
      (S6) crossing-before   (L_H=0 crossing-u u_c < min window u_center)
      (S7) crossing-between  (u_c between the two window u_centers)
      (S8) crossing-after    (u_c > max window u_center)

2.  Benchmark.  Diabatic interaction-picture IVP with Richardson
    extrapolation in the cutoff T,
        P_bench  =  (8 P(2T) - P(T)) / 7.
    Default T = 120, 2T = 240, rtol = 1e-13, atol = 1e-14.  Target ~1e-9.

3.  Persisted dataset.  list[Record] saved to ``validation_dataset.pkl``;
    each record carries (eps, gam, a, x, stratum, P_bench, q_ij[3],
    window_periods[2], separation_width_ratio).

4.  Public API.  ``load_dataset()``, ``validate(formula, dataset)`` and
    ``summarize(residuals)``.

5.  Self-tests built in:
      (a) identity formula   P_pred = I,
      (b) Brundobler-Elser only: the two extreme-survival diagonal entries
          from the q_ij products, all other entries 0,
      (c) incoherent LZ grid from ``assay.grid.grid_P``.

All sampling is seeded; running the generator twice yields identical data.

CLI
---
    python -m assay.validation_matrix build [--samples-per-stratum N]
    python -m assay.validation_matrix selftest
    python -m assay.validation_matrix richardson-demo
"""

from __future__ import annotations

import argparse
import os
import pickle
import sys
import time
from dataclasses import dataclass, asdict, field
from typing import Callable, Dict, Iterable, List, Optional, Tuple

import numpy as np

from .geometry import Params, Geometry
from .ip import propagate_ad_ip
from .actions import all_window_actions
from .grid import grid_P, lz_parameters, crossing_order


# ===========================================================================
#  Record dataclass
# ===========================================================================
@dataclass
class Record:
    """
    One stratified validation sample.

    Fields
    ------
    eps, gam, a, x      : raw model parameters.
    stratum             : stratum tag (S1..S8).
    P_bench             : Richardson-extrapolated IP transition matrix (3x3).
    P_T, P_2T           : the two raw IP transitions feeding the Richardson.
    q_ij                : dict[(i,j) -> q_ij = exp(-2 pi Gamma_ij)] for all pairs.
    window_periods      : tuple[complex, complex], I_X for each window in u-order.
    sep_width_ratio     : Q4 root pair separation / window width (min over pairs).
    notes               : free dict for stratum-specific diagnostics.
    """

    eps: Tuple[float, float, float]
    gam: Tuple[float, float, float]
    a: Tuple[float, float, float]
    x: int
    stratum: str
    P_bench: np.ndarray
    P_T: np.ndarray
    P_2T: np.ndarray
    q_ij: Dict[Tuple[int, int], float]
    window_periods: Tuple[complex, complex]
    sep_width_ratio: float
    notes: dict = field(default_factory=dict)


# ===========================================================================
#  Benchmark with Richardson extrapolation
# ===========================================================================
def benchmark_P(geo: Geometry, x: int = 0, T: float = 120.0,
                rtol: float = 1e-13, atol: float = 1e-14
                ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Richardson-extrapolated transition matrix.

    Returns (P_richardson, P_T, P_2T).  P[x,j] = |U[j,x]|^2 where
    U is the gauged adiabatic IP one-pass propagator.
    """
    U_T = propagate_ad_ip(geo, -T, T, x, rtol=rtol, atol=atol)
    P_T = np.abs(U_T.T) ** 2

    U_2T = propagate_ad_ip(geo, -2.0 * T, 2.0 * T, x, rtol=rtol, atol=atol)
    P_2T = np.abs(U_2T.T) ** 2

    P_rich = (8.0 * P_2T - P_T) / 7.0
    return P_rich, P_T, P_2T


# ===========================================================================
#  Auxiliary diagnostics for each sample
# ===========================================================================
def _separation_width_ratio(geo: Geometry) -> float:
    """
    Min over the two conjugate-pair Q4 windows of
        (pair midpoint imaginary spread) / (real separation to nearest pole).
    A small value flags near-degenerate Q4 root windows.
    """
    wins = geo.windows()
    if len(wins) != 2:
        return 0.0
    ratios = []
    for w in wins:
        q = w["roots"]
        sep = float(abs(q[0] - q[1]))          # conjugate-pair separation
        # nearest pole to the pair centre
        centre_re = float(np.real(0.5 * (q[0] + q[1])))
        pole_dist = float(np.min(np.abs(centre_re - geo.eps)))
        if pole_dist <= 0.0:
            ratios.append(0.0)
        else:
            ratios.append(sep / pole_dist)
    return float(min(ratios))


def _window_periods(geo: Geometry) -> Tuple[complex, complex]:
    acts = all_window_actions(geo)
    # u-order
    wins = sorted(geo.windows(), key=lambda w: w["u_center"])
    out = []
    for w in wins:
        # match action by window centre
        best = min(acts, key=lambda a: abs(a["centre"] - 0.5 * (w["roots"][0] + w["roots"][1])))
        out.append(complex(best["I_X"]))
    return out[0], out[1]


def _u_centers(geo: Geometry) -> Tuple[float, float]:
    wins = sorted(geo.windows(), key=lambda w: w["u_center"])
    return wins[0]["u_center"], wins[1]["u_center"]


def _u_crossing(geo: Geometry) -> float:
    """u at which the linear functional L_H(lambda) vanishes, mapped through u(lambda)."""
    lam_star = geo.LH_root()
    p_val = float(np.polyval(geo.p, lam_star))
    n_val = float(np.polyval(geo.n, lam_star))
    if abs(p_val) < 1e-15:
        return float("inf")
    return n_val / p_val


def _w4_real_distance(geo: Geometry) -> float:
    """Min |Im(W4 root)| over the four (or fewer) roots of W4(lambda)."""
    roots = geo.W4_roots()
    if len(roots) == 0:
        return float("inf")
    return float(np.min(np.abs(np.imag(roots))))


# ===========================================================================
#  Stratified sample generator
# ===========================================================================
class StratumRejection(Exception):
    pass


def _draw_one(rng: np.random.Generator, stratum: str, x: int = 0
              ) -> Params:
    """
    Draw a candidate Params satisfying the stratum's constraints.
    Raises StratumRejection if the candidate fails the constraints.
    """
    # ordered eps with min-gap 0.3
    while True:
        e = np.sort(rng.uniform(-1.2, 1.2, size=3))
        if np.min(np.diff(e)) >= 0.3:
            break

    # gamma: random signs, magnitudes in [0.4, 1.5]
    g = rng.uniform(0.4, 1.5, size=3) * rng.choice([-1.0, 1.0], size=3)

    if stratum == "S1":
        a = rng.uniform(-2.0, 2.0, size=3)
        if np.max(np.abs(a)) < 0.5:
            raise StratumRejection("a too small for generic")
    elif stratum == "S2":
        amax = rng.uniform(0.1, 0.5)
        # distribute three slopes within +-amax
        a = rng.uniform(-amax, amax, size=3)
        if np.max(np.abs(a)) < 0.05:
            raise StratumRejection("|a| effectively zero")
    elif stratum == "S3":
        amax = rng.uniform(4.0, 8.0)
        a = rng.uniform(-amax, amax, size=3)
        # ensure at least one slope is near amax in magnitude
        if np.max(np.abs(a)) < 3.5:
            raise StratumRejection("too small for large-slope stratum")
    elif stratum in ("S4", "S5", "S6", "S7", "S8"):
        a = rng.uniform(-2.0, 2.0, size=3)
        if np.max(np.abs(a)) < 0.5:
            raise StratumRejection("too quiet for stress stratum")
    else:
        raise ValueError(f"unknown stratum {stratum!r}")

    # ensure slopes are non-degenerate
    if np.min(np.abs(np.array([a[0]-a[1], a[0]-a[2], a[1]-a[2]]))) < 0.05:
        raise StratumRejection("slopes too close")

    return Params(eps=tuple(map(float, e)),
                  gam=tuple(map(float, g)),
                  a=tuple(map(float, a)),
                  x=x)


def _accept_stratum(par: Params, stratum: str
                    ) -> Tuple[bool, dict]:
    """Apply the post-construction acceptance test for stress strata."""
    geo = Geometry(par)
    wins = geo.windows()
    if len(wins) != 2:
        return False, {"reason": "windows!=2"}
    # Q4 must have two genuinely conjugate pairs (small conj_defect)
    if any(w["conj_defect"] > 1e-6 for w in wins):
        return False, {"reason": "Q4 not 2 conj pairs"}

    sw = _separation_width_ratio(geo)
    w4dist = _w4_real_distance(geo)
    uc = _u_crossing(geo)
    uA, uB = _u_centers(geo)

    info = {"sep_width_ratio": sw, "w4_real_dist": w4dist,
            "u_crossing": uc, "u_centerA": uA, "u_centerB": uB}

    if stratum == "S4":
        if not (sw < 0.15):
            return False, {**info, "reason": "sw>=0.15"}
    elif stratum == "S5":
        if not (w4dist < 0.10):
            return False, {**info, "reason": "no W4 near real axis"}
    elif stratum == "S6":
        if not (uc < min(uA, uB)):
            return False, {**info, "reason": "u_crossing not before"}
    elif stratum == "S7":
        if not (min(uA, uB) < uc < max(uA, uB)):
            return False, {**info, "reason": "u_crossing not between"}
    elif stratum == "S8":
        if not (uc > max(uA, uB)):
            return False, {**info, "reason": "u_crossing not after"}

    return True, info


STRATA = ("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8")
DEFAULT_SAMPLES_PER_STRATUM = 60


# ===========================================================================
#  Single-sample benchmark worker
# ===========================================================================
def _build_record_from_params(par: Params, stratum: str,
                              T: float, rtol: float, atol: float,
                              info: dict) -> Record:
    geo = Geometry(par)
    P_bench, P_T, P_2T = benchmark_P(geo, x=par.x, T=T, rtol=rtol, atol=atol)
    q = lz_parameters(geo)["q"]
    IXa, IXb = _window_periods(geo)
    sw = _separation_width_ratio(geo)
    return Record(
        eps=par.eps, gam=par.gam, a=par.a, x=par.x,
        stratum=stratum,
        P_bench=P_bench, P_T=P_T, P_2T=P_2T,
        q_ij={k: float(v) for k, v in q.items()},
        window_periods=(IXa, IXb),
        sep_width_ratio=sw,
        notes=info,
    )


def _draw_and_benchmark(args):
    """Worker used by multiprocessing.Pool.imap_unordered (legacy path)."""
    seed, stratum, T, rtol, atol = args
    rng = np.random.default_rng(seed)
    # try up to 5000 candidates to satisfy the stratum
    for _ in range(5000):
        try:
            par = _draw_one(rng, stratum)
        except StratumRejection:
            continue
        try:
            ok, info = _accept_stratum(par, stratum)
        except Exception:
            continue
        if not ok:
            continue
        try:
            rec = _build_record_from_params(par, stratum, T, rtol, atol, info)
            return seed, rec, None
        except Exception as e:
            return seed, None, repr(e)
    return seed, None, "no candidate accepted in 5000 tries"


# ---------------------------------------------------------------------------
#  Two-phase build: (i) generate Params upfront (cheap); (ii) bench in parallel
# ---------------------------------------------------------------------------
def _generate_params_for_stratum(stratum: str, n_samples: int,
                                 rng: np.random.Generator,
                                 max_attempts: int = 200_000
                                 ) -> List[Tuple[Params, dict]]:
    """
    Draw ``n_samples`` parameter sets satisfying the stratum's acceptance
    test.  Returns a list of (Params, info_dict) tuples.
    """
    accepted: List[Tuple[Params, dict]] = []
    tried = 0
    while len(accepted) < n_samples and tried < max_attempts:
        tried += 1
        try:
            par = _draw_one(rng, stratum)
        except StratumRejection:
            continue
        try:
            ok, info = _accept_stratum(par, stratum)
        except Exception:
            continue
        if not ok:
            continue
        accepted.append((par, info))
    if len(accepted) < n_samples:
        # we did our best -- caller decides what to do
        pass
    return accepted


def _bench_worker(args):
    """Worker computing the Richardson-extrapolated benchmark for one Params."""
    idx, eps, gam, a, x, stratum, T, rtol, atol, info = args
    try:
        par = Params(eps=tuple(eps), gam=tuple(gam), a=tuple(a), x=int(x))
        rec = _build_record_from_params(par, stratum, T, rtol, atol, dict(info))
        return idx, rec, None
    except Exception as e:
        return idx, None, repr(e)


# ===========================================================================
#  Top-level dataset builder
# ===========================================================================
DATASET_PATH = os.path.join(os.path.dirname(__file__), "validation_dataset.pkl")


def build_dataset(samples_per_stratum: int = DEFAULT_SAMPLES_PER_STRATUM,
                  T: float = 120.0,
                  rtol: float = 1e-13, atol: float = 1e-14,
                  seed: int = 20260523,
                  nproc: Optional[int] = None,
                  out_path: str = DATASET_PATH,
                  verbose: bool = True,
                  progress_path: Optional[str] = None) -> List[Record]:
    """
    Build the stratified validation dataset and persist it.

    The seed determines all sample parameters deterministically.  The build is
    two-phase:
      (1) generate ``samples_per_stratum`` Params per stratum on the main
          process (cheap rejection sampling).
      (2) run the Richardson-extrapolated IP benchmark in parallel
          (``multiprocessing.Pool``).
    """
    import multiprocessing as mp

    # ---- Phase 1: generate all Params upfront (sequential, fast) -----
    if verbose:
        print(f"[validation_matrix] phase 1: drawing "
              f"{samples_per_stratum} samples per stratum",
              flush=True)
    plan: List[Tuple[Params, str, dict]] = []
    short = {}
    for si, stratum in enumerate(STRATA):
        # deterministic per-stratum sub-seed (Python's hash() is salt-randomised
        # by default, so derive a stable seed from (seed, stratum-index))
        sub = np.random.default_rng(int(seed) + 1_000_000 * (si + 1))
        accepted = _generate_params_for_stratum(stratum, samples_per_stratum, sub)
        if verbose:
            print(f"   {stratum}: drew {len(accepted)}/{samples_per_stratum}",
                  flush=True)
        if len(accepted) < samples_per_stratum:
            short[stratum] = samples_per_stratum - len(accepted)
        for par, info in accepted:
            plan.append((par, stratum, info))
    if verbose and short:
        print(f"   short of target: {short}", flush=True)

    # ---- Phase 2: benchmark in parallel ----
    if nproc is None:
        nproc = max(1, (mp.cpu_count() or 1) - 2)
    if verbose:
        print(f"[validation_matrix] phase 2: benchmarking {len(plan)} samples "
              f"on {nproc} workers, T={T}, rtol={rtol}, atol={atol}",
              flush=True)

    tasks = [
        (i, p.eps, p.gam, p.a, p.x, s, T, rtol, atol, info)
        for i, (p, s, info) in enumerate(plan)
    ]

    records: List[Record] = []
    failures: List[Tuple[int, str]] = []
    t0 = time.time()
    if nproc == 1:
        for ta in tasks:
            idx, rec, err = _bench_worker(ta)
            if err is None and rec is not None:
                records.append(rec)
            else:
                failures.append((idx, err or "no record"))
            done = len(records) + len(failures)
            if verbose:
                el = time.time() - t0
                rate = done / el
                eta = (len(tasks) - done) / max(rate, 1e-9)
                print(f"  [{done:4d}/{len(tasks)}]  "
                      f"stratum={plan[idx][1] if idx < len(plan) else '?'}  "
                      f"elapsed {el/60:5.1f}min  ETA {eta/60:5.1f}min",
                      flush=True)
    else:
        with mp.get_context("spawn").Pool(processes=nproc) as pool:
            done = 0
            last_ckpt = 0
            for idx, rec, err in pool.imap_unordered(_bench_worker, tasks):
                done += 1
                if err is None and rec is not None:
                    records.append(rec)
                else:
                    failures.append((idx, err or "no record"))
                if verbose and done % 5 == 0:
                    el = time.time() - t0
                    rate = done / el
                    eta = (len(tasks) - done) / max(rate, 1e-9)
                    print(f"  [{done:4d}/{len(tasks)}]  "
                          f"elapsed {el/60:5.1f}min  ETA {eta/60:5.1f}min",
                          flush=True)
                    if progress_path:
                        try:
                            with open(progress_path, "a") as f:
                                f.write(f"{time.time():.0f} {done}/{len(tasks)} "
                                        f"{el:.1f}s\n")
                        except Exception:
                            pass
                # checkpoint every 25 completions
                if done - last_ckpt >= 25:
                    last_ckpt = done
                    try:
                        save_dataset(records, out_path + ".tmp")
                    except Exception:
                        pass

    # group records by stratum for legibility
    records.sort(key=lambda r: (r.stratum, str(r.eps) + str(r.a)))
    save_dataset(records, out_path)
    if verbose:
        per = {s: 0 for s in STRATA}
        for r in records:
            per[r.stratum] += 1
        print(f"[validation_matrix] DONE in {(time.time()-t0)/60:.1f} min",
              flush=True)
        for s in STRATA:
            print(f"   {s}: {per[s]:4d} samples", flush=True)
        print(f"   {len(failures)} failures", flush=True)
        if failures[:3]:
            for f in failures[:3]:
                print(f"     failure example: {f}", flush=True)
        print(f"   dataset saved to {out_path}", flush=True)
    return records


# ===========================================================================
#  Persistence
# ===========================================================================
def save_dataset(records: List[Record], path: str = DATASET_PATH) -> None:
    with open(path, "wb") as f:
        pickle.dump(records, f, protocol=pickle.HIGHEST_PROTOCOL)


def load_dataset(path: str = DATASET_PATH) -> List[Record]:
    with open(path, "rb") as f:
        return pickle.load(f)


# ===========================================================================
#  Validation API
# ===========================================================================
def validate(formula: Callable[[Record], np.ndarray],
             dataset: Optional[List[Record]] = None
             ) -> Dict[str, np.ndarray]:
    """
    Apply ``formula(record) -> 3x3 P_pred`` to every record; group residuals
    by stratum.

    Returns
    -------
    {stratum -> 1D numpy array of max-norm residuals |P_pred - P_bench|}.
    """
    if dataset is None:
        dataset = load_dataset()
    out: Dict[str, List[float]] = {s: [] for s in STRATA}
    for r in dataset:
        try:
            P_pred = np.asarray(formula(r), dtype=float)
            if P_pred.shape != (3, 3):
                raise ValueError(f"formula returned shape {P_pred.shape}")
            res = float(np.max(np.abs(P_pred - r.P_bench)))
        except Exception as exc:
            # treat formula failure as 'infinity' residual
            res = float("inf")
            out.setdefault("_errors", []).append((r.stratum, repr(exc)))
        out[r.stratum].append(res)
    return {k: np.asarray(v) for k, v in out.items() if k != "_errors"}


def summarize(residual_dict: Dict[str, np.ndarray]) -> Dict[str, Dict[str, float]]:
    """
    Per-stratum residual statistics: count, mean, p50, p90, p99, max.
    """
    out: Dict[str, Dict[str, float]] = {}
    for s, arr in residual_dict.items():
        a = np.asarray(arr)
        if a.size == 0:
            out[s] = {"count": 0}
            continue
        out[s] = {
            "count": int(a.size),
            "mean": float(np.mean(a)),
            "p50": float(np.percentile(a, 50)),
            "p90": float(np.percentile(a, 90)),
            "p99": float(np.percentile(a, 99)),
            "max": float(np.max(a)),
        }
    return out


# ===========================================================================
#  Self-test formulas
# ===========================================================================
def identity_formula(rec: Record) -> np.ndarray:
    """Baseline: P_pred = I_3.  Tells you how far P_bench is from trivial."""
    return np.eye(3)


def be_only_formula(rec: Record) -> np.ndarray:
    """
    Pure Brundobler-Elser formula in the IP adiabatic endpoint basis.

    The two extreme-slope diabatic levels survive with probability
    $\prod_{m\neq n} q_{nm}$, $q_{ij}=\exp(-2\pi\Gamma_{ij})$.

    The IP harness uses the *exterior-first* spectral labelling: adiabatic
    channel 0 follows the spectral eigenvalue $\lambda_0$ exterior to the
    diagonal poles $\epsilon_j$, which at $u\to\pm\infty$ corresponds to the
    eigenvector with the *largest* energy.  Concretely
    \begin{align*}
       \text{adia 0 at }u\to-\infty &= \text{diabatic }\arg\min a,\\
       \text{adia 0 at }u\to+\infty &= \text{diabatic }\arg\max a,
    \end{align*}
    so the endpoint permutations are
    $\pi_{\mathrm{in}} = \arg\mathrm{sort}(a)$ and
    $\pi_{\mathrm{out}} = \arg\mathrm{sort}(-a)$.

    All non-BE entries are set to 0.
    """
    a = np.array(rec.a)
    N = 3
    pi_in = np.argsort(a)          # diabatic = pi_in[adiabatic incoming]
    pi_out = np.argsort(-a)        # diabatic = pi_out[adiabatic outgoing]
    n_max = int(np.argmax(a))       # steepest diabatic level
    n_min = int(np.argmin(a))       # shallowest diabatic level
    P = np.zeros((N, N))
    for n in (n_max, n_min):
        q_prod = 1.0
        for m in range(N):
            if m == n:
                continue
            pair = (n, m) if n < m else (m, n)
            q_prod *= rec.q_ij[pair]
        # adiabatic indices where diabatic n appears at each end
        i_in = int(np.where(pi_in == n)[0][0])
        i_out = int(np.where(pi_out == n)[0][0])
        P[i_in, i_out] = q_prod
    return P


def grid_formula(rec: Record) -> np.ndarray:
    """
    Incoherent LZ-grid formula (``assay.grid.grid_P``) reindexed into the
    IP adiabatic endpoint basis used by the benchmark.

    ``grid_P`` returns the transition probabilities in the *diabatic* basis.
    With the IP exterior-first spectral labelling
    (adia 0 = largest-$E$ eigenvector at each end),
    the endpoint permutations are
    $\pi_{\mathrm{in}} = \arg\mathrm{sort}(a)$ and
    $\pi_{\mathrm{out}} = \arg\mathrm{sort}(-a)$,
    so $P_{\mathrm{adia}}[i,j] = P_{\mathrm{dia}}[\pi_{\mathrm{in}}[i],\pi_{\mathrm{out}}[j]]$.
    """
    par = Params(eps=rec.eps, gam=rec.gam, a=rec.a, x=rec.x)
    geo = Geometry(par)
    P_dia = grid_P(geo)
    a = np.array(par.a)
    pi_in = np.argsort(a)
    pi_out = np.argsort(-a)
    return P_dia[np.ix_(pi_in, pi_out)]


# ===========================================================================
#  Reporting helpers
# ===========================================================================
def _print_summary_table(name: str, summary: Dict[str, Dict[str, float]]) -> None:
    print(f"\n--- {name} ---")
    hdr = f"  {'stratum':<8} {'n':>4} {'mean':>11} {'p50':>11} " \
          f"{'p90':>11} {'p99':>11} {'max':>11}"
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for s in STRATA:
        st = summary.get(s, {"count": 0})
        if st["count"] == 0:
            print(f"  {s:<8} {0:>4}   (no samples)")
            continue
        print(f"  {s:<8} {st['count']:>4} {st['mean']:>11.3e} {st['p50']:>11.3e} "
              f"{st['p90']:>11.3e} {st['p99']:>11.3e} {st['max']:>11.3e}")


def selftest_report(dataset: Optional[List[Record]] = None) -> Dict[str, Dict[str, Dict[str, float]]]:
    """
    Run the three baseline self-test formulas and print a report.
    """
    if dataset is None:
        dataset = load_dataset()
    out = {}
    for name, fml in [("Identity (P=I)", identity_formula),
                      ("BE-only (extreme survivals)", be_only_formula),
                      ("Incoherent LZ grid", grid_formula)]:
        res = validate(fml, dataset)
        smy = summarize(res)
        _print_summary_table(name, smy)
        out[name] = smy
    return out


# ===========================================================================
#  Richardson convergence demo (for the LaTeX note)
# ===========================================================================
def richardson_convergence_demo(seed: int = 20260523,
                                T_grid: Iterable[float] = (30, 60, 120),
                                T_ref: Optional[float] = None,
                                ) -> Dict[str, dict]:
    """
    On two reference samples, demonstrate convergence of P(T) and the
    Richardson combination as T grows.  Returns a dict suitable for the
    note's table.

    ``T_ref`` defaults to ``2 * max(T_grid)``.
    """
    rng = np.random.default_rng(seed)
    # two well-controlled samples
    samples = []
    # canonical generic sample
    samples.append(("canonical",
                    Params(eps=(-1.0, 0.0, 1.0), gam=(1.0, 1.0, 1.0),
                           a=(-1.0, 0.5, 2.0), x=0)))
    # a random S1
    for _ in range(50):
        try:
            par = _draw_one(rng, "S1")
            ok, _ = _accept_stratum(par, "S1")
            if ok:
                samples.append(("random S1", par))
                break
        except StratumRejection:
            continue

    if T_ref is None:
        T_ref = 2.0 * max(T_grid)

    out = {}
    for name, par in samples:
        geo = Geometry(par)
        rows = []
        U_ref = propagate_ad_ip(geo, -T_ref, T_ref, par.x, rtol=1e-13, atol=1e-14)
        P_ref = np.abs(U_ref.T) ** 2
        for T in T_grid:
            U_T = propagate_ad_ip(geo, -T, T, par.x, rtol=1e-13, atol=1e-14)
            P_T = np.abs(U_T.T) ** 2
            U_2T = propagate_ad_ip(geo, -2*T, 2*T, par.x, rtol=1e-13, atol=1e-14)
            P_2T = np.abs(U_2T.T) ** 2
            P_R = (8*P_2T - P_T) / 7.0
            rows.append({
                "T": float(T),
                "err_T": float(np.max(np.abs(P_T - P_ref))),
                "err_2T": float(np.max(np.abs(P_2T - P_ref))),
                "err_R": float(np.max(np.abs(P_R - P_ref))),
            })
        out[name] = {
            "eps": par.eps, "gam": par.gam, "a": par.a, "x": par.x,
            "T_ref": float(T_ref),
            "rows": rows,
        }
    return out


# ===========================================================================
#  CLI
# ===========================================================================
def _cli_build(args):
    build_dataset(samples_per_stratum=args.samples_per_stratum,
                  T=args.T, seed=args.seed, nproc=args.nproc,
                  out_path=args.out, verbose=True)


def _cli_selftest(args):
    selftest_report()


def _cli_richardson(args):
    rep = richardson_convergence_demo()
    for name, blk in rep.items():
        print(f"\n--- Richardson demo: {name} ---")
        print(f"  eps={blk['eps']}, gam={blk['gam']}, a={blk['a']}")
        print(f"  reference T_ref = {blk['T_ref']}")
        print(f"  {'T':>6} {'err P(T)':>13} {'err P(2T)':>13} {'err R':>13}")
        for r in blk["rows"]:
            print(f"  {r['T']:>6.0f} {r['err_T']:>13.3e} "
                  f"{r['err_2T']:>13.3e} {r['err_R']:>13.3e}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="cmd", required=True)

    pb = sub.add_parser("build", help="build the dataset")
    pb.add_argument("--samples-per-stratum", type=int,
                    default=DEFAULT_SAMPLES_PER_STRATUM)
    pb.add_argument("--T", type=float, default=120.0)
    pb.add_argument("--seed", type=int, default=20260523)
    pb.add_argument("--nproc", type=int, default=None)
    pb.add_argument("--out", type=str, default=DATASET_PATH)
    pb.set_defaults(func=_cli_build)

    ps = sub.add_parser("selftest", help="run baseline self-test formulas")
    ps.set_defaults(func=_cli_selftest)

    pr = sub.add_parser("richardson-demo",
                        help="show Richardson convergence on two samples")
    pr.set_defaults(func=_cli_richardson)

    a = p.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
