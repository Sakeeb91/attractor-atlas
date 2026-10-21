"""Full certification battery of the Trinacria attractor (D3, chiral dihedral).

Mirrors scripts/verify_mobula.py and is stricter where the D3 plan asks for it.
Every gate threshold below was fixed before the run; the pass/fail table at the
end of the JSON lists each criterion with its measured value.

  B1  Lyapunov spectrum: tangent-space QR, dt = 0.005, 20k transient steps,
      600k steps (T = 3000) from DEFAULT_STATE, 10-segment SEM.
      Gate: lambda1 >= 0.15, |lambda2| <= 0.01, lambda3 < 0, sum < 0.
  B2  dt-halving arm: lambda1 at dt in {0.01, 0.005, 0.0025, 0.00125} over the
      SAME physical window T = 2000 after a T = 100 transient, same IC,
      10-segment SEM each. Gate: successive relative change <= 10 percent,
      finest two within 5 percent, and no monotone drift across the four dt
      with |lambda1(0.01) - lambda1(0.00125)| above 2x the combined segment
      SEM. Spot check at dt = 0.000625 over T = 500 (reported, not gated).
  B3  Spectrum-divergence consistency: |sum lambda - mean div| <= 3 percent of
      |mean div|, mean divergence from the B5 orbit every 100 steps.
  B4  Kaplan-Yorke dimension from B1: 2 < D_KY < 3 (beauty band 2.05-2.5).
  B5  Boundedness: 4M steps at dt = 0.004 after 40k transient steps; finite,
      max|X| < 10, mean divergence < 0, |centroid_i| / extent_i <= 0.01 (the
      origin is the only D3-fixed point).
  B6  Initial-condition consistency: 5 ICs within +-0.3 of random points ON
      the attractor (the basin is not all of [-2, 2]^3, see B9), QR 150k steps
      at dt = 0.005. Gate: all lambda1 > 0.10, (max - min) / mean <= 25 percent.
  B7  Equilibria: analytic on-axis triple + 600-seed fsolve; classified;
      partitioned into D3-orbits. Gate: orbit sizes in {1, 2, 3, 6}, every
      orbit complete, sorted spectra agree within an orbit to 1e-6, every
      z = 0 equilibrium off the origin is fixed by a half-turn of D3, the
      on-axis pair has opposite rotation sense J_xy - J_yx.
  B8  Symmetry residuals: extent-normalized nearest-neighbour RMS (cKDTree,
      base = orbit[::66] ~ 60k points, 5k query points off the base) for the
      five non-identity D3 elements, the self floor (first half vs second
      half, matched 60k bases, both directions) and five controls that are
      not symmetries (z-flip, y-mirror, C6, inversion, 37-degree turn).
      Gate: max group <= 1.5 x floor and <= 0.015; min control >= 3x max group;
      C6 and y-mirror >= 5x.
  B9  Single D3-invariant attractor versus a kappa-conjugate pair:
      (a) basin protocol shared with the other specimens: 48 ICs uniform in
          [-2, 2]^3 (rng seed 0), dt = 0.005, escape |X| > 50, bounded count
          at T = 120; bounded orbits continue to T = 300 and land when their
          T in (200, 300] cloud has rho_id <= 0.02 against the B5 reference.
          Strict gate 48/48 (pre-registered); the gate accepted by the user
          for this specimen is "every bounded orbit lands on the symmetric
          set", with the escape fraction disclosed. A second draw of 480 ICs
          (rng seed 1, same protocol) gives the bounded fraction with a
          binomial standard error.
      (b) near-attractor basin: 48 ICs within +-0.3 of attractor points, same
          landing rule. Gate: 48/48 bounded and landed.
      (c) image launch from x0, kappa x0, rho x0 (1M steps each at dt = 0.004
          after 40k transient): early equivariance over the first 2000 steps
          <= 1e-8; cross-cloud residuals <= 1.5 x floor (gated); fingerprint
          distances against the 0.0549 floor (reported; the v1 gate in
          results/d3_selection.json explains why it is not gated).
      (d) z statistics of the B5 orbit: |P(z>0) - 0.5| <= 0.01, KS(z, -z) on
          10k decimated points <= 0.02, and the field-based angular velocity
          above / below the mid-plane equal and opposite within 5 percent of
          the mean speed (v1) or within 3 batch-means standard errors (v2).
  B10 Novelty: fingerprint(pts[::8]) against results/catalog_fingerprints.json
      (nearest + 3 runners-up) and against the four other atlas specimens.
      Interpretation only: distance must exceed 2 x 0.0549 = 0.110 to be
      interpretable and the median NN spacing 0.115 to be called distinct.

Writes results/trinacria_verification.json (partial JSON after each stage).

    python3 scripts/verify_trinacria.py            # full battery, ~15 min (4 cores)
    python3 scripts/verify_trinacria.py --quick    # smoke test, ~1 min
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import pathlib
import sys
import time

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import ks_2samp

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from atlas import family_d  # noqa: E402
from atlas.fingerprint import distance, fingerprint  # noqa: E402
from trinacria.dynamics import (  # noqa: E402
    ALPHA, BETA, DEFAULT_STATE, DELTA, EPS, GAMMA, LAM, NU, OMEGA,
    divergence, group_matrices, trajectory, velocity,
)
from trinacria.equilibria import (  # noqa: E402
    classify, find_equilibria, group_orbits, rotation_sense,
)
from trinacria.lyapunov import kaplan_yorke_dimension, lyapunov_spectrum  # noqa: E402

PARAMS = (ALPHA, BETA, OMEGA, GAMMA, DELTA, NU, LAM, EPS)
OUT = ROOT / "results" / "trinacria_verification.json"
FP_FLOOR = 0.0549          # re-simulation noise floor of the fingerprint distance
NN_SPACING = 0.1149        # median nearest-neighbour spacing inside the catalog
ESCAPE = 50.0

ELEMENT_NAMES = ["identity", "rho", "rho2", "kappa", "rho_kappa", "rho2_kappa"]


def _rotz(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


GROUP = dict(zip(ELEMENT_NAMES, group_matrices()))
KAPPA = GROUP["kappa"]
RHO = GROUP["rho"]
CONTROLS = {
    "z_flip": np.diag([1.0, 1.0, -1.0]),
    "y_mirror": np.diag([1.0, -1.0, 1.0]),
    "C6": _rotz(np.pi / 3),
    "inversion": -np.eye(3),
    "rot37": _rotz(np.deg2rad(37.0)),
}


def r(v, k=5):
    return None if v is None else round(float(v), k)


def rlist(v, k=5):
    return [r(x, k) for x in np.asarray(v).ravel()]


def write(res):
    OUT.write_text(json.dumps(res, indent=2) + "\n")


# ---------------------------------------------------------------------------
# Batched integration (real form, vectorized over rows) for launches.
# ---------------------------------------------------------------------------

def vel_batch(S):
    x, y, z = S[:, 0], S[:, 1], S[:, 2]
    c = ALPHA * z * z - BETA
    d2 = x * x - y * y
    return np.stack([
        c * x - OMEGA * z * y + GAMMA * d2 + 2 * DELTA * z * x * y,
        c * y + OMEGA * z * x - 2 * GAMMA * x * y + DELTA * z * d2,
        NU * z - z ** 3 - LAM * z * (x * x + y * y) + EPS * (3 * x * x * y - y ** 3),
    ], axis=1)


def rk4_batch(S, dt):
    k1 = vel_batch(S)
    k2 = vel_batch(S + 0.5 * dt * k1)
    k3 = vel_batch(S + 0.5 * dt * k2)
    k4 = vel_batch(S + dt * k3)
    return S + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)


def complex_basin_orbits(S0, dt, T_bound, T_total, T_sample_from, every=10):
    """Batched complex RK4 of atlas.family_d (n = 3) from rows of S0.

    Returns (bounded_at_T_bound, bounded_at_T_total, samples (B, n, 3)).
    Escaped rows (non-finite or |X| > ESCAPE, checked every 10 steps) are
    frozen at NaN.
    """
    S0 = np.asarray(S0, float)
    B = len(S0)
    P = np.tile(np.asarray(PARAMS, float), (B, 1))
    W = S0[:, 0] + 1j * S0[:, 1]
    Z = S0[:, 2].astype(float).copy()
    alive = np.ones(B, bool)
    n_bound = int(round(T_bound / dt))
    n_total = int(round(T_total / dt))
    n_from = int(round(T_sample_from / dt))
    bounded_T = None
    pts = []
    with np.errstate(all="ignore"):
        for i in range(n_total):
            W, Z = family_d._rk4(W, Z, P, 3, dt)
            if i % 10 == 0 or i == n_total - 1:
                bad = (~np.isfinite(W)) | (~np.isfinite(Z)) | \
                      (np.maximum(np.abs(W), np.abs(Z)) > ESCAPE)
                if bad.any():
                    alive &= ~bad
                    W[~alive] = 0.0
                    Z[~alive] = 0.0
            if i + 1 == n_bound:
                bounded_T = alive.copy()
            if i >= n_from and (i - n_from) % every == 0:
                X = np.stack([W.real, W.imag, Z], 1)
                X[~alive] = np.nan
                pts.append(X)
    pts = np.array(pts).transpose(1, 0, 2)
    pts[~alive] = np.nan
    return bounded_T, alive.copy(), pts


# ---------------------------------------------------------------------------
# Worker jobs (pure-Python tangent QR, one process each).
# ---------------------------------------------------------------------------

def _qr_job(kw):
    t0 = time.time()
    spec, seg = lyapunov_spectrum(params=PARAMS, segments=10, **kw)
    return {"spectrum": spec.tolist(), "segments": seg.tolist(),
            "seconds": time.time() - t0}


def _atlas_fp_job(name):
    import importlib
    mod = importlib.import_module(f"{name}.dynamics")
    pts = mod.trajectory()
    return name, fingerprint(pts[::8]).tolist(), float(np.ptp(pts, axis=0).max())


def seg_sem(seg):
    seg = np.asarray(seg)
    return float(seg[:, 0].std(ddof=1) / np.sqrt(len(seg)))


# ---------------------------------------------------------------------------
# Symmetry residuals (B8).
# ---------------------------------------------------------------------------

def residual_battery(pts, base_stride=66, n_query=5000, seed=0):
    """Extent-normalized NN-RMS residuals. The cloud is NOT re-centred: the
    origin is the only D3-fixed point, so an off-centre cloud is penalised."""
    q = pts[np.all(np.isfinite(pts), axis=1)]
    N = len(q)
    extent = float(np.ptp(q, axis=0).max())
    rng = np.random.default_rng(seed)
    base = q[::base_stride]
    tree = cKDTree(base)
    off = np.flatnonzero(np.arange(N) % base_stride != 0)
    qi = q[rng.choice(off, size=min(n_query, len(off)), replace=False)]

    def rms(tree_, X):
        d, _ = tree_.query(X)
        return float(np.sqrt(np.mean(d ** 2)) / extent)

    half = N // 2
    A, Bh = q[:half], q[half:]
    sa = max(1, base_stride // 2)  # matched base size (~60k) inside each half
    ta, tb = cKDTree(A[::sa]), cKDTree(Bh[::sa])
    qa = A[rng.choice(len(A), size=n_query, replace=False)]
    qb = Bh[rng.choice(len(Bh), size=n_query, replace=False)]
    f_ab, f_ba = rms(tb, qa), rms(ta, qb)
    floor = 0.5 * (f_ab + f_ba)

    grp, dets = {}, {}
    for k, m in GROUP.items():
        dets[k] = float(np.linalg.det(m))
        assert abs(dets[k] - 1.0) < 1e-12
        if k != "identity":
            grp[k] = rms(tree, qi @ m.T)
    ident = rms(tree, qi)
    ctl = {k: rms(tree, qi @ m.T) for k, m in CONTROLS.items()}
    gmax = max(grp.values())
    out = {
        "n_points": N, "n_base": len(base), "n_query": len(qi), "extent": r(extent, 4),
        "centering": "origin (the only D3-fixed point); no re-centring",
        "determinants": {k: r(v, 12) for k, v in dets.items()},
        "group": {k: r(v) for k, v in grp.items()},
        "identity_same_base": r(ident),
        "floor_half_vs_half": r(floor), "floor_directions": [r(f_ab), r(f_ba)],
        "controls": {k: r(v) for k, v in ctl.items()},
        "max_group": r(gmax), "min_control": r(min(ctl.values())),
        "group_over_floor": r(gmax / floor, 3),
        "min_ratio": r(min(ctl.values()) / gmax, 3),
        "control_ratios": {k: r(v / gmax, 3) for k, v in ctl.items()},
    }
    out["pass_group_vs_floor"] = bool(gmax <= 1.5 * floor and gmax <= 0.015)
    out["pass_controls_3x"] = bool(min(ctl.values()) >= 3 * gmax)
    out["pass_c6_ymirror_5x"] = bool(ctl["C6"] >= 5 * gmax and ctl["y_mirror"] >= 5 * gmax)
    out["pass"] = out["pass_group_vs_floor"] and out["pass_controls_3x"] and out["pass_c6_ymirror_5x"]
    return out, base, extent, floor


def landing(samples, ref_tree, extent, thr=0.02):
    rid, rk = [], []
    for c in samples:
        if not np.all(np.isfinite(c)):
            rid.append(None); rk.append(None)
            continue
        d, _ = ref_tree.query(c)
        rid.append(r(np.sqrt(np.mean(d ** 2)) / extent))
        d, _ = ref_tree.query(c @ KAPPA.T)
        rk.append(r(np.sqrt(np.mean(d ** 2)) / extent))
    landed = [v is not None and v <= thr for v in rid]
    return rid, rk, landed


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    Q = args.quick
    t_start = time.time()
    sc = 0.02 if Q else 1.0  # step-count scale for the smoke test

    def n(k):
        return max(2000, int(round(k * sc)))

    res = {
        "system": "Trinacria attractor",
        "symmetry": "D3 chiral dihedral (Schoenflies), order 6, generators "
                    "rho:(w,z)->(e^{2 pi i/3} w, z) and kappa:(w,z)->(conj(w),-z) "
                    "= (x,y,z)->(x,-y,-z), a half-turn about x",
        "parameters": {"alpha": ALPHA, "beta": BETA, "omega": OMEGA, "gamma": GAMMA,
                       "delta": DELTA, "nu": NU, "lam": LAM, "eps": EPS},
        "default_state": list(DEFAULT_STATE),
        "quick": Q,
        "selection": "nominee #8 of results/d3_candidates.json, accepted by the user "
                     "with the relaxed basin gate 'every bounded orbit lands on the "
                     "symmetric set' (results/d3_selection.json)",
    }

    # Batched real form must equal the package velocity.
    rng = np.random.default_rng(3)
    X = rng.uniform(-1.5, 1.5, (300, 3))
    err = float(np.abs(vel_batch(X) - np.array([velocity(x) for x in X])).max())
    assert err < 1e-12, err
    W = X[:, 0] + 1j * X[:, 1]
    dW, dZ = family_d._deriv(W, X[:, 2], np.tile(np.asarray(PARAMS), (300, 1)), 3)
    err_c = float(np.abs(np.stack([dW.real, dW.imag, dZ], 1) - vel_batch(X)).max())
    assert err_c < 1e-12, err_c
    res["integrator_crosscheck"] = {"real_batch_vs_package": err, "complex_family_vs_real": err_c}

    # ------------------------------------------------------------------
    # Launch the QR jobs and the atlas fingerprints in a worker pool while
    # the main process integrates the long orbit.
    # ------------------------------------------------------------------
    print("short orbit for attractor-point ICs...", flush=True)
    pre = trajectory(n_steps=n(50_000), dt=0.004, transient=n(40_000))
    rng_ic = np.random.default_rng(11)
    ic_pts = pre[rng_ic.choice(len(pre), 5, replace=False)] + rng_ic.uniform(-0.3, 0.3, (5, 3))

    DTS = [0.01, 0.005, 0.0025, 0.00125]
    jobs = {"main": dict(n_steps=n(600_000), dt=0.005, transient=n(20_000),
                         state0=DEFAULT_STATE)}
    for dt in DTS:
        jobs[f"dt_{dt}"] = dict(n_steps=n(int(round(2000 / dt))), dt=dt,
                                transient=n(int(round(100 / dt))), state0=DEFAULT_STATE)
    jobs["dt_0.000625"] = dict(n_steps=n(800_000), dt=0.000625,
                               transient=n(160_000), state0=DEFAULT_STATE)
    for i, s0 in enumerate(ic_pts):
        jobs[f"ic_{i}"] = dict(n_steps=n(150_000), dt=0.005, transient=n(20_000),
                               state0=tuple(s0))
    order = sorted(jobs, key=lambda k: -(jobs[k]["n_steps"] + jobs[k]["transient"]))
    pool = mp.get_context("fork").Pool(processes=3)
    handles = {k: pool.apply_async(_qr_job, (jobs[k],)) for k in order}
    atlas_handles = [pool.apply_async(_atlas_fp_job, (nm,))
                     for nm in ("aurelia", "naiad", "cassiopea", "mobula")]

    # ------------------------------------------------------------------
    # B5 boundedness on the long orbit (main process).
    # ------------------------------------------------------------------
    print("B5 long orbit (4M steps)...", flush=True)
    t0 = time.time()
    pts = trajectory(n_steps=n(4_000_000), dt=0.004, transient=n(40_000))
    finite = bool(np.all(np.isfinite(pts)))
    maxabs = float(np.abs(pts).max())
    divs = np.array([divergence(p) for p in pts[::100]])
    mean_div = float(divs.mean())
    nb = 20
    L = len(divs) // nb
    div_se = float(np.array([divs[k * L:(k + 1) * L].mean() for k in range(nb)]).std(ddof=1) / np.sqrt(nb))
    cen = pts.mean(0)
    ext = np.ptp(pts, axis=0)
    coe = np.abs(cen) / ext
    ev = np.sort(np.linalg.eigvalsh(np.cov(pts.T)))[::-1]
    res["bounded_4M_steps"] = bool(finite and maxabs < 10.0)
    res["long_orbit"] = {
        "n_steps": len(pts), "dt": 0.004, "transient_steps": n(40_000),
        "finite": finite, "max_abs": r(maxabs, 4),
        "centroid": rlist(cen, 5), "extent": rlist(ext, 4),
        "centroid_over_extent": rlist(coe, 5),
        "iso_min_over_max_span": r(ext.min() / ext.max(), 4),
        "pca_e3_over_e1": r(ev[2] / ev[0], 4),
        "seconds": r(time.time() - t0, 1),
    }
    res["extent_min"] = rlist(pts.min(0), 4)
    res["extent_max"] = rlist(pts.max(0), 4)
    res["mean_divergence"] = r(mean_div, 5)
    res["mean_divergence_batch_se"] = r(div_se, 5)
    res["centroid"] = {"mean": rlist(cen, 5), "over_extent": rlist(coe, 5),
                       "threshold": 0.01, "pass": bool(np.all(coe <= 0.01))}
    print(f"  bounded={res['bounded_4M_steps']} max|X|={maxabs:.3f} mean div={mean_div:.4f} "
          f"centroid/extent={np.round(coe, 4)}", flush=True)
    write(res)

    # ------------------------------------------------------------------
    # B8 symmetry residuals on the long orbit.
    # ------------------------------------------------------------------
    print("B8 symmetry residuals...", flush=True)
    stride = 66 if not Q else 2
    sym, ref_base, extent, floor = residual_battery(pts, base_stride=stride,
                                                    n_query=5000 if not Q else 1000)
    res["symmetry_residuals"] = sym
    print(f"  group={sym['group']} floor={sym['floor_half_vs_half']} "
          f"controls={sym['controls']} min ratio={sym['min_ratio']}", flush=True)
    write(res)
    ref_tree = cKDTree(ref_base)

    # ------------------------------------------------------------------
    # B9 single attractor tests.
    # ------------------------------------------------------------------
    print("B9(d) z statistics and counter-rotation...", flush=True)
    z = pts[:, 2]
    zz = z[:: max(1, len(z) // 10_000)]
    ks = float(ks_2samp(zz, -zz).statistic)
    th = np.unwrap(np.arctan2(pts[::10, 1], pts[::10, 0]))
    dth = np.diff(th) / (10 * 0.004)
    zm = 0.5 * (z[::10][1:] + z[::10][:-1])
    om_up_v1, om_dn_v1 = float(dth[zm > 0].mean()), float(dth[zm < 0].mean())
    sub = pts[::4]
    V = vel_batch(sub)
    rr2 = sub[:, 0] ** 2 + sub[:, 1] ** 2
    omf = (sub[:, 0] * V[:, 1] - sub[:, 1] * V[:, 0]) / rr2
    lz = sub[:, 0] * V[:, 1] - sub[:, 1] * V[:, 0]
    up = sub[:, 2] > 0
    Lb = len(sub) // nb

    def bm(vals, mask=None):
        v = []
        for k in range(nb):
            sl = slice(k * Lb, (k + 1) * Lb)
            x = vals[sl] if mask is None else vals[sl][mask[sl]]
            v.append(x.mean())
        v = np.array(v)
        return float(v.mean()), float(v.std(ddof=1) / np.sqrt(nb))

    pu = bm(up.astype(float))
    wu, wd = bm(omf, up), bm(omf, ~up)
    lu, ld = bm(lz, up), bm(lz, ~up)
    se_w = float(np.hypot(wu[1], wd[1]))
    mismatch = abs(wu[0] + wd[0])
    mean_speed = 0.5 * (abs(wu[0]) + abs(wd[0]))
    sgn = np.sign(z)
    flips = np.flatnonzero(np.diff(sgn) != 0)
    runs = np.diff(np.concatenate([[0], flips + 1, [len(z)]]))
    p_pos = float(np.mean(z > 0))
    zs = {
        "p_positive": r(p_pos, 5), "p_positive_batch_se": r(pu[1], 5),
        "ks_z_vs_minus_z": r(ks, 5), "ks_n_points": len(zz),
        "omega_above": [r(wu[0]), r(wu[1])], "omega_below": [r(wd[0]), r(wd[1])],
        "omega_estimator": "field-based (x ydot - y xdot)/(x^2+y^2), batch means (20 batches), [mean, se]",
        "omega_mismatch_over_mean_speed": r(mismatch / mean_speed, 4),
        "omega_mismatch_over_se": r(mismatch / se_w, 3),
        "Lz_above": [r(lu[0]), r(lu[1])], "Lz_below": [r(ld[0]), r(ld[1])],
        "dtheta_dt_above_v1": r(om_up_v1), "dtheta_dt_below_v1": r(om_dn_v1),
        "longest_dwell_fraction": r(runs.max() / len(z), 5),
        "longest_dwell_time": r(runs.max() * 0.004, 2),
        "n_sign_changes": int(len(flips)),
        "pass_p_positive": bool(abs(p_pos - 0.5) <= 0.01),
        "pass_ks": bool(ks <= 0.02),
        "pass_counter_rotation_v1": bool(wu[0] * wd[0] < 0 and mismatch <= 0.05 * mean_speed),
        "pass_counter_rotation": bool(wu[0] * wd[0] < 0 and (mismatch <= 0.05 * mean_speed
                                                            or mismatch <= 3 * se_w)),
    }
    res["single_attractor_test"] = {"z_symmetry": zs}
    print(f"  P(z>0)={p_pos:.4f} KS={ks:.4f} omega above {wu} below {wd}", flush=True)
    write(res)

    print("B9(a) basin protocol (48 ICs in [-2,2]^3)...", flush=True)
    rng0 = np.random.default_rng(0)
    S0 = rng0.uniform(-2, 2, (48, 3))
    Tq = 0.1 if Q else 1.0
    bT, bEnd, samp = complex_basin_orbits(S0, 0.005, 120.0 * Tq, 300.0 * Tq, 200.0 * Tq)
    rid, rk, landed = landing(samp, ref_tree, extent)
    nb_120, nb_end = int(bT.sum()), int(bEnd.sum())
    n_land = int(sum(landed))
    res["single_attractor_test"]["basin"] = {
        "protocol": "48 ICs uniform in [-2,2]^3, rng seed 0, dt 0.005, batched complex RK4 "
                    "(atlas.family_d), escape |X|>50; bounded count at T=120; landing on "
                    "T in (200,300], every 10th step, rho_id <= 0.02 vs the B5 reference",
        "n_ic": 48, "n_bounded_T120": nb_120, "n_bounded_T300": nb_end,
        "n_escaped_T300": 48 - nb_end, "n_landed": n_land,
        "bounded_fraction": r(nb_end / 48, 4),
        "landed_fraction_of_bounded": r(n_land / nb_end, 4) if nb_end else None,
        "rho_id": rid, "rho_kappa": rk, "threshold": 0.02,
        "comparison_same_protocol": {"Aurelia": "48/48", "Naiad": "48/48",
                                     "Cassiopea": "43/48", "Mobula": "41/48"},
        "pass_strict_48_of_48": bool(n_land == 48),
        "pass_multistable_90pct": bool(n_land >= 0.9 * 48),
        "pass_relaxed_every_bounded_lands": bool(nb_end > 0 and n_land == nb_end),
    }
    print(f"  bounded T120 {nb_120}/48, T300 {nb_end}/48, landed {n_land}", flush=True)

    # The 48-IC count is a small sample; a 480-IC draw (rng seed 1, same
    # protocol) pins the bounded fraction down with a binomial standard error.
    n_big = 96 if Q else 480
    S2 = np.random.default_rng(1).uniform(-2, 2, (n_big, 3))
    bT2, bEnd2, samp2 = complex_basin_orbits(S2, 0.005, 120.0 * Tq, 300.0 * Tq, 200.0 * Tq)
    rid2, _, landed2 = landing(samp2, ref_tree, extent)
    nb2, nl2 = int(bEnd2.sum()), int(sum(landed2))
    f2 = nb2 / n_big
    worst = max((v for v in rid2 if v is not None), default=None)
    res["single_attractor_test"]["basin_large_sample"] = {
        "protocol": f"{n_big} ICs uniform in [-2,2]^3, rng seed 1, otherwise identical to the basin protocol",
        "n_ic": n_big, "n_bounded_T120": int(bT2.sum()), "n_bounded_T300": nb2, "n_landed": nl2,
        "bounded_fraction": r(f2, 4), "bounded_fraction_binomial_se": r(np.sqrt(f2 * (1 - f2) / n_big), 4),
        "max_rho_id_of_bounded": worst,
        "every_bounded_lands": bool(nb2 > 0 and nl2 == nb2),
    }
    res["single_attractor_test"]["basin"]["reconciliation"] = (
        "the selection stage (results/d3_selection.json) reported 29/48 bounded from a different "
        "IC draw (rng seed 7 after reference subsampling, dt 0.004, T = 300); the shared "
        "comparison protocol used for Aurelia/Naiad/Cassiopea/Mobula (rng seed 0, dt 0.005, "
        "T = 120) gives the n_bounded_T120 count recorded here, and the 480-IC sample in "
        "basin_large_sample gives the bounded fraction with its standard error")
    print(f"  large sample: bounded {nb2}/{n_big} ({f2:.3f}), landed {nl2}, max rho_id {worst}", flush=True)

    print("B9(b) near-attractor basin...", flush=True)
    rng1 = np.random.default_rng(5)
    S1 = pts[rng1.choice(len(pts), 48, replace=False)] + rng1.uniform(-0.3, 0.3, (48, 3))
    bT1, bEnd1, samp1 = complex_basin_orbits(S1, 0.005, 120.0 * Tq, 300.0 * Tq, 200.0 * Tq)
    rid1, rk1, landed1 = landing(samp1, ref_tree, extent)
    n1 = int(sum(landed1))
    res["single_attractor_test"]["basin_near_attractor"] = {
        "protocol": "48 ICs = random B5 orbit points + uniform(-0.3,0.3)^3, rng seed 5; same "
                    "integration and landing rule as the basin protocol",
        "n_ic": 48, "n_bounded_T300": int(bEnd1.sum()), "n_landed": n1,
        "rho_id": rid1, "rho_kappa": rk1,
        "pass": bool(n1 == 48),
    }
    print(f"  near-attractor bounded {int(bEnd1.sum())}/48 landed {n1}", flush=True)
    write(res)

    print("B9(c) image launch x0, kappa x0, rho x0 (3 x 1M steps)...", flush=True)
    x0 = np.asarray(DEFAULT_STATE, float)
    S = np.stack([x0, KAPPA @ x0, RHO @ x0])
    dt = 0.004
    early = []
    n_early = 2000
    n_tr, n_run = n(40_000), n(1_000_000)
    keep = []
    for i in range(n_tr + n_run):
        if i < n_early:
            early.append(S.copy())
        S = rk4_batch(S, dt)
        if i >= n_tr and (i - n_tr) % 8 == 0:
            keep.append(S.copy())
    early = np.array(early)
    eq_k = float(np.abs(early[:, 1] - early[:, 0] @ KAPPA.T).max())
    eq_r = float(np.abs(early[:, 2] - early[:, 0] @ RHO.T).max())
    clouds = np.array(keep).transpose(1, 0, 2)
    alive3 = [bool(np.all(np.isfinite(c)) and np.abs(c).max() < 10) for c in clouds]
    names = ["x0", "kappa_x0", "rho_x0"]
    fps = [fingerprint(c) for c in clouds]
    rng2 = np.random.default_rng(9)
    fpd, crs = {}, {}
    nbase = len(ref_base)
    for i in range(3):
        ti = cKDTree(clouds[i][rng2.choice(len(clouds[i]), min(nbase, len(clouds[i])), replace=False)])
        for j in range(3):
            if i == j:
                continue
            qj = clouds[j][rng2.choice(len(clouds[j]), min(5000, len(clouds[j])), replace=False)]
            d, _ = ti.query(qj)
            crs[f"{names[j]}->{names[i]}"] = r(np.sqrt(np.mean(d ** 2)) / extent)
            if i < j:
                fpd[f"{names[i]}|{names[j]}"] = r(distance(fps[i], fps[j]), 4)
    half = len(clouds[0]) // 2
    self_fp = r(distance(fingerprint(clouds[0][:half]), fingerprint(clouds[0][half:])), 4)
    cr_lim = 1.5 * floor
    res["single_attractor_test"]["image_launch"] = {
        "x0": list(DEFAULT_STATE), "dt": dt, "transient_steps": n_tr, "n_steps": n_run,
        "sample_every": 8, "bounded": alive3,
        "early_equivariance_kappa": eq_k, "early_equivariance_rho": eq_r,
        "early_window_steps": n_early,
        "cloud_residuals": crs, "cloud_residual_limit": r(cr_lim),
        "fingerprint_distances": fpd, "fingerprint_floor": FP_FLOOR,
        "fingerprint_self_halves_x0": self_fp,
        "pass_equivariance": bool(eq_k <= 1e-8 and eq_r <= 1e-8),
        "pass_cloud_residuals": bool(max(crs.values()) <= cr_lim),
        "pass_fingerprint_v1": bool(all(v < FP_FLOOR for v in fpd.values())),
    }
    res["single_attractor_test"]["image_launch"]["pass"] = bool(
        all(alive3) and res["single_attractor_test"]["image_launch"]["pass_equivariance"]
        and res["single_attractor_test"]["image_launch"]["pass_cloud_residuals"])
    print(f"  eq kappa {eq_k:.2e} rho {eq_r:.2e} cloud res {crs} fp {fpd} self {self_fp}", flush=True)
    write(res)

    # ------------------------------------------------------------------
    # B7 equilibria.
    # ------------------------------------------------------------------
    print("B7 equilibria...", flush=True)
    eqs = find_equilibria(PARAMS, n_seeds=600 if not Q else 100)
    res["equilibria"] = []
    for e in eqs:
        eigs, label = classify(e, PARAMS)
        res["equilibria"].append({
            "point": [round(float(v), 6) for v in e],
            "radius_xy": r(np.hypot(e[0], e[1]), 6),
            "angle_deg": r(np.degrees(np.arctan2(e[1], e[0])), 4) if np.hypot(e[0], e[1]) > 1e-9 else None,
            "eigenvalues": [[round(float(v.real), 5), round(float(v.imag), 5)] for v in eigs],
            "type": label,
            "rotation_sense_Jxy_minus_Jyx": r(rotation_sense(e, PARAMS), 6),
            "residual_norm": float(np.linalg.norm(velocity(e))),
        })
    orbs = group_orbits(eqs, PARAMS)
    orbit_out = []
    for o in orbs:
        rep = eqs[o["members"][0]]
        maps = {}
        for j in o["members"]:
            maps[str(j)] = [k for k, g in GROUP.items() if np.linalg.norm(g @ rep - eqs[j]) < 1e-6]
        orbit_out.append({
            "members": o["members"], "size": o["size"], "complete": o["complete"],
            "label": o["label"],
            "representative": [round(float(v), 6) for v in rep],
            "spectrum": [[round(float(v.real), 5), round(float(v.imag), 5)] for v in o["spectrum"]],
            "spectrum_gap": o["spectrum_gap"],
            "elements_mapping_representative_to_member": maps,
        })
    on_axis = [e for e in eqs if np.hypot(e[0], e[1]) < 1e-9]
    triple = (len(on_axis) == 3 and any(np.linalg.norm(e) < 1e-12 for e in on_axis)
              and any(abs(e[2] - np.sqrt(NU)) < 1e-12 for e in on_axis)
              and any(abs(e[2] + np.sqrt(NU)) < 1e-12 for e in on_axis))
    pair = [e for e in on_axis if abs(e[2]) > 1e-9]
    sense = [rotation_sense(e, PARAMS) for e in pair]
    spec_pair_gap = float(np.max(np.abs(np.sort_complex(classify(pair[0])[0])
                                        - np.sort_complex(classify(pair[1])[0])))) if len(pair) == 2 else None
    half_turns = [GROUP[k] for k in ("kappa", "rho_kappa", "rho2_kappa")]
    plane = [e for e in eqs if abs(e[2]) < 1e-9 and np.hypot(e[0], e[1]) > 1e-9]
    plane_fixed = [min(float(np.linalg.norm(h @ e - e)) for h in half_turns) for e in plane]
    origin_eigs = np.sort(np.real(classify(np.zeros(3))[0]))
    eqchk = {
        "n_equilibria": len(eqs),
        "on_axis_triple_present": bool(triple),
        "orbit_sizes": [o["size"] for o in orbs],
        "all_orbits_complete": bool(all(o["complete"] for o in orbs)),
        "orbit_sizes_valid": bool(all(o["size"] in (1, 2, 3, 6) for o in orbs)),
        "max_spectrum_gap_within_orbit": float(max(o["spectrum_gap"] for o in orbs)),
        "n_plane_equilibria_off_origin": len(plane),
        "plane_equilibria_angles_deg": [r(np.degrees(np.arctan2(e[1], e[0])), 6) for e in plane],
        "plane_equilibria_max_half_turn_defect": max(plane_fixed) if plane_fixed else None,
        "on_axis_pair_rotation_sense": [r(s, 6) for s in sense],
        "on_axis_pair_spectrum_gap": spec_pair_gap,
        "origin_eigenvalues": rlist(origin_eigs, 6),
        "origin_note": "saddle with the double real eigenvalue -beta (rotation rate omega*z vanishes "
                       "at z=0) and nu along z; not a saddle-focus, so the homoclinic story of the "
                       "C_n specimens does not transfer",
    }
    eqchk["pass"] = bool(triple and eqchk["all_orbits_complete"] and eqchk["orbit_sizes_valid"]
                         and eqchk["max_spectrum_gap_within_orbit"] <= 1e-6
                         and (not plane_fixed or max(plane_fixed) <= 1e-6)
                         and len(sense) == 2 and sense[0] * sense[1] < 0
                         and spec_pair_gap is not None and spec_pair_gap <= 1e-6)
    res["equilibrium_orbits"] = orbit_out
    res["equilibrium_checks"] = eqchk
    print(f"  {len(eqs)} equilibria, orbit sizes {eqchk['orbit_sizes']}, gap "
          f"{eqchk['max_spectrum_gap_within_orbit']:.1e}, pass {eqchk['pass']}", flush=True)
    write(res)

    # ------------------------------------------------------------------
    # B10 novelty.
    # ------------------------------------------------------------------
    print("B10 novelty...", flush=True)
    cat = json.loads((ROOT / "results" / "catalog_fingerprints.json").read_text())
    fp = fingerprint(pts[::8])
    dists = sorted((float(distance(fp, np.array(v["fingerprint"]))), k) for k, v in cat.items())
    fp500 = fingerprint(pts[: len(pts) // 8: 8])
    nov = {
        "catalog_size": len(cat),
        "nearest": dists[0][1], "distance": r(dists[0][0], 4),
        "next": [[k, r(d, 4)] for d, k in dists[1:4]],
        "noise_floor": FP_FLOOR, "interpretable_above": round(2 * FP_FLOOR, 4),
        "median_nn_spacing": NN_SPACING,
        "interpretable": bool(dists[0][0] > 2 * FP_FLOOR),
        "distinct": bool(dists[0][0] > NN_SPACING),
        "nn_rank_unstable": bool(dists[1][0] - dists[0][0] < FP_FLOOR),
        "self_sampling_check_first_eighth_vs_full": r(distance(fp, fp500), 4),
        "hedge": "catalog-relative shape distance against 135 dysts systems, not a literature claim",
    }
    res["novelty_vs_dysts"] = nov
    print(f"  nearest {nov['nearest']} {nov['distance']} next {nov['next']}", flush=True)

    # ------------------------------------------------------------------
    # Collect the worker results.
    # ------------------------------------------------------------------
    print("waiting for QR jobs...", flush=True)
    qr = {k: h.get() for k, h in handles.items()}
    atl = {}
    for h in atlas_handles:
        nm, f, _ = h.get()
        atl[nm] = f
    pool.close(); pool.join()

    va = {nm.capitalize(): r(distance(fp, np.array(f)), 4) for nm, f in atl.items()}
    res["novelty_vs_atlas"] = {
        "method": "fingerprint of each specimen's default trajectory() (500k steps at dt 0.004) "
                  "[::8] vs the Trinacria 4M-step orbit [::8]",
        "distances": va,
        "nearest": min(va, key=va.get), "nearest_distance": min(va.values()),
        "aurelia_above_0.110": bool(va["Aurelia"] > 2 * FP_FLOOR),
    }
    print(f"  atlas distances {va}", flush=True)

    m = qr["main"]
    spec = np.array(m["spectrum"])
    ky = kaplan_yorke_dimension(spec)
    res["lyapunov_spectrum"] = rlist(spec, 5)
    res["lyapunov_segment_sem_lambda1"] = r(seg_sem(m["segments"]), 5)
    res["lyapunov_settings"] = {"dt": 0.005, "n_steps": jobs["main"]["n_steps"],
                                "transient_steps": jobs["main"]["transient"],
                                "state0": list(DEFAULT_STATE), "segments": 10}
    res["lyapunov_sum"] = r(spec.sum(), 5)
    res["kaplan_yorke_dimension"] = r(ky, 5)
    print(f"  spectrum {np.round(spec, 5)} D_KY {ky:.4f}", flush=True)

    conv = {}
    for dt in DTS + [0.000625]:
        j = qr[f"dt_{dt}"]
        sp = j["spectrum"]
        conv[str(dt)] = [r(sp[0]), r(sp[1]), r(sp[2]), r(seg_sem(j["segments"]))]
    l1 = np.array([conv[str(dt)][0] for dt in DTS])
    sem = np.array([conv[str(dt)][3] for dt in DTS])
    succ = [abs(l1[i + 1] - l1[i]) / abs(l1[i + 1]) for i in range(len(DTS) - 1)]
    finest = abs(l1[-1] - l1[-2]) / abs(l1[-1])
    d = np.diff(l1)
    monotone = bool(np.all(d > 0) or np.all(d < 0))
    drift = abs(l1[0] - l1[-1])
    drift_lim = 2 * float(np.hypot(sem[0], sem[-1]))
    spot = conv["0.000625"]
    res["dt_convergence"] = conv
    res["dt_convergence_summary"] = {
        "format": "{dt: [lambda1, lambda2, lambda3, segment_sem_lambda1]}; T = 2000 after T = 100 "
                  "transient from DEFAULT_STATE; spot check dt 0.000625 over T = 500",
        "successive_relative_change": [r(v, 4) for v in succ],
        "finest_two_relative_change": r(finest, 4),
        "monotone": monotone, "drift_coarsest_to_finest": r(drift),
        "drift_limit_2x_combined_sem": r(drift_lim),
        "spot_check_minus_finest": r(spot[0] - l1[-1]),
        "spot_check_combined_sem": r(np.hypot(spot[3], sem[-1])),
        "pass_successive_10pct": bool(max(succ) <= 0.10),
        "pass_finest_two_5pct": bool(finest <= 0.05),
        "pass_no_monotone_drift": bool(not (monotone and drift > drift_lim)),
    }
    print(f"  dt arm {conv}", flush=True)

    lles = [r(qr[f"ic_{i}"]["spectrum"][0]) for i in range(5)]
    lles_a = np.array(lles)
    spread = float((lles_a.max() - lles_a.min()) / lles_a.mean())
    res["largest_exponent_across_ics"] = lles
    res["ic_consistency"] = {
        "ics": [rlist(s, 5) for s in ic_pts],
        "protocol": "5 ICs = random points of a short orbit from DEFAULT_STATE + uniform(-0.3,0.3)^3 "
                    "(rng seed 11); QR 150k steps at dt 0.005 after 20k transient steps",
        "spread_over_mean": r(spread, 4),
        "pass": bool(lles_a.min() > 0.10 and spread <= 0.25),
    }
    print(f"  ICs {lles} spread {spread:.3f}", flush=True)

    gap = abs(spec.sum() - mean_div)
    res["sum_lambda_vs_divergence"] = {
        "sum_lambda": r(spec.sum()), "mean_divergence": r(mean_div),
        "abs_difference": r(gap), "relative": r(gap / abs(mean_div), 4),
        "pass": bool(gap <= 0.03 * abs(mean_div)),
    }
    res["qr_job_seconds"] = {k: r(v["seconds"], 1) for k, v in qr.items()}

    # ------------------------------------------------------------------
    # Pass/fail table.
    # ------------------------------------------------------------------
    sa = res["single_attractor_test"]
    dc = res["dt_convergence_summary"]

    def row(cid, crit, value, ok, kind="gate"):
        return {"id": cid, "criterion": crit, "value": value, "pass": bool(ok), "kind": kind}

    T = [
        row("B1a", "lambda1 >= 0.15 (dt 0.005, 600k steps)", res["lyapunov_spectrum"][0], spec[0] >= 0.15),
        row("B1b", "|lambda2| <= 0.01", res["lyapunov_spectrum"][1], abs(spec[1]) <= 0.01),
        row("B1c", "lambda3 < 0 and sum < 0", [res["lyapunov_spectrum"][2], res["lyapunov_sum"]],
            spec[2] < 0 and spec.sum() < 0),
        row("B2a", "dt halving 0.01->0.00125: successive change <= 10%", dc["successive_relative_change"],
            dc["pass_successive_10pct"]),
        row("B2b", "finest two dt within 5%", dc["finest_two_relative_change"], dc["pass_finest_two_5pct"]),
        row("B2c", "no monotone drift beyond 2x combined segment SEM",
            {"monotone": monotone, "drift": dc["drift_coarsest_to_finest"],
             "limit": dc["drift_limit_2x_combined_sem"]}, dc["pass_no_monotone_drift"]),
        row("B3", "|sum lambda - mean div| <= 3% of |mean div|", res["sum_lambda_vs_divergence"]["relative"],
            res["sum_lambda_vs_divergence"]["pass"]),
        row("B4a", "2 < D_KY < 3", res["kaplan_yorke_dimension"], 2 < ky < 3),
        row("B4b", "D_KY in beauty band 2.05-2.5", res["kaplan_yorke_dimension"], 2.05 <= ky <= 2.5, "beauty"),
        row("B5a", "4M-step orbit finite and max|X| < 10", res["long_orbit"]["max_abs"], res["bounded_4M_steps"]),
        row("B5b", "mean divergence < 0", res["mean_divergence"], mean_div < 0),
        row("B5c", "|centroid_i|/extent_i <= 0.01", res["centroid"]["over_extent"], res["centroid"]["pass"]),
        row("B6", "5 near-attractor ICs: all lambda1 > 0.10, spread <= 25%",
            {"lambda1": lles, "spread": res["ic_consistency"]["spread_over_mean"]},
            res["ic_consistency"]["pass"]),
        row("B7", "equilibria: on-axis triple, D3-orbits of size 1/2/3/6, spectra equal within orbits "
                  "(1e-6), plane points on half-turn axes, on-axis pair counter-rotating",
            {"orbit_sizes": eqchk["orbit_sizes"], "gap": eqchk["max_spectrum_gap_within_orbit"],
             "sense": eqchk["on_axis_pair_rotation_sense"]}, eqchk["pass"]),
        row("B8a", "max D3 residual <= 1.5 x floor and <= 0.015",
            {"max_group": sym["max_group"], "floor": sym["floor_half_vs_half"]}, sym["pass_group_vs_floor"]),
        row("B8b", "every control >= 3x max D3 residual", sym["min_ratio"], sym["pass_controls_3x"]),
        row("B8c", "C6 and y-mirror controls >= 5x",
            [sym["control_ratios"]["C6"], sym["control_ratios"]["y_mirror"]], sym["pass_c6_ymirror_5x"]),
        row("B9a-strict", "48/48 random ICs in [-2,2]^3 land on the symmetric set (pre-registered)",
            f"{sa['basin']['n_landed']}/48 landed, {sa['basin']['n_escaped_T300']} escaped",
            sa["basin"]["pass_strict_48_of_48"], "gate (relaxed by user decision)"),
        row("B9a-relaxed", "every bounded orbit lands on the symmetric set (user-accepted gate)",
            f"{sa['basin']['n_landed']}/{sa['basin']['n_bounded_T300']}",
            sa["basin"]["pass_relaxed_every_bounded_lands"]),
        row("B9a-large", "every bounded orbit of the 480-IC sample lands on the symmetric set",
            f"{sa['basin_large_sample']['n_landed']}/{sa['basin_large_sample']['n_bounded_T300']} "
            f"(bounded fraction {sa['basin_large_sample']['bounded_fraction']} "
            f"+- {sa['basin_large_sample']['bounded_fraction_binomial_se']})",
            sa["basin_large_sample"]["every_bounded_lands"]),
        row("B9b", "48/48 near-attractor ICs bounded and landed",
            f"{sa['basin_near_attractor']['n_landed']}/48", sa["basin_near_attractor"]["pass"]),
        row("B9c-1", "image launch: early equivariance <= 1e-8 (kappa, rho)",
            [sa["image_launch"]["early_equivariance_kappa"], sa["image_launch"]["early_equivariance_rho"]],
            sa["image_launch"]["pass_equivariance"]),
        row("B9c-2", "image launch: cross-cloud residuals <= 1.5 x floor",
            {"max": max(sa["image_launch"]["cloud_residuals"].values()),
             "limit": sa["image_launch"]["cloud_residual_limit"]},
            sa["image_launch"]["pass_cloud_residuals"]),
        row("B9c-3", "image launch: fingerprint distances < 0.0549 (v1, reported only)",
            sa["image_launch"]["fingerprint_distances"], sa["image_launch"]["pass_fingerprint_v1"], "report"),
        row("B9d-1", "|P(z>0) - 0.5| <= 0.01", zs["p_positive"], zs["pass_p_positive"]),
        row("B9d-2", "KS(z, -z) <= 0.02", zs["ks_z_vs_minus_z"], zs["pass_ks"]),
        row("B9d-3", "counter-rotation: omega above/below opposite, within 5% of mean speed or 3 SE",
            {"above": zs["omega_above"], "below": zs["omega_below"],
             "mismatch_over_speed": zs["omega_mismatch_over_mean_speed"],
             "mismatch_over_se": zs["omega_mismatch_over_se"]}, zs["pass_counter_rotation"]),
        row("B9d-3v1", "counter-rotation within 5% of mean speed (v1, strict)",
            zs["omega_mismatch_over_mean_speed"], zs["pass_counter_rotation_v1"], "report"),
        row("B10a", "nearest catalog distance > 0.110 (interpretable)", nov["distance"], nov["interpretable"],
            "interpretation"),
        row("B10b", "nearest catalog distance > 0.115 (distinct)", nov["distance"], nov["distinct"],
            "interpretation"),
        row("B10c", "distance to Aurelia > 0.110 before shape-novelty wording",
            res["novelty_vs_atlas"]["distances"]["Aurelia"], res["novelty_vs_atlas"]["aurelia_above_0.110"],
            "interpretation"),
    ]
    gates = [t for t in T if t["kind"] == "gate"]
    res["pass_fail_table"] = T
    res["verdict"] = {
        "all_gates_pass": bool(all(t["pass"] for t in gates)),
        "failed_gates": [t["id"] for t in gates if not t["pass"]],
        "failed_other": [t["id"] for t in T if t["kind"] != "gate" and not t["pass"]],
        "basin_disclosure": f"{sa['basin']['n_bounded_T300']}/48 random ICs in [-2,2]^3 stay bounded "
                            f"and {sa['basin']['n_landed']} land on the symmetric set; "
                            f"{sa['basin']['n_escaped_T300']} escape to infinity (shared protocol, rng seed 0); "
                            f"{sa['basin_large_sample']['n_ic']}-IC sample: bounded fraction {sa['basin_large_sample']['bounded_fraction']} "
                            f"+- {sa['basin_large_sample']['bounded_fraction_binomial_se']}, "
                            f"{sa['basin_large_sample']['n_landed']}/{sa['basin_large_sample']['n_bounded_T300']} "
                            f"bounded orbits land on the symmetric set",
        "certified": bool(all(t["pass"] for t in gates)) and not Q,
        "note": "gates of kind 'gate' decide certification; B9a-strict is the pre-registered "
                "48/48 basin gate, replaced for this specimen by the user-accepted B9a-relaxed "
                "gate with the escape fraction disclosed; 'report', 'beauty' and "
                "'interpretation' rows do not decide certification",
    }
    res["runtime_s"] = r(time.time() - t_start, 1)
    write(res)
    print("\nPASS/FAIL")
    for t in T:
        print(f"  [{'PASS' if t['pass'] else 'FAIL'}] {t['id']:12s} {t['kind']:32s} {t['criterion']}")
    print(f"\nverdict: {res['verdict']}")
    print(f"wrote {OUT}  ({res['runtime_s']} s)")


if __name__ == "__main__":
    main()
