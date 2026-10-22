"""Select the canonical D3 (chiral dihedral) specimen from the probe nominees.

The D_n family lives in atlas/family_d.py. A uniform probe of the n = 3 box
(3000 samples, two-orbit Benettin at dt = 0.01 and 0.005) left 13 nominees,
stored with their provenance in results/d3_candidates.json. Two-orbit numbers
at coarse dt are a search signal only (docs/METHODOLOGY.md section 9), so this
script re-certifies every nominee from scratch and picks one by the rules
written below. It does not tune parameters: the canonical values are the
4-decimal probe values of the chosen nominee.

Stages

  0  candidates   load results/d3_candidates.json, or with --regenerate re-run
                  the uniform probe that produced it (seeded, identical output).
  A  screen       all 13 nominees x 2 initial conditions, batched tangent-space
                  QR (full spectrum) at dt = 0.005, 0.0025, 0.00125 over the
                  SAME physical window (transient T = 100, window T = 400).
                  Gate: alive everywhere; lambda1 >= 0.10 at every dt; spread
                  across dt <= 25 percent; seed scatter <= 20 percent;
                  on a separate T = 2000 orbit (dt = 0.0025): bounded, and
                  centroid within 2 percent of the extent on every axis (the
                  QR window of T = 400 is too short for this: slowly switching
                  orbits show 2-3 percent finite-time offsets there, recorded
                  as centroid_over_extent_T_window). Shape statistics, symmetry
                  residuals and a pre-rank by the section 6.5 beauty
                  score (lambda* = 0.38, s = 0.25) times a symmetry-separation
                  factor min(control residual) / max(group residual), capped
                  at 5.
  B  certify      every screen survivor: one continued orbit per dt in
                  {0.005, 0.0025, 0.00125} with T = 2000 (400k, 800k, 1.6M
                  steps), full QR spectrum, 10-segment SEM of lambda1, D_KY,
                  mean divergence along the orbit. Gate: lambda1 >= 0.15;
                  successive relative change <= 10 percent; finest two within
                  5 percent; a MONOTONE trend across the three dt with
                  |lambda1(0.005) - lambda1(0.00125)| above 2x the combined
                  segment SEM is a rejection; |lambda2| <= 0.01 at the finest dt;
                  lambda3 < 0; |sum lambda - mean div| <= 3 percent.
  C  symmetry     on a reference orbit per finalist (dt = 0.004, T = 8000 from
                  (0.3, 0.2, 0.1), every 4th step): centroid within 1 percent
                  of the extent; residuals of all five
                  non-identity D3 elements against the self floor (first half
                  vs second half of the orbit) and five controls that are not
                  symmetries (z-flip, y-mirror, C6, inversion, 37-degree turn),
                  cKDTree nearest neighbours on the extent-normalized cloud;
                  a kappa/rho image launch (x0, kappa x0, rho x0) with an
                  early-time equivariance check and cloud-to-cloud residual and
                  fingerprint distances; a 48-IC basin landing test; z-sign
                  statistics and the counter-rotation of the two half-spaces.
  D  choose       among finalists passing B and C: prefer the cleanest D3
                  symmetry (largest control/group separation) with
                  lambda1 >= 0.2 and D_KY in [2.1, 2.4]; strongest chaos is
                  secondary.

  R  refine       only when no finalist passes and the preferred miss fails
                  the basin gate alone: 256 Gaussian mutations (2 percent of
                  each box span) screened for >= 90 percent bounded ICs in
                  [-2, 2]^3, then two-orbit lambda1 and centroid. A search
                  signal; it never certifies anything by itself.

Gate versions. Stage C records two verdicts. v1 is the gate set fixed before
the run. v2 replaces two estimators that turned out unable to test what they
were meant to test, and keeps every other threshold: the image-launch
fingerprint distance (rotation-invariant, so it measures segment sampling
noise, not symmetry) is reported but not gated, and the counter-rotation
ratio uses the field-based angular velocity with batch-means errors instead of
the unwrapped d(theta) per sample, whose variance is unbounded near the z axis.
Both verdicts are in the JSON.

Writes results/d3_selection.json (partial JSON after each stage).

    python3 scripts/search_d3.py            # full run, ~25-35 min
    python3 scripts/search_d3.py --quick    # smoke test, ~1 min
    python3 scripts/search_d3.py --regenerate   # also re-run the probe (~1 min)
    python3 scripts/search_d3.py --resume       # rerun stages C and D only (~10 min)
    python3 scripts/search_d3.py --finalize     # redo stage D + naming from the JSON
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from atlas import family_d  # noqa: E402
from atlas.fingerprint import distance, fingerprint  # noqa: E402

N = 3
CAND_PATH = ROOT / "results" / "d3_candidates.json"
OUT_PATH = ROOT / "results" / "d3_selection.json"
FP_FLOOR = 0.0549  # re-simulation noise floor of the fingerprint distance

# Naming record (section 7.8). Filled in after the selection; searched on the
# date given, WebSearch only (WebFetch fails on DNS in the build environment).
NAMING = {
    "name": "Trinacria",
    "lineage": "marine: Trinacria Heiberg 1863 is an extinct marine centric diatom "
               "genus (family Hemiaulaceae; type species T. regina, Eocene, "
               "Denmark; AlgaeBase lists it as marine). The genus takes the Greek "
               "name of Sicily, 'three-pointed', whose emblem is the triskelion: "
               "three legs turning one way, reversed when the disc is turned over, "
               "which is what the half-turn kappa does to this flow's two halves",
    "previous_specimens": "Aurelia, Cassiopea (jellies), Naiad (nymph), Mobula (ray)",
    "searched_on": "2026-10-21",
    "tool": "WebSearch (WebFetch fails on DNS in the build environment)",
    "queries": [
        {"query": "\"Trinacria attractor\"",
         "result": "no dynamical-systems use; hits were the Sicilian emblem, the "
                   "dragonfly Orthetrum trinacria, and generic attractor pages"},
        {"query": "\"Trinacria\" chaotic system dynamical system attractor",
         "result": "no source names a chaotic system Trinacria"},
        {"query": "\"Trinacria system\" OR \"Trinacria flow\" Lyapunov chaos",
         "result": "no hit uses either term"},
        {"query": "Trinacria diatom genus marine Heiberg triangular valve",
         "result": "confirms Trinacria Heiberg 1863, extinct Eocene diatom genus, "
                   "Hemiaulaceae, type T. regina, marine (AlgaeBase); the valve "
                   "shape was not confirmed by the snippets, so it is not claimed"},
    ],
    "verdict": "unclaimed as an attractor or chaotic-system name as of the search date",
}


# ---------------------------------------------------------------------------
# Real form of the n = 3 field, batched over parameter rows (B, 8).
# ---------------------------------------------------------------------------

def _cols(P):
    return [P[:, i] for i in range(8)]


def vel(S, P):
    al, be, om, ga, de, nu, lam, ep = _cols(P)
    x, y, z = S[:, 0], S[:, 1], S[:, 2]
    c = al * z * z - be
    dx = c * x - om * z * y + ga * (x * x - y * y) + 2 * de * z * x * y
    dy = c * y + om * z * x - 2 * ga * x * y + de * z * (x * x - y * y)
    dz = nu * z - z ** 3 - lam * z * (x * x + y * y) + ep * (3 * x * x * y - y ** 3)
    return np.stack([dx, dy, dz], axis=1)


def jac(S, P):
    al, be, om, ga, de, nu, lam, ep = _cols(P)
    x, y, z = S[:, 0], S[:, 1], S[:, 2]
    c = al * z * z - be
    J = np.empty((len(S), 3, 3))
    J[:, 0, 0] = c + 2 * ga * x + 2 * de * z * y
    J[:, 0, 1] = -om * z - 2 * ga * y + 2 * de * z * x
    J[:, 0, 2] = 2 * al * z * x - om * y + 2 * de * x * y
    J[:, 1, 0] = om * z - 2 * ga * y + 2 * de * z * x
    J[:, 1, 1] = c - 2 * ga * x - 2 * de * z * y
    J[:, 1, 2] = 2 * al * z * y + om * x + de * (x * x - y * y)
    J[:, 2, 0] = -2 * lam * z * x + 6 * ep * x * y
    J[:, 2, 1] = -2 * lam * z * y + 3 * ep * (x * x - y * y)
    J[:, 2, 2] = nu - 3 * z * z - lam * (x * x + y * y)
    return J


def div(S, P):
    al, be, om, ga, de, nu, lam, ep = _cols(P)
    x, y, z = S[:, 0], S[:, 1], S[:, 2]
    return 2 * (al * z * z - be) + nu - 3 * z * z - lam * (x * x + y * y)


def rk4(S, P, dt):
    k1 = vel(S, P)
    k2 = vel(S + 0.5 * dt * k1, P)
    k3 = vel(S + 0.5 * dt * k2, P)
    k4 = vel(S + dt * k3, P)
    return S + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)


def tangent_rk4(S, Q, P, dt):
    k1 = vel(S, P); j1 = jac(S, P) @ Q
    S2 = S + 0.5 * dt * k1
    k2 = vel(S2, P); j2 = jac(S2, P) @ (Q + 0.5 * dt * j1)
    S3 = S + 0.5 * dt * k2
    k3 = vel(S3, P); j3 = jac(S3, P) @ (Q + 0.5 * dt * j2)
    S4 = S + dt * k3
    k4 = vel(S4, P); j4 = jac(S4, P) @ (Q + dt * j3)
    return (S + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4),
            Q + dt / 6.0 * (j1 + 2 * j2 + 2 * j3 + j4))


def qr_spectrum(P, S0, dt, T_transient, T_window, n_seg=10, keep_every=0,
                escape=50.0):
    """Batched tangent-space QR Lyapunov spectrum on one continued orbit.

    Returns dict with spectrum (B, 3), seg_l1 (B, n_seg), alive (B,),
    mean_div (B,), and, if keep_every > 0, the orbit samples (n, B, 3).
    """
    P = np.asarray(P, float)
    S = np.array(S0, float)
    B = len(S)
    alive = np.ones(B, bool)

    def cull(S):
        bad = ~np.all(np.isfinite(S), axis=1) | (np.abs(S).max(axis=1) > escape)
        alive[bad] = False
        S[bad] = 0.3
        return S

    with np.errstate(all="ignore"):
        nt = int(round(T_transient / dt))
        for i in range(nt):
            S = rk4(S, P, dt)
            if i % 500 == 0:
                S = cull(S)
        S = cull(S)
        nw = int(round(T_window / dt))
        seg_len = nw // n_seg
        Q = np.broadcast_to(np.eye(3), (B, 3, 3)).copy()
        sums = np.zeros((B, 3))
        seg = np.zeros((B, n_seg))
        dsum = np.zeros(B)
        pts = []
        for i in range(nw):
            S, Q = tangent_rk4(S, Q, P, dt)
            Q, R = np.linalg.qr(Q)
            lg = np.log(np.abs(np.diagonal(R, axis1=1, axis2=2)))
            sums += lg
            k = min(i // seg_len, n_seg - 1)
            seg[:, k] += lg[:, 0]
            dsum += div(S, P)
            if keep_every and i % keep_every == 0:
                pts.append(S.copy())
            if i % 500 == 0:
                bad = ~np.all(np.isfinite(S), axis=1) | (np.abs(S).max(axis=1) > escape)
                if bad.any():
                    alive[bad] = False
                    S[bad] = 0.3
                    Q[bad] = np.eye(3)
        bad = ~np.all(np.isfinite(S), axis=1) | (np.abs(S).max(axis=1) > escape)
        alive[bad] = False
    spec = sums / (nw * dt)
    seg_l1 = seg / (seg_len * dt)
    seg_l1[:, -1] = seg[:, -1] / ((nw - seg_len * (n_seg - 1)) * dt)
    out = {"spectrum": np.where(alive[:, None], spec, np.nan),
           "seg_l1": seg_l1, "alive": alive, "mean_div": dsum / nw}
    if keep_every:
        out["points"] = np.array(pts)
    return out


def kaplan_yorke(spec):
    e = np.sort(np.asarray(spec))[::-1]
    cum = 0.0
    for j, lam in enumerate(e):
        if cum + lam < 0.0:
            return j + cum / abs(lam) if j else 0.0
        cum += lam
    return float(len(e))


# ---------------------------------------------------------------------------
# Symmetry statistics.
# ---------------------------------------------------------------------------

def _rotz(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


GROUP = family_d.group_matrices(N)
ELEMENT_NAMES = ["rho", "rho2", "kappa", "rho_kappa", "rho2_kappa"]
KAPPA = np.diag([1.0, -1.0, -1.0])


def named_group():
    r = _rotz(2 * np.pi / 3)
    mats = {"rho": r, "rho2": r @ r, "kappa": KAPPA,
            "rho_kappa": r @ KAPPA, "rho2_kappa": r @ r @ KAPPA}
    for m in mats.values():
        assert abs(np.linalg.det(m) - 1.0) < 1e-12
        assert any(np.allclose(m, g, atol=1e-12) for g in GROUP)
    return mats


CONTROLS = {
    "z_flip": np.diag([1.0, 1.0, -1.0]),
    "y_mirror": np.diag([1.0, -1.0, 1.0]),
    "C6": _rotz(np.pi / 3),
    "inversion": -np.eye(3),
    "rot37": _rotz(np.deg2rad(37.0)),
}


def residuals(pts, n_base=20_000, n_query=4_000, seed=0):
    """Extent-normalized NN-RMS residual of each group element and control.

    The cloud is centred on the ORIGIN (the only D3-fixed point), not on its
    own mean, so an off-centre cloud is penalised rather than hidden.
    """
    q = pts[np.all(np.isfinite(pts), axis=1)]
    extent = float(np.ptp(q, axis=0).max())
    rng = np.random.default_rng(seed)
    half = len(q) // 2
    A, Bh = q[:half], q[half:]

    def sub(X, k):
        return X[rng.choice(len(X), size=min(k, len(X)), replace=False)]

    base = sub(q, n_base)
    tree = cKDTree(base)
    qi = sub(q, n_query)

    def rms(tree_, X):
        d, _ = tree_.query(X)
        return float(np.sqrt(np.mean(d ** 2)) / extent)

    # self floor: half vs half, both directions, matched base size
    ta, tb = cKDTree(sub(A, n_base // 2)), cKDTree(sub(Bh, n_base // 2))
    floor = 0.5 * (rms(tb, sub(A, n_query)) + rms(ta, sub(Bh, n_query)))
    ident = rms(tree, qi)  # identity on the same base: the sampling floor
    grp = {k: rms(tree, qi @ m.T) for k, m in named_group().items()}
    ctl = {k: rms(tree, qi @ m.T) for k, m in CONTROLS.items()}
    gmax = max(grp.values())
    return {
        "group": {k: round(v, 5) for k, v in grp.items()},
        "identity_same_base": round(ident, 5),
        "floor_half_vs_half": round(floor, 5),
        "controls": {k: round(v, 5) for k, v in ctl.items()},
        "max_group": round(gmax, 5),
        "min_control": round(min(ctl.values()), 5),
        "separation": round(min(ctl.values()) / gmax, 3),
        "group_over_floor": round(gmax / floor, 3),
        "c6_ratio": round(ctl["C6"] / gmax, 3),
        "y_mirror_ratio": round(ctl["y_mirror"] / gmax, 3),
    }


def shape_stats(pts):
    q = pts[np.all(np.isfinite(pts), axis=1)]
    ext = np.ptp(q, axis=0)
    cen = q.mean(0)
    ev = np.sort(np.linalg.eigvalsh(np.cov((q - cen).T)))[::-1]
    return {
        "centroid": [round(float(v), 5) for v in cen],
        "extent": [round(float(v), 4) for v in ext],
        "centroid_over_extent": [round(float(abs(c) / e), 5) for c, e in zip(cen, ext)],
        "iso": round(float(ext.min() / ext.max()), 4),
        "e3_over_e1": round(float(ev[2] / ev[0]), 4),
        "max_abs": round(float(np.abs(q).max()), 4),
    }


def beauty(l1, iso, e31, lam_star=0.38, s=0.25):
    return float(np.exp(-(l1 - lam_star) ** 2 / (2 * s * s)) * (0.3 + iso) * (0.3 + e31))


# ---------------------------------------------------------------------------
# Complex batched integration (family_d) for launches and basin tests.
# ---------------------------------------------------------------------------

def complex_orbits(p, S0, dt, T_transient, T_sample, every=10, record_first=0):
    S0 = np.asarray(S0, float)
    B = len(S0)
    P = np.asarray(p, float)
    if P.ndim == 1:
        P = np.tile(P, (B, 1))
    W = S0[:, 0] + 1j * S0[:, 1]
    Z = S0[:, 2].copy()
    early = []
    with np.errstate(all="ignore"):
        for i in range(int(round(T_transient / dt))):
            if i < record_first:
                early.append(np.stack([W.real, W.imag, Z], 1))
            W, Z = family_d._rk4(W, Z, P, N, dt)
        pts = []
        for i in range(int(round(T_sample / dt))):
            W, Z = family_d._rk4(W, Z, P, N, dt)
            if i % every == 0:
                pts.append(np.stack([W.real, W.imag, Z], 1))
    pts = np.array(pts).transpose(1, 0, 2)  # (B, n, 3)
    alive = np.all(np.isfinite(pts), axis=(1, 2)) & (np.nan_to_num(np.abs(pts), nan=1e9).max(axis=(1, 2)) < family_d.ESCAPE)
    return pts, alive, (np.array(early).transpose(1, 0, 2) if record_first else None)


def cloud_residual(ref, X, extent):
    d, _ = cKDTree(ref).query(X)
    return float(np.sqrt(np.mean(d ** 2)) / extent)


def symmetry_tests(p, ref_pts, quick=False, seed=7, floor_hint=None):
    """Image launch, basin landing and z statistics for one parameter vector."""
    rng = np.random.default_rng(seed)
    ref = ref_pts[np.all(np.isfinite(ref_pts), axis=1)]
    extent = float(np.ptp(ref, axis=0).max())
    ref_base = ref[rng.choice(len(ref), size=min(20_000, len(ref)), replace=False)]
    r = _rotz(2 * np.pi / 3)
    out = {}

    # (b) image launch: x0, kappa x0, rho x0
    x0 = np.array([0.31, 0.17, 0.42])
    S0 = np.stack([x0, KAPPA @ x0, r @ x0])
    T_s = 200.0 if quick else 2000.0
    pts, alive, early = complex_orbits(p, S0, 0.004, 100.0, T_s, every=10, record_first=2000)
    eq_k = float(np.abs(early[1] - early[0] @ KAPPA.T).max())
    eq_r = float(np.abs(early[2] - early[0] @ r.T).max())
    fps = [fingerprint(c[::3]) for c in pts]  # ~16.7k points, as for the floor
    names = ["x0", "kappa_x0", "rho_x0"]
    fpd, crs = {}, {}
    for i in range(3):
        for j in range(i + 1, 3):
            key = f"{names[i]}|{names[j]}"
            fpd[key] = round(float(distance(fps[i], fps[j])), 4) if fps[i] is not None and fps[j] is not None else None
            crs[key] = round(cloud_residual(pts[i], pts[j][::4], extent), 5)
    # Over the whole run the image orbits separate (roundoff, amplified by the
    # chaos), so this number is reported, not gated; the early window is the test.
    full_k = float(np.nanmax(np.abs(pts[1] - pts[0] @ KAPPA.T)))
    fp_ok = all(v is not None and v < FP_FLOOR for v in fpd.values())
    cr_ok = max(crs.values()) <= 1.5 * floor_hint if floor_hint else True
    out["image_launch"] = {
        "alive": alive.tolist(),
        "early_equivariance_kappa": eq_k, "early_equivariance_rho": eq_r,
        "full_run_equivariance_kappa": full_k,
        "fingerprint_distances": fpd, "fingerprint_floor": FP_FLOOR,
        "cloud_residuals": crs, "cloud_residual_limit": 1.5 * floor_hint if floor_hint else None,
        # v1 (pre-registered) also required every fingerprint distance below the
        # generic floor. The fingerprint is (nearly) invariant under rotations,
        # and after the early window the three orbits are just different
        # segments of the same attractor, so the distance measures segment-to-
        # segment sampling noise (which for this attractor can sit slightly
        # above the cross-system floor) and cannot detect a symmetry failure;
        # v2 reports it and gates on equivariance plus cloud-to-cloud residuals.
        "pass_v1": bool(alive.all() and eq_k < 1e-8 and eq_r < 1e-8 and fp_ok),
        "pass": bool(alive.all() and eq_k < 1e-8 and eq_r < 1e-8 and cr_ok),
    }

    # (a) basin landing: 48 random ICs in [-2, 2]^3
    n_ic = 16 if quick else 48
    S0 = rng.uniform(-2, 2, (n_ic, 3))
    pts, alive, _ = complex_orbits(p, S0, 0.004, 200.0, 50.0 if quick else 100.0, every=10)
    rid, rk = [], []
    for k in range(n_ic):
        if not alive[k]:
            rid.append(None); rk.append(None)
            continue
        c = pts[k]
        rid.append(round(cloud_residual(ref_base, c, extent), 5))
        rk.append(round(cloud_residual(ref_base, c @ KAPPA.T, extent), 5))
    landed = sum(1 for v in rid if v is not None and v <= 0.02)
    n_alive = int(alive.sum())
    out["basin"] = {"n_ic": n_ic, "n_alive": n_alive, "n_escaped": n_ic - n_alive,
                    "n_landed": landed,
                    "landed_fraction_of_bounded": round(landed / n_alive, 4) if n_alive else None,
                    "rho_id": rid, "rho_kappa": rk, "threshold": 0.02,
                    "pass": landed == n_ic,
                    "pass_multistable": landed >= 0.9 * n_ic}

    # (c) z statistics and counter-rotation on the reference orbit, with
    # batch-means standard errors (20 batches).
    z = ref[:, 2]
    zz = z[:: max(1, len(z) // 10_000)]
    from scipy.stats import ks_2samp
    ks = float(ks_2samp(zz, -zz).statistic)
    th = np.unwrap(np.arctan2(ref[:, 1], ref[:, 0]))
    dth = np.diff(th)
    zm = 0.5 * (z[1:] + z[:-1])
    om_up = float(dth[zm > 0].mean())
    om_dn = float(dth[zm < 0].mean())
    # angular velocity from the vector field: omega = (x ydot - y xdot) / r^2
    V = vel(ref, np.tile(np.asarray(p, float), (len(ref), 1)))
    lz = ref[:, 0] * V[:, 1] - ref[:, 1] * V[:, 0]
    omf = lz / (ref[:, 0] ** 2 + ref[:, 1] ** 2)
    up = z > 0
    nb = 20
    L = len(z) // nb

    def bm(vals, mask=None):
        v = []
        for k in range(nb):
            sl = slice(k * L, (k + 1) * L)
            x = vals[sl] if mask is None else vals[sl][mask[sl]]
            v.append(x.mean())
        v = np.array(v)
        return float(v.mean()), float(v.std(ddof=1) / np.sqrt(nb))

    pu = bm(up.astype(float))
    wu, wd = bm(omf, up), bm(omf, ~up)
    lu, ld = bm(lz, up), bm(lz, ~up)
    se_w = float(np.hypot(wu[1], wd[1]))
    mismatch = abs(wu[0] + wd[0])
    sgn = np.sign(z)
    flips = np.flatnonzero(np.diff(sgn) != 0)
    runs = np.diff(np.concatenate([[0], flips + 1, [len(z)]]))
    zs_ok = abs(np.mean(z > 0) - 0.5) <= 0.01 and ks <= 0.02
    out["z_symmetry"] = {
        "p_positive": round(float(np.mean(z > 0)), 4), "p_positive_se": round(pu[1], 4),
        "ks_z_vs_minus_z": round(ks, 4),
        "dtheta_per_sample_above": om_up, "dtheta_per_sample_below": om_dn,
        "counter_rotation_ratio": round(-om_up / om_dn, 4) if om_dn != 0 else None,
        "omega_above": [round(wu[0], 5), round(wu[1], 5)],
        "omega_below": [round(wd[0], 5), round(wd[1], 5)],
        "omega_ratio": round(-wu[0] / wd[0], 4) if wd[0] != 0 else None,
        "omega_mismatch_over_se": round(mismatch / se_w, 3) if se_w > 0 else None,
        "Lz_above": [round(lu[0], 5), round(lu[1], 5)],
        "Lz_below": [round(ld[0], 5), round(ld[1], 5)],
        "longest_dwell_fraction": round(float(runs.max() / len(z)), 4),
        "n_sign_changes": int(len(flips)),
        # v1 (pre-registered): unwrapped d(theta) per sample, ratio within 5 percent.
        # That estimator has unbounded variance near the z axis (theta jumps by up
        # to pi between samples), so v2 uses the field-based angular velocity with
        # batch-means errors: opposite signs, and the mismatch within 5 percent of
        # the mean speed OR within 3 standard errors.
        "pass_v1": bool(zs_ok and om_up * om_dn < 0 and abs(-om_up / om_dn - 1) <= 0.05),
        "pass": bool(zs_ok and wu[0] * wd[0] < 0
                     and (mismatch <= 0.05 * 0.5 * (abs(wu[0]) + abs(wd[0]))
                          or mismatch <= 3 * se_w)),
    }
    return out


# ---------------------------------------------------------------------------
# Stage R: local basin refinement (only when the preferred miss fails the
# basin gate alone).
# ---------------------------------------------------------------------------

def refine_basin(p0, index, quick=False, n_mut=256, sigma=0.02, n_ic=64, seed=1):
    """Gaussian mutations of p0 (sigma = 2 percent of each box span); for each,
    the fraction of n_ic random ICs in [-2, 2]^3 that stay bounded over T = 120
    (dt = 0.004). Mutants that keep >= 90 percent bounded get a two-orbit
    lambda1 (dt = 0.005) and a centroid check, to see whether the basin can be
    enlarged without losing the centred chaotic attractor. A search signal
    only: any refined nominee would need the full stage A-C run again.
    """
    if quick:
        n_mut, n_ic = 32, 16
    rng = np.random.default_rng(seed)
    span = family_d.PARAM_HI - family_d.PARAM_LO
    M = np.clip(p0 + rng.normal(0, sigma, (n_mut, 8)) * span,
                family_d.PARAM_LO, family_d.PARAM_HI)
    M = np.vstack([p0, M])  # row 0 is the parent itself
    ics = rng.uniform(-2, 2, (n_ic, 3))
    P = np.repeat(M, n_ic, 0)
    S = np.tile(ics, (len(M), 1))
    W = S[:, 0] + 1j * S[:, 1]
    Z = S[:, 2].copy()
    alive = np.ones(len(P), bool)
    with np.errstate(all="ignore"):
        for i in range(int(round((40.0 if quick else 120.0) / 0.004))):
            W, Z = family_d._rk4(W, Z, P, N, 0.004)
            if i % 100 == 0:
                bad = (~np.isfinite(W) | ~np.isfinite(Z) | (np.abs(W) > family_d.ESCAPE)
                       | (np.abs(Z) > family_d.ESCAPE))
                alive &= ~bad
                W[bad] = 0.0
                Z[bad] = 0.0
    cap = alive.reshape(len(M), n_ic).mean(1)
    good = [int(k) for k in np.where(cap[1:] >= 0.9)[0] + 1]
    rows = []
    if good:
        r = family_d.evaluate_batch(M[good], N, dt=0.005, transient=20_000,
                                    lyap_steps=60_000, sample_steps=8_000, sample_every=4)
        for j, k in enumerate(good):
            pts = r["points"][j]
            cen = (float(np.max(np.abs(pts.mean(0)) / np.ptp(pts, 0)))
                   if r["alive"][j] else None)
            rows.append({"p": M[k].round(4).tolist(), "capture": float(cap[k]),
                         "lambda1_two_orbit": (float(r["lle"][j]) if r["alive"][j] else None),
                         "centroid_over_extent": cen,
                         "distance_from_parent": round(float(np.linalg.norm((M[k] - p0) / span)), 4)})
    viable = [x for x in rows if x["lambda1_two_orbit"] is not None
              and x["lambda1_two_orbit"] >= 0.15 and x["centroid_over_extent"] is not None
              and x["centroid_over_extent"] <= 0.02]
    return {"parent_index": index, "n_mutants": n_mut, "sigma_frac": sigma, "n_ic": n_ic,
            "parent_capture": float(cap[0]),
            "capture_quantiles_0_25_50_75_100": np.quantile(cap[1:], [0, .25, .5, .75, 1]).round(4).tolist(),
            "n_capture_ge_0.9": len(good), "rows": rows, "n_viable": len(viable),
            "verdict": ("viable refined nominees found; rerun stages A-C on them"
                        if viable else "no mutant keeps a centred chaotic attractor "
                        "with >= 90 percent capture: the partial basin is structural")}


# ---------------------------------------------------------------------------
# Stage 0: the probe (re-run on request).
# ---------------------------------------------------------------------------

def regenerate_probe():
    """Re-run the uniform probe that produced results/d3_candidates.json.

    3000 uniform samples (rng seed 3), two-orbit Benettin at dt = 0.01 over a
    T = 60 transient and T = 200 window (ICs seed 0); nominees with
    lambda1 > 0.05 are re-run at dt = 0.005 (ICs seed 1) and kept if
    lambda1 > 0.05 and the two estimates agree within 50 percent.
    """
    lo, hi = family_d.PARAM_LO, family_d.PARAM_HI

    def lle(P, dt, Tt, Tw, seed):
        B = len(P)
        rng = np.random.default_rng(seed)
        w = rng.uniform(-.5, .5, B) + 1j * rng.uniform(-.5, .5, B) + (.3 + .2j)
        z = rng.uniform(-.5, .5, B)
        alive = np.ones(B, bool)

        def cull(w, z):
            bad = ~np.isfinite(w) | ~np.isfinite(z) | (np.abs(w) > 40) | (np.abs(z) > 40)
            alive[bad] = False
            return np.where(bad, 0, w), np.where(bad, 0, z)

        with np.errstate(all="ignore"):
            for i in range(int(Tt / dt)):
                w, z = family_d._rk4(w, z, P, N, dt)
                if i % 200 == 0:
                    w, z = cull(w, z)
            w, z = cull(w, z)
            d0 = 1e-8
            ws, zs = w + d0, z.copy()
            acc = np.zeros(B)
            n = int(Tw / dt)
            for i in range(n):
                w, z = family_d._rk4(w, z, P, N, dt)
                ws, zs = family_d._rk4(ws, zs, P, N, dt)
                dw, dz = ws - w, zs - z
                dist = np.sqrt(np.abs(dw) ** 2 + dz ** 2)
                dist = np.where(dist > 0, dist, d0)
                acc += np.log(dist / d0)
                ws, zs = w + dw * d0 / dist, z + dz * d0 / dist
                if i % 200 == 0:
                    w, z = cull(w, z)
        return np.where(alive, acc / (n * dt), np.nan)

    P = np.random.default_rng(N).uniform(lo, hi, (3000, 8))
    l1 = lle(P, 0.01, 60, 200, 0)
    nom = np.where(l1 > 0.05)[0]
    l2 = lle(P[nom], 0.005, 60, 200, 1)
    ok = (l2 > 0.05) & (np.abs(l2 - l1[nom]) / l1[nom] < 0.5)
    surv = [{"p": P[nom][i].round(4).tolist(), "l1_01": float(l1[nom][i]),
             "l1_005": float(l2[i])} for i in np.where(ok)[0]]
    surv.sort(key=lambda s: -s["l1_005"])
    return surv, {"n_samples": 3000, "alive_fraction": float(np.mean(np.isfinite(l1))),
                  "n_nominees_dt01": int(len(nom)), "n_survivors": len(surv)}


# ---------------------------------------------------------------------------

def decide(res, cands):
    """Stage D from the stored stage B and C rows; returns the preferred miss."""
    dts = res["stage_B_certify"]["dts"]
    cert = {int(k): v for k, v in res["stage_B_certify"]["rows"].items()}
    sym = {int(k): v for k, v in res["stage_C_symmetry"]["rows"].items()}
    certified = [i for i in res["finalists"] if cert[i]["pass"]]
    ok = [i for i in certified if sym[i]["pass"]]

    def pref(i):
        if "residuals" not in sym[i]:
            return (-1, 0.0, 0.0)
        L = float(cert[i]["lambda1_by_dt"][-1])
        dky = float(cert[i]["D_KY"])
        in_band = int(L >= 0.2) + int(2.1 <= dky <= 2.4)
        return (in_band, min(float(sym[i]["residuals"]["separation"]), 10.0), L)

    ranked = sorted(ok, key=pref, reverse=True)
    res["ranking"] = [{"index": i, "key": list(map(float, pref(i)))} for i in ranked]
    res.pop("closest_misses", None)
    if ranked:
        i = ranked[0]
        res["chosen"] = {
            "index": i,
            "parameters": dict(zip(family_d.param_names(), cands[i]["p"])),
            "lambda1_by_dt": cert[i]["lambda1_by_dt"],
            "spectrum_finest": cert[i]["dt"][str(dts[-1])]["spectrum"],
            "D_KY": cert[i]["D_KY"],
            "separation": sym[i]["residuals"]["separation"],
            "multistable": sym[i]["multistable"],
            "box_face_note": "values are the 4-decimal probe values, not optimized; "
                             "a parameter near a box face is irrelevant for a selected point",
        }
        res["verdict"] = f"chosen nominee #{i}"
        return None
    res["chosen"] = None
    res["verdict"] = "FAILED: no finalist passed certification and symmetry gates"
    # Closest misses, ranked by the same preference, with every failed gate.
    res["closest_misses"] = [
        {"index": i, "p": cands[i]["p"], "key": list(map(float, pref(i))),
         "lambda1_by_dt": cert[i]["lambda1_by_dt"], "D_KY": cert[i]["D_KY"],
         "separation": sym[i]["residuals"]["separation"] if "residuals" in sym[i] else None,
         "basin": ({k: sym[i]["basin"][k] for k in ("n_ic", "n_escaped", "n_landed")}
                   if "basin" in sym[i] else None),
         "failed_gates_v2": sym[i]["reject_reasons"],
         "failed_gates_v1": sym[i].get("reject_reasons_v1")}
        for i in sorted(certified, key=pref, reverse=True)]
    return res["closest_misses"][0] if res["closest_misses"] else None


def annotate(res):
    """Attach the naming record and, on failure, the conditional recommendation."""
    res["naming"] = {**NAMING, "applied": bool(res.get("chosen"))}
    if not res.get("chosen"):
        res["naming"]["note"] = ("name searched and reserved, not applied: no "
                                 "nominee passed every gate")
        misses = res.get("closest_misses") or []
        if misses and misses[0]["failed_gates_v2"] and all(
                r.startswith("basin") for r in misses[0]["failed_gates_v2"]):
            m = misses[0]
            b = m["basin"]
            ref = res.get("stage_R_refine", {})
            res["conditional_recommendation"] = {
                "index": m["index"], "p": m["p"],
                "status": "NOT certified by this selection",
                "passes": "stage A screen, stage B certification, and every stage C "
                          "gate except the basin gate (v2); v1 also fails the "
                          "fingerprint image-launch and d(theta) estimators",
                "fails": m["failed_gates_v2"],
                "basin_detail": (f"{b['n_ic'] - b['n_escaped']}/{b['n_ic']} ICs in "
                                 f"[-2, 2]^3 stay bounded and all {b['n_landed']} of "
                                 f"them land on the symmetric set; {b['n_escaped']} "
                                 "escape to infinity"),
                "refinement": ref.get("verdict"),
                "to_accept": "relax the basin gate to 'every bounded orbit lands on "
                             "the symmetric set' and state in the README that the "
                             "attractor coexists with escape to infinity; this is a "
                             "decision for the project owner, not for this script",
            }


def write(res):
    OUT_PATH.write_text(json.dumps(res, indent=1, default=_jsonable) + "\n")


def _jsonable(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def rel_spread(v):
    v = np.asarray(v, float)
    return float((v.max() - v.min()) / abs(v.mean()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="small budgets (smoke test)")
    ap.add_argument("--regenerate", action="store_true", help="re-run the probe first")
    ap.add_argument("--resume", action="store_true",
                    help="reuse stages A and B from results/d3_selection.json, run C and D")
    ap.add_argument("--finalize", action="store_true",
                    help="recompute stage D and the naming record from the stored JSON")
    args = ap.parse_args()
    if args.finalize:
        res = json.loads(OUT_PATH.read_text())
        decide(res, json.loads(CAND_PATH.read_text())["candidates"])
        annotate(res)
        write(res)
        return
    q = args.quick
    t0 = time.time()
    res = {"family": "atlas/family_d.py, n = 3 (chiral dihedral D3)",
           "param_names": family_d.param_names(), "quick": q}

    # ---- Stage 0 ----------------------------------------------------------
    cand_doc = json.loads(CAND_PATH.read_text())
    cands = cand_doc["candidates"]
    if args.regenerate:
        surv, meta = regenerate_probe()
        # The nominee SET is reproduced exactly; the two-orbit lambda1 values
        # differ in the second decimal because chaotic two-orbit estimates
        # amplify floating-point ordering differences, which can reorder ties.
        same = {tuple(s["p"]) for s in surv} == {tuple(c["p"]) for c in cands}
        res["probe_regenerated"] = {**meta, "same_nominee_set_as_stored": same,
                                    "nominees": surv}
        print(f"probe regenerated: {meta}, identical={same}", flush=True)
    P = np.array([c["p"] for c in cands], float)
    K = len(P)
    print(f"[{time.time()-t0:6.0f}s] {K} candidates", flush=True)

    if args.resume:
        prev = json.loads(OUT_PATH.read_text())
        if prev.get("quick") != q or "stage_B_certify" not in prev:
            sys.exit("--resume needs a stored run with the same --quick setting "
                     "that completed stage B")
        for key in ("stage_A_screen", "finalists", "stage_B_certify"):
            res[key] = prev[key]
        res["resumed_after_stage_B"] = True
        finalists = prev["finalists"]
        dts = tuple(prev["stage_B_certify"]["dts"])
        cert = {int(k): v for k, v in prev["stage_B_certify"]["rows"].items()}
        PF = P[finalists]
        SF = np.tile([0.3, 0.2, 0.1], (len(finalists), 1))
        print(f"resumed after stage B: finalists {finalists}", flush=True)
    else:
        # ---- Stage A: batched three-dt screen ---------------------------------
        dts = (0.005, 0.0025, 0.00125)
        Tt, Tw = (20.0, 60.0) if q else (100.0, 400.0)
        seeds = (0, 1)
        ics = []
        for sd in seeds:
            rng = np.random.default_rng(100 + sd)
            ics.append(np.column_stack([0.3 + rng.uniform(-0.3, 0.3, K),
                                        0.2 + rng.uniform(-0.3, 0.3, K),
                                        0.1 + rng.uniform(-0.3, 0.3, K)]))
        S0 = np.concatenate(ics)
        PP = np.concatenate([P] * len(seeds))
        l1 = np.full((len(dts), len(seeds), K), np.nan)
        spec_all = np.full((len(dts), len(seeds), K, 3), np.nan)
        alive_all = np.ones(K, bool)
        orbit = None
        for a, dt in enumerate(dts):
            keep = int(round(0.025 / dt)) if dt == 0.0025 else 0
            r = qr_spectrum(PP, S0, dt, Tt, Tw, keep_every=keep)
            spec = r["spectrum"].reshape(len(seeds), K, 3)
            spec_all[a] = spec
            l1[a] = spec[:, :, 0]
            alive_all &= r["alive"].reshape(len(seeds), K).all(0)
            if keep:
                orbit = r["points"][:, :K, :]  # seed 0
            print(f"[{time.time()-t0:6.0f}s] screen dt={dt}: lambda1 seed0 "
                  f"{np.round(l1[a, 0], 3).tolist()}", flush=True)

        # Shape and symmetry statistics need a longer orbit than the QR window:
        # orbits that switch slowly between the kappa-related halves show
        # finite-time centroid offsets of 2-3 percent of the extent at T = 400
        # that shrink below 1 percent by T = 4000 (nominees #1, #2, #3, #8, #12).
        T_long = 200.0 if q else 2000.0
        lpts, lalive, _ = complex_orbits(P, ics[0], 0.0025, 100.0, T_long, every=20)
        print(f"[{time.time()-t0:6.0f}s] screen long orbits T={T_long}: alive "
              f"{int(lalive.sum())}/{K}", flush=True)

        screen = []
        for k in range(K):
            row = {"index": k, "p": cands[k]["p"],
                   "probe_l1": [round(cands[k]["l1_01"], 4), round(cands[k]["l1_005"], 4)],
                   "alive": bool(alive_all[k]),
                   "lambda1": {str(dt): [round(float(v), 4) for v in l1[a, :, k]]
                               for a, dt in enumerate(dts)},
                   "lambda2_dt00125_seed0": round(float(spec_all[2, 0, k, 1]), 4)}
            reasons = []
            if not alive_all[k] or not np.all(np.isfinite(l1[:, :, k])):
                reasons.append("escaped or non-finite")
            else:
                m = l1[:, :, k].mean(1)
                row["dt_spread"] = round(rel_spread(m), 4)
                row["seed_scatter"] = round(float(np.max(np.abs(l1[:, 0, k] - l1[:, 1, k])
                                                         / np.abs(l1[:, :, k].mean(1)))), 4)
                if np.min(l1[:, :, k]) < 0.10:
                    reasons.append(f"lambda1 min {np.min(l1[:, :, k]):.3f} < 0.10")
                if row["dt_spread"] > 0.25:
                    reasons.append(f"dt spread {row['dt_spread']:.2f} > 0.25")
                if row["seed_scatter"] > 0.20:
                    reasons.append(f"seed scatter {row['seed_scatter']:.2f} > 0.20")
                short = orbit[:, k, :]
                row["centroid_over_extent_T_window"] = shape_stats(short)["centroid_over_extent"]
                pts = lpts[k]
                if not lalive[k]:
                    reasons.append(f"escaped on the T = {T_long} shape orbit")
                elif np.all(np.isfinite(pts)):
                    sh = shape_stats(pts)
                    row["shape"] = sh
                    if max(sh["centroid_over_extent"]) > 0.02:
                        reasons.append(f"off-centre {max(sh['centroid_over_extent']):.3f} > 0.02 "
                                       "(symmetry-broken)")
                    rs = residuals(pts, n_base=8000 if q else 12_000, n_query=2000 if q else 4000)
                    row["residuals"] = rs
                    l1m = float(l1[1, :, k].mean())
                    row["beauty"] = round(beauty(l1m, sh["iso"], sh["e3_over_e1"]), 4)
                    row["score"] = round(row["beauty"] * min(rs["separation"], 5.0), 4)
            row["pass"] = not reasons
            row["reject_reasons"] = reasons
            screen.append(row)
            print(f"  #{k:2d} pass={row['pass']} score={row.get('score')} {reasons}", flush=True)
        res["stage_A_screen"] = {"dts": dts, "T_transient": Tt, "T_window": Tw,
                                 "seeds": len(seeds),
                                 "shape_orbit": {"dt": 0.0025, "T": T_long, "every": 20,
                                                 "ic": "seed-0 IC of the QR screen"},
                                 "rows": screen}
        write(res)

        survivors = sorted([r for r in screen if r["pass"]], key=lambda r: -r["score"])
        finalists = [r["index"] for r in survivors]
        res["finalists"] = finalists
        print(f"[{time.time()-t0:6.0f}s] finalists {finalists}", flush=True)
        if not finalists:
            res["chosen"] = None
            res["verdict"] = "FAILED: no nominee passed the three-dt screen"
            write(res)
            print(res["verdict"])
            return

        # ---- Stage B: long continued-orbit certification ----------------------
        Tw_b = 200.0 if q else 2000.0
        PF = P[finalists]
        SF = np.tile([0.3, 0.2, 0.1], (len(finalists), 1))
        cert = {i: {"dt": {}} for i in finalists}
        for dt in dts:
            r = qr_spectrum(PF, SF, dt, 100.0, Tw_b)
            for j, i in enumerate(finalists):
                seg = r["seg_l1"][j]
                cert[i]["dt"][str(dt)] = {
                    "alive": bool(r["alive"][j]),
                    "spectrum": [round(float(v), 5) for v in r["spectrum"][j]],
                    "segment_sem_l1": round(float(seg.std(ddof=1) / np.sqrt(len(seg))), 5),
                    "mean_div": round(float(r["mean_div"][j]), 5),
                }
            print(f"[{time.time()-t0:6.0f}s] certify dt={dt}: "
                  f"{[cert[i]['dt'][str(dt)]['spectrum'] for i in finalists]}", flush=True)
        for i in finalists:
            c = cert[i]
            reasons = []
            if not all(v["alive"] for v in c["dt"].values()):
                reasons.append("escaped under the long integrator")
            else:
                L = np.array([c["dt"][str(dt)]["spectrum"][0] for dt in dts])
                sem = np.array([c["dt"][str(dt)]["segment_sem_l1"] for dt in dts])
                ch = np.abs(np.diff(L)) / np.abs(L[:-1])
                fin = c["dt"][str(dts[-1])]
                spec = np.array(fin["spectrum"])
                c["lambda1_by_dt"] = L.round(5).tolist()
                c["successive_rel_change"] = ch.round(4).tolist()
                mono = bool(np.all(np.diff(L) > 0) or np.all(np.diff(L) < 0))
                drift = float(abs(L[-1] - L[0]))
                c["monotone"] = mono
                c["drift_over_2sem"] = round(drift / (2 * np.sqrt(sem[0] ** 2 + sem[-1] ** 2)), 3)
                c["D_KY"] = round(kaplan_yorke(spec), 4)
                c["sum_vs_div_rel"] = round(abs(spec.sum() - fin["mean_div"]) / abs(fin["mean_div"]), 4)
                if L.min() < 0.15:
                    reasons.append(f"lambda1 {L.min():.3f} < 0.15")
                if ch.max() > 0.10:
                    reasons.append(f"successive change {ch.max():.3f} > 0.10")
                if ch[-1] > 0.05:
                    reasons.append(f"finest two differ {ch[-1]:.3f} > 0.05")
                if mono and c["drift_over_2sem"] > 1.0:
                    reasons.append(f"monotone drift {drift:.4f} > 2x SEM")
                if abs(spec[1]) > 0.01:
                    reasons.append(f"|lambda2| {abs(spec[1]):.4f} > 0.01")
                if spec[2] >= 0:
                    reasons.append("lambda3 >= 0")
                if c["sum_vs_div_rel"] > 0.03:
                    reasons.append(f"sum vs div {c['sum_vs_div_rel']:.3f} > 0.03")
                if not 2.0 < c["D_KY"] < 3.0:
                    reasons.append(f"D_KY {c['D_KY']} outside (2, 3)")
            c["pass"] = not reasons
            c["reject_reasons"] = reasons
            print(f"  #{i} cert pass={c['pass']} {reasons}", flush=True)
        res["stage_B_certify"] = {"dts": dts, "T_transient": 100.0, "T_window": Tw_b,
                                  "state0": [0.3, 0.2, 0.1],
                                  "rows": {str(i): cert[i] for i in finalists}}
        write(res)

    # ---- Stage C: symmetry on the certified orbit --------------------------
    # Reference orbit: dt = 0.004, T = 8000 from the certification IC, every
    # 4th step (2M points), long enough for 1-percent centroid statistics.
    T_ref = 800.0 if q else 8000.0
    # Only finalists that passed stage B are examined.
    certified = [i for i in finalists if cert[i]["pass"]]
    if not certified:
        res["chosen"] = None
        res["verdict"] = "FAILED: no finalist passed the long-orbit certification"
        write(res)
        print(res["verdict"])
        return
    rpts, ralive, _ = complex_orbits(P[certified], SF[:len(certified)], 0.004, 100.0,
                                     T_ref, every=4)
    print(f"[{time.time()-t0:6.0f}s] reference orbits T={T_ref}: alive {ralive.tolist()}",
          flush=True)
    sym = {}
    for j, i in enumerate(certified):
        pts = rpts[j]
        if not ralive[j]:
            sym[i] = {"pass": False, "reject_reasons": ["reference orbit escaped"],
                      "multistable": None}
            continue
        s = {"shape": shape_stats(pts),
             "residuals": residuals(pts, n_base=8000 if q else 20_000,
                                    n_query=2000 if q else 4000)}
        s.update(symmetry_tests(P[i], pts, quick=q,
                                floor_hint=s["residuals"]["floor_half_vs_half"]))
        rs = s["residuals"]
        reasons = []
        if rs["group_over_floor"] > 1.5:
            reasons.append(f"group residual {rs['group_over_floor']}x floor > 1.5")
        if rs["max_group"] > 0.015:
            reasons.append(f"group residual {rs['max_group']} > 0.015")
        if rs["separation"] < 3:
            reasons.append(f"control separation {rs['separation']} < 3")
        if rs["c6_ratio"] < 5 or rs["y_mirror_ratio"] < 5:
            reasons.append("C6 or y-mirror control < 5x")
        if max(s["shape"]["centroid_over_extent"]) > 0.01:
            reasons.append("centroid > 1 percent of extent")
        if not s["image_launch"]["pass"]:
            reasons.append("image launch failed")
        if not s["basin"]["pass_multistable"]:
            reasons.append(f"basin {s['basin']['n_landed']}/{s['basin']['n_ic']} < 90 percent")
        if not s["z_symmetry"]["pass"]:
            reasons.append("z statistics / counter-rotation failed")
        v1 = [r for r in reasons if not r.startswith(("image launch", "z statistics"))]
        if not s["image_launch"]["pass_v1"]:
            v1.append("image launch failed (v1: fingerprint distance >= generic floor)")
        if not s["z_symmetry"]["pass_v1"]:
            v1.append("z statistics / counter-rotation failed (v1: d(theta) estimator)")
        s["pass_v1"] = not v1
        s["reject_reasons_v1"] = v1
        s["pass"] = not reasons
        s["multistable"] = not s["basin"]["pass"]
        s["reject_reasons"] = reasons
        sym[i] = s
        print(f"[{time.time()-t0:6.0f}s] symmetry #{i}: sep {rs['separation']} "
              f"grp/floor {rs['group_over_floor']} basin {s['basin']['n_landed']}/"
              f"{s['basin']['n_ic']} z {s['z_symmetry']['p_positive']} "
              f"pass={s['pass']} {reasons}", flush=True)
    res["stage_C_symmetry"] = {"reference_orbit": {"dt": 0.004, "T": T_ref, "every": 4,
                                                   "state0": [0.3, 0.2, 0.1]},
                               "rows": {str(i): v for i, v in sym.items()}}
    write(res)

    # ---- Stage D: choose ---------------------------------------------------
    best = decide(res, cands)
    # Stage R: if the preferred miss fails ONLY the basin gate, test whether a
    # local refinement can enlarge its basin without losing the attractor.
    if best and best["failed_gates_v2"] and all(
            r.startswith("basin") for r in best["failed_gates_v2"]):
        res["stage_R_refine"] = refine_basin(np.array(best["p"]), best["index"], q)
        write(res)
    res["notes"] = [
        "sum_vs_div_rel in stage B compares the QR log-volume growth with the "
        "divergence averaged along the SAME orbit (a Liouville consistency check "
        "of the tangent integrator); the independent-orbit check belongs to the "
        "certification battery.",
        "stage A residuals use a 12k-point base on a T = 2000 orbit; stage C uses "
        "a 20k-point base on a T = 8000 orbit, so stage C floors are lower.",
    ]
    annotate(res)
    res["runtime_s"] = round(time.time() - t0, 1)
    write(res)
    print(f"[{time.time()-t0:6.0f}s] {res['verdict']}", flush=True)


if __name__ == "__main__":
    main()
