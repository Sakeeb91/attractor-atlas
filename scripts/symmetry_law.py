"""The chaos-vs-symmetry-group law (Phase 3 of the symmetry program).

Matched-budget quality-diversity search over every equivariant family in the
program, with one robust comparable recorded per group: the **maximum largest
Lyapunov exponent that survives long-integrator certification** under an equal
number of candidate evaluations. The batched two-orbit estimate is optimistic
(it can fabricate chaos around stable equilibria at coarse dt), so the top-K of
each group are re-certified with a fine-dt continued-orbit integrator and the
escapees/artifacts dropped before the maximum is taken.

Groups (and the family that realizes each):
  C2..C8  reducible action on R^3 (plane rotation + height): atlas/family.py
  S4      improper rotoreflection (order 4, Schoenflies):    atlas/family_s4.py
  D3      chiral dihedral (order 6, half-turn kappa):        atlas/family_d.py (n=3)
  T,O,I   polyhedral ROTATION groups, IRREDUCIBLE on R^3:    atlas/family_{t,o,i}.py

The headline finding: chaos strength is non-monotonic across the reducible-action
groups (Cn, S4), which host chaos readily through their saddle-focus engine,
while the irreducible-action polyhedral groups (T, O, I) collapse to ~zero
certified chaos -- the irreducibility obstruction (no commuting linear map but a
scalar -> no linear engine). See docs/SYMMETRY_PROGRAM.md.

Fairness caveat (stated honestly): the families have different parameter counts
(Cn 7, S4 8, D3 8, T 5, O/I 4; read from each family's N_PARAMS), so the
chaotic FRACTION is only loosely comparable across
families. The robust comparable is max certified lambda1 under equal evaluation
budget. Null hypothesis (Letellier-Gilmore): in a cover/image construction the
cover's LLE is independent of the covering group; our families are independently
searched, not covers of a shared image, so any group-dependence is a property of
the equivariant design space, not a cover artifact.

Reproducibility note: run_one seeds its generator from hash(label), and Python
randomizes str hashes per process unless PYTHONHASHSEED is set, so a rerun draws a
different sample of the same box. The D3 row (results/law/law_D3.json) was produced
with PYTHONHASHSEED=0 OMP_NUM_THREADS=1; the earlier rows record no hash seed.

Usage:
  PYTHONHASHSEED=0 python scripts/symmetry_law.py <label>   # one group -> results/law/law_<label>.json
  python scripts/symmetry_law.py aggregate   # build results/symmetry_law.csv + figure
"""

from __future__ import annotations

import json
import pathlib
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

ROOT = pathlib.Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
LAW_DIR = RESULTS / "law"  # one JSON row per group: law/law_<label>.json

# Matched budget and evaluator settings (identical for every group).
N_EVAL = 12_000
TOP_K = 24
SEARCH_DT = 0.005
SEARCH_KW = dict(dt=SEARCH_DT, transient=16_000, lyap_steps=50_000,
                 sample_steps=4_000, sample_every=16)
CERT_DT = 0.0025
CHAOS_THRESH = 0.02  # lambda1 above this (and bounded) counts as chaotic

# group label -> (order, kind, extra). The parameter count is not stored here:
# n_params(label) reads N_PARAMS from the family module, so it cannot drift.
GROUPS = {
    "C2": (2, "cn", 2), "C3": (3, "cn", 3), "C4": (4, "cn", 4),
    "C5": (5, "cn", 5), "C6": (6, "cn", 6), "C7": (7, "cn", 7),
    "C8": (8, "cn", 8),
    "S4": (4, "s4", None),
    "D3": (6, "dn", 3),
    "T": (12, "poly", "t"), "O": (24, "poly", "o"), "I": (60, "poly", "i"),
}

_FAMILY_MODULE = {"cn": "family", "s4": "family_s4", "dn": "family_d"}
_POLY_MODULE = {"t": "family_t", "o": "family_o", "i": "family_i"}


def family_module(label):
    """The atlas family module that realizes group `label`."""
    import importlib
    order, kind, extra = GROUPS[label]
    name = _POLY_MODULE[extra] if kind == "poly" else _FAMILY_MODULE[kind]
    return importlib.import_module(f"atlas.{name}")


def n_params(label):
    """Search-box dimension of group `label`, read from its family's N_PARAMS."""
    return int(family_module(label).N_PARAMS)


# ---------------------------------------------------------------- evaluators
def search_group(label, rng):
    """Run the matched-budget batched search; return (batch_lle, params, chaotic_fraction)."""
    order, kind, extra = GROUPS[label]
    if kind == "cn":
        from atlas import family
        P = rng.uniform(family.PARAM_LO, family.PARAM_HI, size=(N_EVAL, family.N_PARAMS))
        res = family.evaluate_batch(P, extra, seed=int(rng.integers(1 << 30)), **SEARCH_KW)
    elif kind == "s4":
        from atlas import family_s4 as fam
        P = rng.uniform(fam.PARAM_LO, fam.PARAM_HI, size=(N_EVAL, fam.N_PARAMS))
        res = fam.evaluate_batch(P, seed=int(rng.integers(1 << 30)), **SEARCH_KW)
    elif kind == "dn":
        from atlas import family_d as fam
        P = rng.uniform(fam.PARAM_LO, fam.PARAM_HI, size=(N_EVAL, fam.N_PARAMS))
        res = fam.evaluate_batch(P, extra, seed=int(rng.integers(1 << 30)), **SEARCH_KW)
    else:
        fam = {"t": "family_t", "o": "family_o", "i": "family_i"}[extra]
        mod = __import__(f"atlas.{fam}", fromlist=["x"])
        P = rng.uniform(mod.PARAM_LO, mod.PARAM_HI, size=(N_EVAL, mod.N_PARAMS))
        res = mod.evaluate_batch(P, seed=int(rng.integers(1 << 30)), **SEARCH_KW)
    lle = res["lle"]
    alive = res["alive"]
    chaotic = np.isfinite(lle) & alive & (lle > CHAOS_THRESH)
    frac = float(chaotic.mean())
    return lle, P, frac


# real-coordinate single-state velocity adapters (for certification)
def velocity_real(label):
    order, kind, extra = GROUPS[label]
    if kind == "cn":
        n = extra

        def f(s, p):
            x, y, z = s
            al, be, om, ga, mu, nu, lam = p
            w = x + 1j * y
            dw = (al * (z - be) + 1j * om) * w + ga * np.conj(w) ** (n - 1)
            dz = mu + nu * z - z**3 - lam * (x * x + y * y)
            return np.array([dw.real, dw.imag, float(np.real(dz))])
        return f
    if kind == "s4":
        from atlas.family_s4 import velocity as v

        def f(s, p):
            x, y, z = s
            dw, dz = v(x + 1j * y, z, p)
            return np.array([dw.real, dw.imag, float(np.real(dz))])
        return f
    if kind == "dn":
        from atlas.family_d import velocity as v
        n = extra

        def f(s, p):
            x, y, z = s
            dw, dz = v(x + 1j * y, z, p, n)
            return np.array([dw.real, dw.imag, float(np.real(dz))])
        return f
    fam = {"t": "family_t", "o": "family_o", "i": "family_i"}[extra]
    mod = __import__(f"atlas.{fam}", fromlist=["velocity"])
    return lambda s, p: mod.velocity(s, p)


def _rk4(f, s, p, dt):
    k1 = f(s, p); k2 = f(s + 0.5 * dt * k1, p); k3 = f(s + 0.5 * dt * k2, p); k4 = f(s + dt * k3, p)
    return s + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)


def long_lle(f, p, dt=CERT_DT, transient=40_000, steps=120_000, s0=(0.3, 0.2, 0.1)):
    """Continued-orbit two-orbit Benettin LLE at fine dt; None if it escapes."""
    s = np.array(s0, float)
    for _ in range(transient):
        s = _rk4(f, s, p, dt)
        if not np.all(np.isfinite(s)) or np.max(np.abs(s)) > 1e3:
            return None
    s2 = s + np.array([1e-8, 0, 0]); d0 = 1e-8; ls = 0.0
    for _ in range(steps):
        s = _rk4(f, s, p, dt); s2 = _rk4(f, s2, p, dt)
        d = np.linalg.norm(s2 - s); d = max(d, 1e-16)
        ls += np.log(d / d0); s2 = s + (s2 - s) * (d0 / d)
        if not np.all(np.isfinite(s)) or np.max(np.abs(s)) > 1e3:
            return None
    return ls / (steps * dt)


def spectrum_dky(f, p, dt=CERT_DT, transient=30_000, steps=120_000, s0=(0.3, 0.2, 0.1)):
    """Full 3-exponent spectrum via finite-difference-Jacobian QR; Kaplan-Yorke dim."""
    def jac(s):
        h = 1e-6
        J = np.empty((3, 3))
        for k in range(3):
            e = np.zeros(3); e[k] = h
            J[:, k] = (f(s + e, p) - f(s - e, p)) / (2 * h)
        return J
    s = np.array(s0, float)
    for _ in range(transient):
        s = _rk4(f, s, p, dt)
    Q = np.eye(3); sums = np.zeros(3)
    for _ in range(steps):
        J = jac(s)
        k1 = f(s, p); j1 = J @ Q
        s2 = s + 0.5 * dt * k1; J2 = jac(s2); k2 = f(s2, p); j2 = J2 @ (Q + 0.5 * dt * j1)
        s3 = s + 0.5 * dt * k2; J3 = jac(s3); k3 = f(s3, p); j3 = J3 @ (Q + 0.5 * dt * j2)
        s4 = s + dt * k3; J4 = jac(s4); k4 = f(s4, p); j4 = J4 @ (Q + dt * j3)
        s = s + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        Q = Q + dt / 6.0 * (j1 + 2 * j2 + 2 * j3 + j4)
        Q, R = np.linalg.qr(Q); sums += np.log(np.abs(np.diag(R)))
    spec = np.sort(sums / (steps * dt))[::-1]
    cum = 0.0; ky = 0.0
    for j, l in enumerate(spec):
        if cum + l < 0:
            ky = j + cum / abs(l) if j else 0.0
            break
        cum += l
    else:
        ky = float(len(spec))
    return spec.tolist(), float(ky)


def run_one(label):
    t0 = time.time()
    rng = np.random.default_rng(abs(hash(label)) % (1 << 31))
    lle, P, frac = search_group(label, rng)
    order, kind, extra = GROUPS[label]
    f = velocity_real(label)
    # top-K by batch lle, certify with the long integrator
    finite = np.where(np.isfinite(lle), lle, -1e9)
    topk = np.argsort(-finite)[:TOP_K]
    certified = []
    for i in topk:
        if not np.isfinite(lle[i]):
            continue
        l = long_lle(f, P[i])
        if l is not None and l > CHAOS_THRESH:
            certified.append((l, P[i].tolist()))
    certified.sort(key=lambda t: -t[0])
    out = {
        "group": label, "order": order, "n_params": n_params(label), "n_eval": N_EVAL,
        "chaotic_fraction_batch": round(frac, 5),
        "n_certified_of_topK": len(certified),
        "max_lambda1_certified": round(certified[0][0], 5) if certified else 0.0,
        "best_params": certified[0][1] if certified else None,
        "seconds": round(time.time() - t0, 1),
    }
    if certified:
        spec, ky = spectrum_dky(f, np.array(certified[0][1]))
        out["spectrum_at_max"] = [round(v, 4) for v in spec]
        out["D_KY_at_max"] = round(ky, 4)
    else:
        out["spectrum_at_max"] = None
        out["D_KY_at_max"] = None
    LAW_DIR.mkdir(parents=True, exist_ok=True)
    (LAW_DIR / f"law_{label}.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"{label}: maxλ1={out['max_lambda1_certified']} D_KY={out['D_KY_at_max']} "
          f"frac={frac:.4f} certK={len(certified)} ({out['seconds']}s)")
    return out


def aggregate():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = []
    for label in GROUPS:
        p = LAW_DIR / f"law_{label}.json"
        if p.exists():
            rows.append(json.loads(p.read_text()))
    rows.sort(key=lambda r: (r["order"], r["group"]))
    # CSV
    cols = ["group", "order", "n_params", "n_eval", "max_lambda1_certified",
            "D_KY_at_max", "chaotic_fraction_batch", "n_certified_of_topK"]
    lines = [",".join(cols)]
    for r in rows:
        lines.append(",".join(str(r.get(c, "")) for c in cols))
    (RESULTS / "symmetry_law.csv").write_text("\n".join(lines) + "\n")
    print("wrote results/symmetry_law.csv")
    for line in lines:
        print(" ", line)

    # figure: max certified lambda1 vs group, ordered by group order
    fig, ax = plt.subplots(figsize=(10, 6), facecolor="white")
    labels = [r["group"] for r in rows]
    lam = [r["max_lambda1_certified"] for r in rows]
    orders = [r["order"] for r in rows]
    kinds = ["poly" if r["group"] in ("T", "O", "I") else
             ("s4" if r["group"] == "S4" else
              ("dn" if r["group"].startswith("D") else "cn")) for r in rows]
    colors = {"cn": "#1f77b4", "s4": "#9467bd", "dn": "#e07b39", "poly": "#d62728"}
    x = np.arange(len(rows))
    bars = ax.bar(x, lam, color=[colors[k] for k in kinds], edgecolor="black", linewidth=0.6)
    for xi, r in zip(x, rows):
        ax.annotate(f"|G|={r['order']}", (xi, max(r['max_lambda1_certified'], 0) + 0.008),
                    ha="center", fontsize=8, color="#444")
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel(r"max certified $\lambda_1$  (equal eval budget)")
    ax.set_title("Chaos strength vs symmetry group\n"
                 "blue = $C_n$ · purple = $S_4$ rotoreflection · orange = $D_3$ dihedral (all reducible)"
                 " · red = polyhedral (irreducible)")
    ax.axhline(0, color="black", lw=0.8)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=colors["cn"], label=r"$C_n$ rotation (reducible)"),
                       Patch(color=colors["s4"], label=r"$S_4$ rotoreflection"),
                       Patch(color=colors["dn"], label=r"$D_3$ dihedral (chiral)"),
                       Patch(color=colors["poly"], label="polyhedral T/O/I (irreducible)")],
              loc="upper right")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(ROOT / "gallery" / "symmetry_law.png", dpi=160)
    print("wrote gallery/symmetry_law.png")


if __name__ == "__main__":
    arg = sys.argv[1] if len(sys.argv) > 1 else "aggregate"
    if arg == "aggregate":
        aggregate()
    elif arg == "all":
        for label in GROUPS:
            run_one(label)
        aggregate()
    else:
        run_one(arg)
