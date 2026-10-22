"""Shared basin protocol applied to all five specimens.

Protocol (the one quoted in results/trinacria_verification.json):
  ICs uniform in [-2, 2]^3, rng seed 0 (48 ICs) and rng seed 1 (480 ICs),
  classical RK4 at dt = 0.005 on each package's own velocity (vectorized
  over ICs), escape when a coordinate is non-finite or |X|_inf > 50
  (checked every 10 steps). Bounded counts are taken at T = 120 and T = 300.
  A bounded orbit LANDS on the specimen's attractor when its cloud over
  T in (200, 300] (every 10th step) has extent-normalized NN-RMS distance
  rho_id <= 0.02 to the reference cloud trajectory()[::8] of that package.

Writes results/basin_comparison.json.

    python3 scripts/basin_comparison.py      # about 2 min
"""

from __future__ import annotations

import importlib
import json
import pathlib
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

OUT = ROOT / "results" / "basin_comparison.json"
SPECIMENS = ["aurelia", "naiad", "cassiopea", "mobula", "trinacria"]
DT, ESCAPE = 0.005, 50.0
T_BOUND, T_TOTAL, T_FROM = 120.0, 300.0, 200.0
THRESHOLD = 0.02


def run(mod, S0, ref_tree, extent):
    vel = mod.velocity
    S = np.asarray(S0, float).T.copy()
    B = S.shape[1]
    alive = np.ones(B, bool)
    n_bound, n_total, n_from = (int(round(t / DT)) for t in (T_BOUND, T_TOTAL, T_FROM))
    b120 = None
    cloud = []
    with np.errstate(all="ignore"):
        for i in range(n_total):
            k1 = vel(S)
            k2 = vel(S + 0.5 * DT * k1)
            k3 = vel(S + 0.5 * DT * k2)
            k4 = vel(S + DT * k3)
            S = S + DT / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
            if i % 10 == 0 or i == n_total - 1:
                bad = ~np.all(np.isfinite(S), axis=0) | (np.abs(S).max(axis=0) > ESCAPE)
                alive &= ~bad
                S[:, ~alive] = 0.0
            if i + 1 == n_bound:
                b120 = alive.copy()
            if i >= n_from and (i - n_from) % 10 == 0:
                cloud.append(S.T.copy())
    cloud = np.array(cloud).transpose(1, 0, 2)
    rho = []
    for j in range(B):
        if not alive[j]:
            rho.append(None)
            continue
        d, _ = ref_tree.query(cloud[j])
        rho.append(round(float(np.sqrt(np.mean(d ** 2)) / extent), 5))
    landed = sum(1 for v in rho if v is not None and v <= THRESHOLD)
    nb = int(alive.sum())
    return {
        "n_ic": B, "n_bounded_T120": int(b120.sum()), "n_bounded_T300": nb,
        "n_landed": int(landed),
        "bounded_fraction": round(nb / B, 4),
        "bounded_fraction_binomial_se": round(float(np.sqrt(nb / B * (1 - nb / B) / B)), 4),
        "max_rho_id_of_bounded": max((v for v in rho if v is not None), default=None),
    }


def main():
    t0 = time.time()
    out = {
        "protocol": __doc__.split("Writes")[0].strip(),
        "dt": DT, "escape": ESCAPE, "T_bound": T_BOUND, "T_total": T_TOTAL,
        "landing_window": [T_FROM, T_TOTAL], "landing_threshold": THRESHOLD,
        "specimens": {},
    }
    S48 = np.random.default_rng(0).uniform(-2, 2, (48, 3))
    S480 = np.random.default_rng(1).uniform(-2, 2, (480, 3))
    for name in SPECIMENS:
        mod = importlib.import_module(f"{name}.dynamics")
        ref = mod.trajectory()
        extent = float(np.ptp(ref, axis=0).max())
        tree = cKDTree(ref[::8])
        a = run(mod, S48, tree, extent)
        b = run(mod, S480, tree, extent)
        out["specimens"][name] = {"seed0_48": a, "seed1_480": b}
        print(f"{name:10s} 48: bounded {a['n_bounded_T120']}/{a['n_bounded_T300']} landed {a['n_landed']}  "
              f"480: bounded {b['n_bounded_T300']} ({b['bounded_fraction']}) landed {b['n_landed']}",
              flush=True)
    out["runtime_s"] = round(time.time() - t0, 1)
    OUT.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
