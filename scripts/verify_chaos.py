"""Full chaos verification of the Aurelia attractor.

Computes and records everything needed to certify a strange attractor:
  * Lyapunov spectrum (long Benettin/QR run) and Kaplan-Yorke dimension
  * equilibria with eigenvalues and stability type
  * boundedness over a long orbit, attractor extent
  * time-averaged divergence (dissipativity)
  * sensitivity to the largest-exponent estimate across initial conditions

Writes results/aurelia_verification.json and prints a summary.
"""

from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from aurelia.dynamics import A, B, C, divergence, trajectory
from aurelia.equilibria import classify, find_equilibria
from aurelia.lyapunov import kaplan_yorke_dimension, lyapunov_spectrum


def main():
    results = {"system": "Aurelia attractor", "parameters": {"a": A, "b": B, "c": C}}

    print("Lyapunov spectrum (long run)...")
    spec = lyapunov_spectrum(n_steps=600_000, dt=0.005)
    ky = kaplan_yorke_dimension(spec)
    results["lyapunov_spectrum"] = [round(v, 5) for v in spec]
    results["kaplan_yorke_dimension"] = round(ky, 5)
    print(f"  spectrum = {np.round(spec, 5)}   D_KY = {ky:.4f}")

    print("Initial-condition robustness of the largest exponent...")
    rng = np.random.default_rng(11)
    lles = []
    for _ in range(5):
        s0 = rng.uniform(-0.5, 0.5, 3) + [0.0, 0.0, 0.5]
        lle = lyapunov_spectrum(n_steps=150_000, dt=0.005, state0=s0)[0]
        lles.append(round(float(lle), 5))
    results["largest_exponent_across_ics"] = lles
    print(f"  LLE across 5 random ICs: {lles}")

    print("Equilibria...")
    eqs = find_equilibria()
    results["equilibria"] = []
    for e in eqs:
        eigs, label = classify(e)
        results["equilibria"].append(
            {
                "point": [round(v, 6) for v in e],
                "eigenvalues": [[round(v.real, 5), round(v.imag, 5)] for v in eigs],
                "type": label,
            }
        )
        print(f"  {np.round(e, 4)}  {label}  eigs={np.round(eigs, 4)}")

    print("Long-orbit boundedness and dissipativity (4 million steps)...")
    pts = trajectory(n_steps=4_000_000, dt=0.004, transient=30_000)
    bounded = bool(np.all(np.isfinite(pts)) and np.max(np.abs(pts)) < 10.0)
    mean_div = float(np.mean([divergence(p) for p in pts[::100]]))
    results["bounded_4M_steps"] = bounded
    results["extent_min"] = [round(v, 4) for v in pts.min(0)]
    results["extent_max"] = [round(v, 4) for v in pts.max(0)]
    results["mean_divergence"] = round(mean_div, 5)
    print(f"  bounded={bounded}  extent=[{np.round(pts.min(0),3)}, {np.round(pts.max(0),3)}]")
    print(f"  time-averaged divergence = {mean_div:.4f} (negative => dissipative)")

    out = pathlib.Path(__file__).resolve().parents[1] / "results" / "aurelia_verification.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(results, indent=2) + "\n")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
