"""Lyapunov convergence study with error bars for the three C_n specimens (Aurelia, Naiad, Cassiopea).

For each of Aurelia, Naiad, and Cassiopea this computes the full Lyapunov
spectrum (Benettin tangent-flow + QR) on a grid:

    dt       in {0.01, 0.005, 0.0025}
    n_steps  in {200_000, 400_000, 800_000}   (read off as checkpoints of
                                               a single 800k-step run)
    10 random initial conditions per (system, dt), fixed across dt

and reports mean +/- std (across initial conditions) for lambda_1, the full
spectrum, and the Kaplan-Yorke dimension, plus segment-based error bars
(the 800k run split into 10 segments) and the lambda_2 -> 0 consistency
check that every autonomous flow must pass.

The integrator is the same classical RK4 tangent scheme as the per-system
``lyapunov`` modules, re-implemented in scalarized form and JIT-compiled
with numba so the full 90-run grid finishes in seconds. Without numba it
falls back to pure Python (slow but identical).

Outputs:
    results/convergence.json
    gallery/convergence.png
"""

from __future__ import annotations

import json
import math
import pathlib
import sys
import time

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

try:
    from numba import njit
except ImportError:  # pragma: no cover - documented fallback, ~100x slower
    def njit(*args, **kwargs):
        if args and callable(args[0]):
            return args[0]
        return lambda f: f

# ---------------------------------------------------------------------------
# Systems. Parameters padded to length 7; sys_id selects the vector field.
#   0 Aurelia    (a, b, c)                      = (1.4, 0.5, 1.9)
#   1 Naiad      (a, b, w, g, mu, nu, lam)      = (1.2, 0.68, 3.7, 1.7, 1.2, 2.0, 2.1)
#   2 Cassiopea  (a, b, w, g, mu, nu, lam)      = (2.0, 0.85, 1.8, 0.9, 1.28, 1.8, 2.8)
# ---------------------------------------------------------------------------

SYSTEMS = [
    ("aurelia", np.array([1.4, 0.5, 1.9, 0.0, 0.0, 0.0, 0.0])),
    ("naiad", np.array([1.2, 0.68, 3.7, 1.7, 1.2, 2.0, 2.1])),
    ("cassiopea", np.array([2.0, 0.85, 1.8, 0.9, 1.28, 1.8, 2.8])),
]

DTS = (0.01, 0.005, 0.0025)
CHECKPOINTS = np.array([200_000, 400_000, 800_000], dtype=np.int64)
N_IC = 10
N_SEGMENTS = 10
TRANSIENT_TIME = 100.0
ESCAPE = 50.0
SEED = 2026


@njit(cache=True)
def _f(sys_id, p, x, y, z):
    if sys_id == 0:
        a, b, c = p[0], p[1], p[2]
        return (
            a * (z - b) * x - c * y + b * (x * x - y * y),
            c * x + a * (z - b) * y - 2.0 * b * x * y,
            a * (1.0 + z) - z * z * z - c * (x * x + y * y),
        )
    elif sys_id == 1:
        a, b, w, g, mu, nu, lam = p[0], p[1], p[2], p[3], p[4], p[5], p[6]
        return (
            (a * (z - b) + g) * x - w * y,
            w * x + (a * (z - b) - g) * y,
            mu + nu * z - z * z * z - lam * (x * x + y * y),
        )
    else:
        a, b, w, g, mu, nu, lam = p[0], p[1], p[2], p[3], p[4], p[5], p[6]
        return (
            a * (z - b) * x - w * y + g * (x * x * x - 3.0 * x * y * y),
            w * x + a * (z - b) * y + g * (y * y * y - 3.0 * x * x * y),
            mu + nu * z - z * z * z - lam * (x * x + y * y),
        )


@njit(cache=True)
def _j(sys_id, p, x, y, z):
    """Jacobian entries (row-major) of the flow at (x, y, z)."""
    if sys_id == 0:
        a, b, c = p[0], p[1], p[2]
        return (
            a * (z - b) + 2.0 * b * x, -c - 2.0 * b * y, a * x,
            c - 2.0 * b * y, a * (z - b) - 2.0 * b * x, a * y,
            -2.0 * c * x, -2.0 * c * y, a - 3.0 * z * z,
        )
    elif sys_id == 1:
        a, b, w, g, mu, nu, lam = p[0], p[1], p[2], p[3], p[4], p[5], p[6]
        return (
            a * (z - b) + g, -w, a * x,
            w, a * (z - b) - g, a * y,
            -2.0 * lam * x, -2.0 * lam * y, nu - 3.0 * z * z,
        )
    else:
        a, b, w, g, mu, nu, lam = p[0], p[1], p[2], p[3], p[4], p[5], p[6]
        return (
            a * (z - b) + 3.0 * g * (x * x - y * y), -w - 6.0 * g * x * y, a * x,
            w - 6.0 * g * x * y, a * (z - b) + 3.0 * g * (y * y - x * x), a * y,
            -2.0 * lam * x, -2.0 * lam * y, nu - 3.0 * z * z,
        )


@njit(cache=True)
def _rk4_base(sys_id, p, x, y, z, dt):
    k1x, k1y, k1z = _f(sys_id, p, x, y, z)
    k2x, k2y, k2z = _f(sys_id, p, x + 0.5 * dt * k1x, y + 0.5 * dt * k1y, z + 0.5 * dt * k1z)
    k3x, k3y, k3z = _f(sys_id, p, x + 0.5 * dt * k2x, y + 0.5 * dt * k2y, z + 0.5 * dt * k2z)
    k4x, k4y, k4z = _f(sys_id, p, x + dt * k3x, y + dt * k3y, z + dt * k3z)
    return (
        x + dt / 6.0 * (k1x + 2.0 * k2x + 2.0 * k3x + k4x),
        y + dt / 6.0 * (k1y + 2.0 * k2y + 2.0 * k3y + k4y),
        z + dt / 6.0 * (k1z + 2.0 * k2z + 2.0 * k3z + k4z),
    )


@njit(cache=True)
def _matmul3(j, Q, out):
    """out = J @ Q with J given as a row-major 9-tuple."""
    j00, j01, j02, j10, j11, j12, j20, j21, j22 = j
    for c in range(3):
        q0, q1, q2 = Q[0, c], Q[1, c], Q[2, c]
        out[0, c] = j00 * q0 + j01 * q1 + j02 * q2
        out[1, c] = j10 * q0 + j11 * q1 + j12 * q2
        out[2, c] = j20 * q0 + j21 * q1 + j22 * q2


@njit(cache=True)
def _spectrum_run(sys_id, p, x0, y0, z0, dt, n_total, transient, checkpoints, n_seg):
    """One long tangent run. Returns (ok, spectra_at_checkpoints, segment_rates)."""
    cp_out = np.zeros((checkpoints.shape[0], 3))
    seg_len = n_total // n_seg
    seg = np.zeros((n_seg, 3))

    x, y, z = x0, y0, z0
    for _ in range(transient):
        x, y, z = _rk4_base(sys_id, p, x, y, z, dt)
    if not (math.isfinite(x) and math.isfinite(y) and math.isfinite(z)):
        return False, cp_out, seg
    if abs(x) > ESCAPE or abs(y) > ESCAPE or abs(z) > ESCAPE:
        return False, cp_out, seg

    Q = np.eye(3)
    K1 = np.empty((3, 3))
    K2 = np.empty((3, 3))
    K3 = np.empty((3, 3))
    K4 = np.empty((3, 3))
    Qt = np.empty((3, 3))
    sums = np.zeros(3)
    ci = 0

    for i in range(n_total):
        # Base-flow RK4 stages (shared with the tangent stages).
        k1x, k1y, k1z = _f(sys_id, p, x, y, z)
        x2, y2, z2 = x + 0.5 * dt * k1x, y + 0.5 * dt * k1y, z + 0.5 * dt * k1z
        k2x, k2y, k2z = _f(sys_id, p, x2, y2, z2)
        x3, y3, z3 = x + 0.5 * dt * k2x, y + 0.5 * dt * k2y, z + 0.5 * dt * k2z
        k3x, k3y, k3z = _f(sys_id, p, x3, y3, z3)
        x4, y4, z4 = x + dt * k3x, y + dt * k3y, z + dt * k3z
        k4x, k4y, k4z = _f(sys_id, p, x4, y4, z4)

        # Tangent RK4: dQ/dt = J(s(t)) Q evaluated at the same stage states.
        _matmul3(_j(sys_id, p, x, y, z), Q, K1)
        for r in range(3):
            for c in range(3):
                Qt[r, c] = Q[r, c] + 0.5 * dt * K1[r, c]
        _matmul3(_j(sys_id, p, x2, y2, z2), Qt, K2)
        for r in range(3):
            for c in range(3):
                Qt[r, c] = Q[r, c] + 0.5 * dt * K2[r, c]
        _matmul3(_j(sys_id, p, x3, y3, z3), Qt, K3)
        for r in range(3):
            for c in range(3):
                Qt[r, c] = Q[r, c] + dt * K3[r, c]
        _matmul3(_j(sys_id, p, x4, y4, z4), Qt, K4)

        for r in range(3):
            for c in range(3):
                Q[r, c] += dt / 6.0 * (K1[r, c] + 2.0 * K2[r, c] + 2.0 * K3[r, c] + K4[r, c])
        x = x + dt / 6.0 * (k1x + 2.0 * k2x + 2.0 * k3x + k4x)
        y = y + dt / 6.0 * (k1y + 2.0 * k2y + 2.0 * k3y + k4y)
        z = z + dt / 6.0 * (k1z + 2.0 * k2z + 2.0 * k3z + k4z)

        # Modified Gram-Schmidt QR; diagonal of R = column norms after projection.
        r0 = math.sqrt(Q[0, 0] ** 2 + Q[1, 0] ** 2 + Q[2, 0] ** 2)
        for r in range(3):
            Q[r, 0] /= r0
        d01 = Q[0, 0] * Q[0, 1] + Q[1, 0] * Q[1, 1] + Q[2, 0] * Q[2, 1]
        for r in range(3):
            Q[r, 1] -= d01 * Q[r, 0]
        r1 = math.sqrt(Q[0, 1] ** 2 + Q[1, 1] ** 2 + Q[2, 1] ** 2)
        for r in range(3):
            Q[r, 1] /= r1
        d02 = Q[0, 0] * Q[0, 2] + Q[1, 0] * Q[1, 2] + Q[2, 0] * Q[2, 2]
        d12 = Q[0, 1] * Q[0, 2] + Q[1, 1] * Q[1, 2] + Q[2, 1] * Q[2, 2]
        for r in range(3):
            Q[r, 2] -= d02 * Q[r, 0] + d12 * Q[r, 1]
        r2 = math.sqrt(Q[0, 2] ** 2 + Q[1, 2] ** 2 + Q[2, 2] ** 2)
        for r in range(3):
            Q[r, 2] /= r2

        l0, l1, l2 = math.log(r0), math.log(r1), math.log(r2)
        sums[0] += l0
        sums[1] += l1
        sums[2] += l2
        s_idx = i // seg_len
        if s_idx >= n_seg:
            s_idx = n_seg - 1
        seg[s_idx, 0] += l0
        seg[s_idx, 1] += l1
        seg[s_idx, 2] += l2

        if ci < checkpoints.shape[0] and i + 1 == checkpoints[ci]:
            t = (i + 1) * dt
            cp_out[ci, 0] = sums[0] / t
            cp_out[ci, 1] = sums[1] / t
            cp_out[ci, 2] = sums[2] / t
            ci += 1

        if i % 1024 == 0:
            if not (math.isfinite(x) and math.isfinite(y) and math.isfinite(z)):
                return False, cp_out, seg
            if abs(x) > ESCAPE or abs(y) > ESCAPE or abs(z) > ESCAPE:
                return False, cp_out, seg

    return True, cp_out, seg / (seg_len * dt)


def kaplan_yorke(spec):
    exps = np.sort(np.asarray(spec))[::-1]
    cum = 0.0
    for j, lam in enumerate(exps):
        if cum + lam < 0.0:
            return 0.0 if j == 0 else j + cum / abs(lam)
        cum += lam
    return float(len(exps))


def main():
    rng = np.random.default_rng(SEED)
    results = {
        "method": "Benettin tangent-flow RK4 + modified Gram-Schmidt QR",
        "dts": list(DTS),
        "n_steps": [int(c) for c in CHECKPOINTS],
        "n_initial_conditions": N_IC,
        "n_segments": N_SEGMENTS,
        "transient_time": TRANSIENT_TIME,
        "seed": SEED,
        "systems": {},
    }

    t_start = time.time()
    for sys_id, (name, params) in enumerate(SYSTEMS):
        ics = np.column_stack(
            [
                rng.uniform(-0.6, 0.6, N_IC),
                rng.uniform(-0.6, 0.6, N_IC),
                rng.uniform(-0.3, 0.8, N_IC),
            ]
        )
        sys_block = {"parameters": params.tolist(), "grid": {}, "escaped_runs": 0}

        for dt in DTS:
            transient = int(round(TRANSIENT_TIME / dt))
            spectra = np.full((N_IC, len(CHECKPOINTS), 3), np.nan)
            seg_rates = np.full((N_IC, N_SEGMENTS, 3), np.nan)
            for k in range(N_IC):
                ok, cp, seg = _spectrum_run(
                    sys_id, params, ics[k, 0], ics[k, 1], ics[k, 2],
                    dt, int(CHECKPOINTS[-1]), transient, CHECKPOINTS, N_SEGMENTS,
                )
                if ok:
                    spectra[k] = cp
                    seg_rates[k] = seg
                else:
                    sys_block["escaped_runs"] += 1

            for j, n in enumerate(CHECKPOINTS):
                sp = spectra[:, j, :]
                valid = ~np.isnan(sp[:, 0])
                sp = sp[valid]
                dky = np.array([kaplan_yorke(s) for s in sp])
                sys_block["grid"][f"dt={dt}|n={int(n)}"] = {
                    "lambda1_mean": float(sp[:, 0].mean()),
                    "lambda1_std": float(sp[:, 0].std(ddof=1)),
                    "lambda2_mean": float(sp[:, 1].mean()),
                    "lambda2_std": float(sp[:, 1].std(ddof=1)),
                    "lambda3_mean": float(sp[:, 2].mean()),
                    "lambda3_std": float(sp[:, 2].std(ddof=1)),
                    "dky_mean": float(dky.mean()),
                    "dky_std": float(dky.std(ddof=1)),
                    "n_valid_ic": int(valid.sum()),
                }
            # Segment-based error bar on the full-length run: std of the
            # per-segment lambda_1 over the 10 segments, averaged across ICs,
            # divided by sqrt(n_segments) to estimate the long-run error.
            seg1 = seg_rates[:, :, 0]
            seg_std = np.nanstd(seg1, axis=1, ddof=1)
            sys_block["grid_segment_error"] = sys_block.get("grid_segment_error", {})
            sys_block["grid_segment_error"][f"dt={dt}"] = {
                "lambda1_segment_std_mean": float(np.nanmean(seg_std)),
                "lambda1_segment_sem": float(np.nanmean(seg_std) / np.sqrt(N_SEGMENTS)),
            }

        results["systems"][name] = sys_block
        print(f"[{name}] done at t+{time.time() - t_start:.1f}s "
              f"(escaped runs: {sys_block['escaped_runs']})")

    # Headline numbers: dt = 0.005 (the certification dt), n = 800k.
    print("\nHeadline (dt=0.005, n=800k, mean +/- std over "
          f"{N_IC} initial conditions):")
    for name in results["systems"]:
        g = results["systems"][name]["grid"]["dt=0.005|n=800000"]
        seg = results["systems"][name]["grid_segment_error"]["dt=0.005"]
        results["systems"][name]["headline"] = {
            "lambda1": g["lambda1_mean"],
            "lambda1_std_ic": g["lambda1_std"],
            "lambda1_sem_segments": seg["lambda1_segment_sem"],
            "lambda2": g["lambda2_mean"],
            "lambda3": g["lambda3_mean"],
            "dky": g["dky_mean"],
            "dky_std_ic": g["dky_std"],
        }
        print(f"  {name:10s} lambda1 = {g['lambda1_mean']:+.4f} +/- {g['lambda1_std']:.4f}"
              f"   lambda2 = {g['lambda2_mean']:+.5f}"
              f"   D_KY = {g['dky_mean']:.4f} +/- {g['dky_std']:.4f}")

    out_json = ROOT / "results" / "convergence.json"
    out_json.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {out_json}")

    plot(results)


def plot(results):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.style.use("dark_background")
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), facecolor="#0d0d14")
    colors = {0.01: "#f4c95d", 0.005: "#6fd3c7", 0.0025: "#c792ea"}
    certified = {"aurelia": (0.2327, 2.155), "naiad": (0.296, 2.110), "cassiopea": (0.525, 2.218)}

    for col, name in enumerate(results["systems"]):
        grid = results["systems"][name]["grid"]
        ax_l, ax_d = axes[0, col], axes[1, col]
        for dt in results["dts"]:
            ns = results["n_steps"]
            l1 = [grid[f"dt={dt}|n={n}"]["lambda1_mean"] for n in ns]
            l1e = [grid[f"dt={dt}|n={n}"]["lambda1_std"] for n in ns]
            dk = [grid[f"dt={dt}|n={n}"]["dky_mean"] for n in ns]
            dke = [grid[f"dt={dt}|n={n}"]["dky_std"] for n in ns]
            ax_l.errorbar(ns, l1, yerr=l1e, marker="o", ms=4, capsize=3,
                          lw=1.2, color=colors[dt], label=f"dt={dt}")
            ax_d.errorbar(ns, dk, yerr=dke, marker="s", ms=4, capsize=3,
                          lw=1.2, color=colors[dt], label=f"dt={dt}")
        ax_l.axhline(certified[name][0], color="w", lw=0.8, ls="--", alpha=0.5)
        ax_d.axhline(certified[name][1], color="w", lw=0.8, ls="--", alpha=0.5)
        ax_l.set_title(name.capitalize(), fontsize=12)
        ax_l.set_xscale("log")
        ax_d.set_xscale("log")
        ax_d.set_xlabel("integration steps")
        if col == 0:
            ax_l.set_ylabel(r"$\lambda_1$")
            ax_d.set_ylabel(r"$D_{KY}$")
        for ax in (ax_l, ax_d):
            ax.set_facecolor("#0d0d14")
            ax.grid(alpha=0.15)
    axes[0, 0].legend(fontsize=8, framealpha=0.2)
    fig.suptitle("Lyapunov convergence: mean ± std over 10 initial conditions; "
                 "dashed = certified value", fontsize=11)
    fig.tight_layout()
    out = ROOT / "gallery" / "convergence.png"
    fig.savefig(out, dpi=160, facecolor=fig.get_facecolor())
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
