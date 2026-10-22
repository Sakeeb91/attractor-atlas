"""Shilnikov homoclinic hunt for Aurelia, Naiad, and Cassiopea (gap 4).

Each system has a single on-axis saddle-focus whose Jacobian is block
diagonal: a 2-D unstable spiral plane (exactly the plane z = z*) and a 1-D
strongly stable axis direction. A homoclinic orbit must leave along the
unstable plane and return along the stable direction, so we:

  1. shoot trajectories from a small circle of radius r0 in the unstable
     eigenplane, scanning the launch phase theta over one symmetry sector
     (the C_n symmetry makes the rest redundant);
  2. integrate each shot and record the closest return to the equilibrium
     after the trajectory first leaves an r_far-neighborhood;
  3. sweep the rotation-rate parameter (c for Aurelia, w_rot for the
     siblings) across the chaotic band, then zoom twice on the minimum of
     the return distance D(param) = min over theta and time.

A parameter where D drops to ~0 is numerical evidence of a homoclinic
(Shilnikov) bifurcation nearby. The saddle index nu = rho/|lambda_3| < 1
(Shilnikov's inequality, in the time-reversed convention appropriate for a
2-D unstable manifold) is recorded alongside; for all three C_n specimens
(Aurelia, Naiad, Cassiopea) it is parameter-independent along the swept
direction because rho and lambda_3 do not involve the rotation rate.

Outputs:
    results/homoclinic_hunt.json
    gallery/homoclinic.png
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
except ImportError:  # pragma: no cover
    def njit(*args, **kwargs):
        if args and callable(args[0]):
            return args[0]
        return lambda f: f

# (name, padded params, swept index, sweep range, symmetry order n)
SWEEPS = [
    ("aurelia", np.array([1.4, 0.5, 1.9, 0.0, 0.0, 0.0, 0.0]), 2, (0.80, 2.55), 3),
    ("naiad", np.array([1.2, 0.68, 3.7, 1.7, 1.2, 2.0, 2.1]), 2, (1.80, 5.00), 2),
    ("cassiopea", np.array([2.0, 0.85, 1.8, 0.9, 1.28, 1.8, 2.8]), 2, (0.30, 3.00), 4),
]
PARAM_NAMES = {"aurelia": "c", "naiad": "w_rot", "cassiopea": "w_rot"}

DT = 0.0025
T_MAX = 100.0
R0 = 1e-3
R_FAR = 0.5
N_THETA = 48
N_COARSE = 121
ESCAPE = 50.0


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
def _shoot(sys_id, p, eqz, theta, dt, n_steps):
    """Launch from the unstable eigenplane; return (d_min, t_min) of the
    closest approach to the equilibrium after first leaving R_FAR."""
    x = R0 * math.cos(theta)
    y = R0 * math.sin(theta)
    z = eqz
    left = False
    d_min = 1e9
    t_min = -1.0
    for i in range(n_steps):
        k1x, k1y, k1z = _f(sys_id, p, x, y, z)
        k2x, k2y, k2z = _f(sys_id, p, x + 0.5 * dt * k1x, y + 0.5 * dt * k1y, z + 0.5 * dt * k1z)
        k3x, k3y, k3z = _f(sys_id, p, x + 0.5 * dt * k2x, y + 0.5 * dt * k2y, z + 0.5 * dt * k2z)
        k4x, k4y, k4z = _f(sys_id, p, x + dt * k3x, y + dt * k3y, z + dt * k3z)
        x += dt / 6.0 * (k1x + 2.0 * k2x + 2.0 * k3x + k4x)
        y += dt / 6.0 * (k1y + 2.0 * k2y + 2.0 * k3y + k4y)
        z += dt / 6.0 * (k1z + 2.0 * k2z + 2.0 * k3z + k4z)
        if abs(x) > ESCAPE or abs(y) > ESCAPE or abs(z) > ESCAPE:
            break
        dz = z - eqz
        d = math.sqrt(x * x + y * y + dz * dz)
        if not left:
            if d > R_FAR:
                left = True
        elif d < d_min:
            d_min = d
            t_min = (i + 1) * dt
    return d_min, t_min


@njit(cache=True)
def _scan_param(sys_id, p_base, idx, values, eqzs, thetas, dt, n_steps):
    """D(param) = min over theta of the closest-return distance."""
    out_d = np.empty(values.shape[0])
    out_th = np.empty(values.shape[0])
    out_t = np.empty(values.shape[0])
    p = p_base.copy()
    for j in range(values.shape[0]):
        p[idx] = values[j]
        best_d = 1e9
        best_th = 0.0
        best_t = -1.0
        for k in range(thetas.shape[0]):
            d, t = _shoot(sys_id, p, eqzs[j], thetas[k], dt, n_steps)
            if d < best_d:
                best_d = d
                best_th = thetas[k]
                best_t = t
        out_d[j] = best_d
        out_th[j] = best_th
        out_t[j] = best_t
    return out_d, out_th, out_t


def axis_equilibrium(name, p):
    """The single real root of the on-axis z-equation."""
    if name == "aurelia":
        roots = np.roots([-1.0, 0.0, p[0], p[0]])      # a(1+z) - z^3
    else:
        roots = np.roots([-1.0, 0.0, p[5], p[4]])      # mu + nu z - z^3
    real = [r.real for r in roots if abs(r.imag) < 1e-10]
    return float(max(real))


def linearization(name, p, zs):
    """(rho, omega, lambda3) at the on-axis saddle-focus."""
    a, b = p[0], p[1]
    if name == "aurelia":
        rho, omega, lam3 = a * (zs - b), p[2], a - 3.0 * zs ** 2
    elif name == "naiad":
        w, g = p[2], p[3]
        rho = a * (zs - b)
        omega = math.sqrt(max(w * w - g * g, 0.0))
        lam3 = p[5] - 3.0 * zs ** 2
    else:
        rho, omega, lam3 = a * (zs - b), p[2], p[5] - 3.0 * zs ** 2
    return rho, omega, lam3


def hunt(sys_id, name, p_base, idx, lo, hi, n_sym):
    thetas = np.linspace(0.0, 2.0 * np.pi / n_sym, N_THETA, endpoint=False)
    n_steps = int(T_MAX / DT)

    stages = []
    values = np.linspace(lo, hi, N_COARSE)
    coarse = None
    for stage in range(3):
        eqzs = np.array([axis_equilibrium(name, _with(p_base, idx, v)) for v in values])
        d, th, t = _scan_param(sys_id, p_base, idx, values, eqzs, thetas, DT, n_steps)
        j = int(np.argmin(d))
        stages.append({"values": values, "d": d})
        if stage == 0:
            coarse = (values.copy(), d.copy())
        step = values[1] - values[0]
        lo2 = max(lo, values[j] - 2 * step)
        hi2 = min(hi, values[j] + 2 * step)
        best = (float(values[j]), float(d[j]), float(th[j]), float(t[j]))
        values = np.linspace(lo2, hi2, 41)
    return coarse, stages, best


def _with(p, idx, v):
    q = p.copy()
    q[idx] = v
    return q


def best_orbit(sys_id, name, p_base, idx, best):
    """Re-integrate the minimizing shot, recording the orbit for the figure."""
    v, _, theta, t_min = best
    p = _with(p_base, idx, v)
    zs = axis_equilibrium(name, p)
    n_steps = int(min(t_min + 5.0, T_MAX) / DT)
    s = np.array([R0 * math.cos(theta), R0 * math.sin(theta), zs])
    pts = np.empty((n_steps // 4 + 1, 3))
    k = 0
    for i in range(n_steps):
        k1 = np.array(_f(sys_id, p, s[0], s[1], s[2]))
        k2 = np.array(_f(sys_id, p, *(s + 0.5 * DT * k1)))
        k3 = np.array(_f(sys_id, p, *(s + 0.5 * DT * k2)))
        k4 = np.array(_f(sys_id, p, *(s + DT * k3)))
        s = s + DT / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        if i % 4 == 0:
            pts[k] = s
            k += 1
    return pts[:k], zs


def main():
    t0 = time.time()
    results = {
        "method": "shooting from the unstable eigenplane (radius 1e-3, "
                  f"{N_THETA} phases over one C_n sector), RK4 dt={DT}, "
                  f"T_max={T_MAX}, closest return after leaving r={R_FAR}",
        "systems": {},
    }
    figdata = {}
    for sys_id, (name, p_base, idx, (lo, hi), n_sym) in enumerate(SWEEPS):
        zs0 = axis_equilibrium(name, p_base)
        rho, omega, lam3 = linearization(name, p_base, zs0)
        nu_idx = rho / abs(lam3)
        coarse, stages, best = hunt(sys_id, name, p_base, idx, lo, hi, n_sym)
        orbit, zs_best = best_orbit(sys_id, name, p_base, idx, best)

        # D at the canonical parameter, exactly.
        thetas = np.linspace(0.0, 2.0 * np.pi / n_sym, N_THETA, endpoint=False)
        d_canon, _, _ = _scan_param(
            sys_id, p_base, idx, np.array([p_base[idx]]),
            np.array([zs0]), thetas, DT, int(T_MAX / DT))
        # Interior local minima of the coarse sweep below 0.02 (the
        # near-homoclinic "window ladder").
        cv, cd = coarse
        windows = [
            {"value": round(float(cv[i]), 5), "d": round(float(cd[i]), 6)}
            for i in range(1, len(cv) - 1)
            if cd[i] < cd[i - 1] and cd[i] < cd[i + 1] and cd[i] < 0.02
        ]

        results["systems"][name] = {
            "swept_parameter": PARAM_NAMES[name],
            "canonical_value": float(p_base[idx]),
            "sweep_range": [lo, hi],
            "equilibrium_z_canonical": round(zs0, 5),
            "rho_unstable": round(rho, 5),
            "omega_canonical": round(omega, 5),
            "lambda3_stable": round(lam3, 5),
            "saddle_index_nu": round(nu_idx, 4),
            "shilnikov_inequality_nu_lt_1": bool(nu_idx < 1.0),
            "best": {
                PARAM_NAMES[name]: round(best[0], 6),
                "min_return_distance": round(best[1], 6),
                "launch_phase_theta": round(best[2], 5),
                "time_of_closest_pass": round(best[3], 3),
            },
            "min_return_at_canonical": round(float(d_canon[0]), 6),
            "near_homoclinic_windows_interior": windows,
            "coarse_sweep": {
                "values": [round(float(v), 5) for v in coarse[0]],
                "min_return_distance": [round(float(d), 6) for d in coarse[1]],
            },
        }
        figdata[name] = (coarse, stages, best, orbit, zs_best)
        print(f"[{name}] nu={nu_idx:.3f} best {PARAM_NAMES[name]}={best[0]:.5f} "
              f"D_min={best[1]:.2e} (t+{time.time() - t0:.0f}s)")

    out = ROOT / "results" / "homoclinic_hunt.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"wrote {out}")
    plot(results, figdata)


def plot(results, figdata):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.style.use("dark_background")
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8), facecolor="#0d0d14")
    colors = {"aurelia": "#f4c95d", "naiad": "#6fd3c7", "cassiopea": "#c792ea"}

    for col, name in enumerate(figdata):
        (cv, cd), stages, best, orbit, zs = figdata[name]
        r = results["systems"][name]
        c = colors[name]
        ax = axes[0, col]
        ax.set_facecolor("#0d0d14")
        ax.semilogy(cv, cd, color=c, lw=1.1)
        for st in stages[1:]:
            ax.semilogy(st["values"], st["d"], ".", color="w", ms=2, alpha=0.5)
        ax.semilogy([best[0]], [best[1]], "*", color="w", ms=12, zorder=5)
        ax.axvline(r["canonical_value"], color="w", lw=0.8, ls="--", alpha=0.4)
        ax.set_title(f"{name.capitalize()}  (nu={r['saddle_index_nu']}, "
                     f"D$_{{min}}$={best[1]:.1e} at {r['swept_parameter']}={best[0]:.4f})",
                     fontsize=10)
        ax.set_xlabel(r["swept_parameter"])
        if col == 0:
            ax.set_ylabel("closest return to saddle-focus")
        ax.grid(alpha=0.15)

        ax2 = axes[1, col]
        ax2.set_facecolor("#0d0d14")
        ax2.plot(orbit[:, 0], orbit[:, 2], color=c, lw=0.45, alpha=0.9)
        ax2.plot([0], [zs], "x", color="w", ms=9, mew=2, zorder=5)
        ax2.set_xlabel("x")
        if col == 0:
            ax2.set_ylabel("z")
        ax2.set_title("the minimizing shot (x-z), x = saddle-focus", fontsize=9)

    fig.suptitle("Homoclinic hunt: unstable-manifold shooting, closest return vs rotation "
                 "rate; dashed = canonical parameter", fontsize=11)
    fig.tight_layout()
    out = ROOT / "gallery" / "homoclinic.png"
    fig.savefig(out, dpi=160, facecolor=fig.get_facecolor())
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
