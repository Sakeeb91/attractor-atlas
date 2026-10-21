"""Equilibria of the Trinacria system, their stability, and their D3-orbits.

On the axis x = y = 0 the planar equations vanish and the z-equation reduces to
z*(nu - z^2) = 0, because kappa forbids a constant term in dz/dt (dz/dt must be
odd under (x, y, z) -> (x, -y, -z)). The on-axis equilibria are therefore the
triple (0,0,0) and (0,0,+-sqrt(nu)); the two non-zero ones are exchanged by
kappa and have identical spectra with opposite rotation sense.

Off-axis equilibria are hunted by multi-start root finding. By equivariance
the image of an equilibrium under any of the six group elements is again an
equilibrium with a conjugate Jacobian, so the set splits into D3-orbits:

    size 1   the origin (the only point fixed by all of D3)
    size 2   the on-axis pair (fixed by the rotations, swapped by kappa)
    size 3   points on a half-turn axis in the plane z = 0
    size 6   generic points

``group_orbits`` performs this partition and checks that spectra agree within
each orbit. The origin is a saddle with the double real eigenvalue -beta (the
rotation rate omega*z vanishes at z = 0), not a saddle-focus, so the
homoclinic narrative of the C_n specimens does not carry over.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import fsolve

from .dynamics import (
    ALPHA, BETA, DELTA, EPS, GAMMA, LAM, NU, OMEGA, group_matrices, jacobian,
    velocity,
)

PARAMS = (ALPHA, BETA, OMEGA, GAMMA, DELTA, NU, LAM, EPS)


def find_equilibria(params=PARAMS, n_seeds=600, seed=0, tol=1e-10, box=3.0,
                    complete_orbits=True):
    """Locate all equilibria (on-axis roots + multi-start search off-axis).

    Seeds are uniform in [-box, box]^3; fsolve runs with xtol 1e-13, a root is
    accepted when ||F|| < tol and deduplicated at distance 1e-6. Equivariance
    maps every root to a root, so the six D3 images of each root found are
    then added (each polished by fsolve and re-checked against tol). A
    D3-orbit is therefore complete as soon as one member is hit, which makes
    the result independent of the seed draw (a bare 600-seed search misses a
    member of a size-6 orbit for some seeds). ``complete_orbits=False``
    returns the bare multi-start result, for comparison.
    """
    nu = params[5]
    found = [np.array([0.0, 0.0, 0.0])]
    if nu > 0:
        found.append(np.array([0.0, 0.0, np.sqrt(nu)]))
        found.append(np.array([0.0, 0.0, -np.sqrt(nu)]))

    rng = np.random.default_rng(seed)
    for s0 in rng.uniform(-box, box, size=(n_seeds, 3)):
        sol, _, ier, _ = fsolve(
            lambda s: velocity(s, *params), s0, full_output=True, xtol=1e-13
        )
        if ier != 1 or np.linalg.norm(velocity(sol, *params)) > tol:
            continue
        if all(np.linalg.norm(sol - f) > 1e-6 for f in found):
            found.append(sol)

    for root in (list(found) if complete_orbits else []):
        for g in group_matrices():
            img = g @ root
            sol, _, ier, _ = fsolve(
                lambda s: velocity(s, *params), img, full_output=True, xtol=1e-13
            )
            if np.linalg.norm(velocity(sol, *params)) > tol or np.linalg.norm(sol - img) > 1e-6:
                continue
            if all(np.linalg.norm(sol - f) > 1e-6 for f in found):
                found.append(sol)

    found.sort(key=lambda s: (round(s[2], 8), round(np.hypot(s[0], s[1]), 8),
                              np.arctan2(s[1], s[0])))
    return found


def classify(eq, params=PARAMS):
    """Eigenvalues at an equilibrium and a stability label."""
    eigs = np.linalg.eigvals(jacobian(eq, *params))
    re = np.real(eigs)
    has_complex = np.any(np.abs(np.imag(eigs)) > 1e-9)
    if np.all(re < 0):
        label = "stable focus" if has_complex else "stable node"
    elif np.all(re > 0):
        label = "unstable focus" if has_complex else "unstable node"
    else:
        label = "saddle-focus" if has_complex else "saddle"
    return eigs, label


def rotation_sense(eq, params=PARAMS):
    """J_xy - J_yx at an equilibrium (sign = sense of rotation about z).

    Equals -2*omega*z, so the on-axis pair (0,0,+-sqrt(nu)) carry opposite
    signs: the counter-rotation that kappa enforces.
    """
    J = jacobian(eq, *params)
    return float(J[0, 1] - J[1, 0])


def _sorted_spectrum(eigs):
    eigs = np.asarray(eigs, dtype=complex)
    return eigs[np.lexsort((np.round(eigs.imag, 9), np.round(eigs.real, 9)))]


def group_orbits(eqs, params=PARAMS, match_tol=1e-6):
    """Partition equilibria into D3-orbits.

    Returns a list of dicts, one per orbit, with keys
      members        list of indices into ``eqs``
      size           orbit size (1, 2, 3 or 6 for a complete orbit)
      complete       True when every group image of every member was found
      spectrum       sorted eigenvalues of the first member
      spectrum_gap   max deviation of any member's sorted spectrum from it
      label          stability label of the first member
    """
    eqs = [np.asarray(e, dtype=float) for e in eqs]
    G = group_matrices()
    unassigned = set(range(len(eqs)))
    orbits = []
    while unassigned:
        i = min(unassigned)
        images = [g @ eqs[i] for g in G]
        members, complete = [], True
        for img in images:
            hits = [j for j in range(len(eqs)) if np.linalg.norm(eqs[j] - img) < match_tol]
            if not hits:
                complete = False
                continue
            if hits[0] not in members:
                members.append(hits[0])
        members.sort()
        unassigned -= set(members)
        ref_eigs, label = classify(eqs[i], params)
        ref = _sorted_spectrum(ref_eigs)
        gap = 0.0
        for j in members:
            other = _sorted_spectrum(classify(eqs[j], params)[0])
            gap = max(gap, float(np.max(np.abs(other - ref))))
        orbits.append({
            "members": members,
            "size": len(members),
            "complete": complete,
            "spectrum": ref,
            "spectrum_gap": gap,
            "label": label,
        })
    return orbits
