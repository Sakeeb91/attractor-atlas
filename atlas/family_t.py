"""The tetrahedral (T) rotation-equivariant search family, evaluated in batches.

This is the first *polyhedral* member of the symmetry program: where the C_n /
S4 families act on the plane w = x+iy with z a separate height coordinate, the
tetrahedral family acts on all of R^3, and the group T (order 12, == A4) acts
**irreducibly**. By Schur, the only linear map commuting with an irreducible
real representation is a scalar, so the linearization at the origin is forced to
be rho*I -- there is **no saddle-focus engine** as in the C_n family. All chaos
must come from the nonlinear equivariant terms; the spin enters only through the
cross term.

The dissipative ansatz spans the *complete* degree-<=3 T-equivariant vector
field (5 dials rho, d, a, b, c), with r^2 = x^2 + y^2 + z^2, the cubic invariant
I_T = x*y*z, and the quartic invariant J = x^4 + y^4 + z^4:

    dX/dt = (rho - d*r^2) X  +  a * grad(I_T)(X)
                             +  b * ( X x grad(I_T)(X) )
                             +  c * (1/4) grad(J)(X)

with

    grad(I_T)       = (y z, z x, x y)                              [degree 2]
    X x grad(I_T)   = ( x(y^2 - z^2), y(z^2 - x^2), z(x^2 - y^2) ) [degree 3]
    (1/4) grad(J)   = (x^3, y^3, z^3)                              [degree 3]

so component-wise:

    dx/dt = (rho - d*r^2) x + a*y z + b*x(y^2 - z^2) + c*x^3
    dy/dt = (rho - d*r^2) y + a*z x + b*y(z^2 - x^2) + c*y^3
    dz/dt = (rho - d*r^2) z + a*x y + b*z(x^2 - y^2) + c*z^3

These five terms are exactly the independent T-equivariant maps of degree <= 3:
the scalar rho*X (forced by irreducibility), the radial cubic r^2*X, the unique
degree-2 equivariant grad(I_T), the degree-3 cross term X x grad(I_T), and the
anisotropic cubic (x^3,y^3,z^3) = (1/4)grad(J). Every term is checked numerically
in tests/test_t_equivariance.py.

The grad(I_T) term is the T-signature: I_T = xyz is anti-invariant under the
octahedral 4-fold c4z, so the a-term breaks O while preserving T, and (being even
under inversion) it also breaks T_h -- the attractor carries the proper group T
exactly. The radial, cross, and (x^3,y^3,z^3) terms are all O-equivariant and
odd under inversion, so they neither restore O nor break chirality on their own;
setting a = b = 0 restores the full octahedral symmetry (the load-bearing test).

THE IRREDUCIBILITY OBSTRUCTION (the central finding of the polyhedral phase).
Because T acts irreducibly on R^3, the only T-commuting linear map is rho*I:
there is no saddle-focus linear engine. An exhaustive search of this family
(degree-3 here, plus degree-4/5 equivariants and non-gradient rotational terms;
random and fixed-point-stability-targeted; ~10^5 evaluations) finds NO dt-robust
strange attractor: bounded orbits relax to
equilibria or limit cycles. Trigonometric "labyrinth" T-equivariant flows
(Thomas-style) yield at most weak (lambda1 ~ 0.04), spatially localized chaos
that breaks the symmetry (a single labyrinth cell, not a T-invariant set). This
is strong numerical evidence -- not a proof -- that low-degree equivariant 3-D
flows under an *irreducible* polyhedral action resist chaos, in contrast to the
*reducible* C_n / S4 actions (atlas/family.py, atlas/family_s4.py), which host
chaos readily via their saddle-focus engine. It is a candidate explanation for
the absence of any polyhedral chaotic FLOW in the literature (Phase 0). See
docs/SYMMETRY_PROGRAM.md and the chaos-vs-symmetry-group law (Phase 3).

Where the evidence lives: docs/METHODOLOGY.md section 10 (10.2 polynomial
campaigns, 10.3 stability-targeted scan, 10.4 labyrinth route) records the
campaigns and their outcomes, and the T/O/I rows of scripts/symmetry_law.py
(results/symmetry_law.csv) are the committed, reproducible matched-budget runs.
The higher-degree (degree-4/5) and non-gradient campaigns were scratch runs and
their scripts were not committed; only this degree-3 family is in the tree.

CAUTION for any further search: the batched two-orbit estimate at coarse dt is
optimistic *to the point of fabricating chaos around stable equilibria* (a
dt=0.01 RK4 artifact that reverses sign by dt=0.0025). Always re-certify at finer
dt with the long integrator and keep only dt-robust positive exponents.

The evaluator mirrors atlas/family.py and atlas/family_s4.py: a vectorized RK4
over a whole batch, a two-orbit Benettin estimate of the largest Lyapunov
exponent, and a decimated orbit sample for fingerprinting.
"""

from __future__ import annotations

import numpy as np

# Parameter order: (rho, d, a, b, c)
_PARAM_NAMES = ("rho", "d", "a", "b", "c")
PARAM_LO = np.array([1.0, 0.10, 1.0, -11.0, -4.0])
PARAM_HI = np.array([12.0, 2.00, 14.0, 4.0, 2.0])
N_PARAMS = 5

ESCAPE = 60.0

# Generators of T (order 12), exported for verify scripts / symmetry residuals.
SIGMA3 = np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]])  # (x,y,z)->(y,z,x)
SIGMA2 = np.array([[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]])  # (x,y,z)->(x,-y,-z)


def param_names():
    """Names of the five parameters, in the order PARAM_LO/PARAM_HI use."""
    return list(_PARAM_NAMES)


def group_elements():
    """The 12 rotation matrices of the tetrahedral group T (closure of the gens)."""
    elems = [np.eye(3)]
    gens = [SIGMA3, SIGMA2]

    def known(M):
        return any(np.allclose(M, E, atol=1e-9) for E in elems)

    changed = True
    while changed:
        changed = False
        for A in list(elems):
            for g in gens:
                M = A @ g
                if not known(M):
                    elems.append(M)
                    changed = True
    return elems


def velocity(state, p):
    """Right-hand side (dx, dy, dz) at a single real state X = (x, y, z).

    ``p`` is a length-4 array-like in the order given by ``param_names()``.
    """
    x, y, z = (float(v) for v in state)
    rho, d, a, b, c = (float(v) for v in p)
    r2 = x * x + y * y + z * z
    rad = rho - d * r2
    return np.array(
        [
            rad * x + a * y * z + b * x * (y * y - z * z) + c * x**3,
            rad * y + a * z * x + b * y * (z * z - x * x) + c * y**3,
            rad * z + a * x * y + b * z * (x * x - y * y) + c * z**3,
        ]
    )


def _deriv(X, Y, Z, p):
    rho, d, a, b, c = p.T
    r2 = X * X + Y * Y + Z * Z
    rad = rho - d * r2
    dX = rad * X + a * Y * Z + b * X * (Y * Y - Z * Z) + c * X**3
    dY = rad * Y + a * Z * X + b * Y * (Z * Z - X * X) + c * Y**3
    dZ = rad * Z + a * X * Y + b * Z * (X * X - Y * Y) + c * Z**3
    return dX, dY, dZ


def _rk4(X, Y, Z, p, dt):
    k1x, k1y, k1z = _deriv(X, Y, Z, p)
    k2x, k2y, k2z = _deriv(X + 0.5 * dt * k1x, Y + 0.5 * dt * k1y, Z + 0.5 * dt * k1z, p)
    k3x, k3y, k3z = _deriv(X + 0.5 * dt * k2x, Y + 0.5 * dt * k2y, Z + 0.5 * dt * k2z, p)
    k4x, k4y, k4z = _deriv(X + dt * k3x, Y + dt * k3y, Z + dt * k3z, p)
    return (
        X + dt / 6.0 * (k1x + 2 * k2x + 2 * k3x + k4x),
        Y + dt / 6.0 * (k1y + 2 * k2y + 2 * k3y + k4y),
        Z + dt / 6.0 * (k1z + 2 * k2z + 2 * k3z + k4z),
    )


def evaluate_batch(
    params,
    dt=0.005,
    transient=16_000,
    lyap_steps=50_000,
    sample_steps=24_000,
    sample_every=16,
    seed=0,
):
    """Evaluate a (B, 5) batch of T parameter vectors.

    Returns dict with:
      lle     (B,) largest Lyapunov exponent (nan where dead)
      alive   (B,) bool: bounded and finite throughout
      points  (B, sample_steps // sample_every, 3) orbit samples (nan where dead)
    """
    p = np.asarray(params, dtype=float)
    B = p.shape[0]
    rng = np.random.default_rng(seed)

    # Break the symmetry of the initial condition so a symmetric (origin-bound)
    # transient does not bias toward the unstable fixed point at 0.
    X = rng.uniform(0.1, 0.6, B)
    Y = rng.uniform(-0.4, 0.4, B)
    Z = rng.uniform(-0.4, 0.4, B)
    alive = np.ones(B, dtype=bool)

    def cull():
        nonlocal alive
        bad = (
            ~np.isfinite(X) | (np.abs(X) > ESCAPE)
            | ~np.isfinite(Y) | (np.abs(Y) > ESCAPE)
            | ~np.isfinite(Z) | (np.abs(Z) > ESCAPE)
        )
        alive &= ~bad
        X[bad] = 0.0
        Y[bad] = 0.0
        Z[bad] = 0.0

    with np.errstate(all="ignore"):
        for i in range(transient):
            X, Y, Z = _rk4(X, Y, Z, p, dt)
            if i % 200 == 0:
                cull()
        cull()

        d0 = 1e-8
        Xs, Ys, Zs = X + d0, Y.copy(), Z.copy()
        log_sum = np.zeros(B)
        for i in range(lyap_steps):
            X, Y, Z = _rk4(X, Y, Z, p, dt)
            Xs, Ys, Zs = _rk4(Xs, Ys, Zs, p, dt)
            dX, dY, dZ = Xs - X, Ys - Y, Zs - Z
            dist = np.sqrt(dX * dX + dY * dY + dZ * dZ)
            dist = np.where(dist > 0, dist, 1e-16)
            log_sum += np.log(dist / d0)
            scale = d0 / dist
            Xs = X + dX * scale
            Ys = Y + dY * scale
            Zs = Z + dZ * scale
            if i % 200 == 0:
                cull()
        cull()
        lle = np.where(alive, log_sum / (lyap_steps * dt), np.nan)

        n_keep = sample_steps // sample_every
        pts = np.full((B, n_keep, 3), np.nan)
        k = 0
        for i in range(sample_steps):
            X, Y, Z = _rk4(X, Y, Z, p, dt)
            if i % sample_every == 0:
                pts[:, k, 0] = X
                pts[:, k, 1] = Y
                pts[:, k, 2] = Z
                k += 1
            if i % 200 == 0:
                cull()
        cull()
        pts[~alive] = np.nan

    return {"lle": lle, "alive": alive, "points": pts}


def sample_params(B, rng):
    """Uniform random parameter vectors inside the search box."""
    return rng.uniform(PARAM_LO, PARAM_HI, size=(B, N_PARAMS))


def mutate_params(parent, rng, sigma_frac=0.08):
    """Gaussian mutation, clipped to the box; sigma_frac is a fraction of each range."""
    span = PARAM_HI - PARAM_LO
    child = parent + rng.normal(0, sigma_frac, size=parent.shape) * span
    return np.clip(child, PARAM_LO, PARAM_HI)
