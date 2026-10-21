"""Lyapunov spectrum and Kaplan-Yorke dimension of the Trinacria attractor.

Standard Benettin tangent-space method: an orthonormal frame is advected by the
linearized flow and re-orthonormalized by QR at every step. The optional
segment split returns the spectrum of each of ``segments`` equal sub-windows
so the caller can form a segment standard error (std / sqrt(segments)), as in
scripts/convergence_study.py. lyapunov_spectrum_batch runs the same scheme on
an ensemble of initial conditions at once (used by the dt-halving arm of
scripts/verify_trinacria.py).
"""

from __future__ import annotations

import numpy as np

from .dynamics import (
    ALPHA, BETA, DELTA, EPS, GAMMA, LAM, NU, OMEGA, DEFAULT_STATE,
    jacobian, rk4_step, velocity,
)

PARAMS = (ALPHA, BETA, OMEGA, GAMMA, DELTA, NU, LAM, EPS)


def _tangent_rk4(s, Q, dt, p):
    def f(s):
        return velocity(s, *p)

    def J(s):
        return jacobian(s, *p)

    k1 = f(s); j1 = J(s) @ Q
    s2 = s + 0.5 * dt * k1
    k2 = f(s2); j2 = J(s2) @ (Q + 0.5 * dt * j1)
    s3 = s + 0.5 * dt * k2
    k3 = f(s3); j3 = J(s3) @ (Q + 0.5 * dt * j2)
    s4 = s + dt * k3
    k4 = f(s4); j4 = J(s4) @ (Q + dt * j3)
    return (
        s + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4),
        Q + (dt / 6.0) * (j1 + 2 * j2 + 2 * j3 + j4),
    )


def lyapunov_spectrum(n_steps=600_000, dt=0.005, transient=20_000,
                      state0=DEFAULT_STATE, params=PARAMS, segments=None):
    """Return the three Lyapunov exponents, largest first.

    With ``segments = k`` (k >= 2; the sub-windows are as equal as integer
    step counts allow) also return a (k, 3) array holding the spectrum of
    each sub-window, as ``(spectrum, segment_spectra)``.
    """
    s = np.asarray(state0, dtype=float)
    for _ in range(transient):
        s = rk4_step(s, dt, *params)
    Q = np.eye(3)
    sums = np.zeros(3)
    if segments:
        bounds = np.linspace(0, n_steps, int(segments) + 1).astype(int)
        seg = np.zeros((int(segments), 3))
        k, seg_start, seg_sum = 0, 0, np.zeros(3)
    for i in range(n_steps):
        s, Q = _tangent_rk4(s, Q, dt, params)
        Q, R = np.linalg.qr(Q)
        inc = np.log(np.abs(np.diag(R)))
        sums += inc
        if segments:
            seg_sum += inc
            if i + 1 == bounds[k + 1]:
                seg[k] = seg_sum / ((i + 1 - seg_start) * dt)
                k, seg_start, seg_sum = k + 1, i + 1, np.zeros(3)
    spectrum = sums / (n_steps * dt)
    if segments:
        return spectrum, seg
    return spectrum



def _velocity_rows(S, p):
    al, be, om, ga, de, nu, lam, ep = p
    x, y, z = S[:, 0], S[:, 1], S[:, 2]
    c = al * z * z - be
    d2 = x * x - y * y
    return np.stack([
        c * x - om * z * y + ga * d2 + 2 * de * z * x * y,
        c * y + om * z * x - 2 * ga * x * y + de * z * d2,
        nu * z - z * z * z - lam * z * (x * x + y * y) + ep * (3 * x * x * y - y * y * y),
    ], axis=1)


def _jacobian_rows(S, p):
    al, be, om, ga, de, nu, lam, ep = p
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


def lyapunov_spectrum_batch(states0, n_steps, dt=0.005, transient=20_000,
                            params=PARAMS, segments=None):
    """Spectra of a batch of initial conditions, integrated side by side.

    The same tangent RK4 + QR-every-step scheme as lyapunov_spectrum,
    vectorized over the rows of ``states0`` (B, 3), so an ensemble of
    independent orbits costs about as much as one pure-Python orbit.
    Returns a (B, 3) array, or ``(spectra, segment_spectra)`` with a
    (B, k, 3) segment array when ``segments = k``. Rows that leave
    |X| < 50 become NaN.
    """
    S = np.array(states0, dtype=float)
    B = len(S)
    p = tuple(float(v) for v in params)

    def f(X):
        return _velocity_rows(X, p)

    for _ in range(transient):
        k1 = f(S)
        k2 = f(S + 0.5 * dt * k1)
        k3 = f(S + 0.5 * dt * k2)
        k4 = f(S + dt * k3)
        S = S + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    Q = np.tile(np.eye(3), (B, 1, 1))
    sums = np.zeros((B, 3))
    if segments:
        bounds = np.linspace(0, n_steps, int(segments) + 1).astype(int)
        seg = np.zeros((B, int(segments), 3))
        k, seg_start, seg_sum = 0, 0, np.zeros((B, 3))
    with np.errstate(all="ignore"):
        for i in range(n_steps):
            k1 = f(S); j1 = _jacobian_rows(S, p) @ Q
            S2 = S + 0.5 * dt * k1
            k2 = f(S2); j2 = _jacobian_rows(S2, p) @ (Q + 0.5 * dt * j1)
            S3 = S + 0.5 * dt * k2
            k3 = f(S3); j3 = _jacobian_rows(S3, p) @ (Q + 0.5 * dt * j2)
            S4 = S + dt * k3
            k4 = f(S4); j4 = _jacobian_rows(S4, p) @ (Q + dt * j3)
            S = S + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
            Q = Q + (dt / 6.0) * (j1 + 2 * j2 + 2 * j3 + j4)
            Q, R = np.linalg.qr(Q)
            inc = np.log(np.abs(np.diagonal(R, axis1=1, axis2=2)))
            sums += inc
            if segments:
                seg_sum += inc
                if i + 1 == bounds[k + 1]:
                    seg[:, k] = seg_sum / ((i + 1 - seg_start) * dt)
                    k, seg_start, seg_sum = k + 1, i + 1, np.zeros((B, 3))
    spectra = sums / (n_steps * dt)
    bad = ~np.all(np.isfinite(S), axis=1) | (np.abs(S).max(axis=1) > 50.0)
    spectra[bad] = np.nan
    if segments:
        seg[bad] = np.nan
        return spectra, seg
    return spectra

def kaplan_yorke_dimension(spectrum):
    """Kaplan-Yorke dimension from a spectrum."""
    exps = np.sort(np.asarray(spectrum))[::-1]
    cum = 0.0
    for j, lam in enumerate(exps):
        if cum + lam < 0.0:
            return j + cum / abs(lam) if j else 0.0
        cum += lam
    return float(len(exps))
