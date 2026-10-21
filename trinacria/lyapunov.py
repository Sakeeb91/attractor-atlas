"""Lyapunov spectrum and Kaplan-Yorke dimension of the Trinacria attractor.

Standard Benettin tangent-space method: an orthonormal frame is advected by the
linearized flow and re-orthonormalized by QR at every step. The optional
segment split returns the spectrum of each of ``segments`` equal sub-windows
so the caller can form a segment standard error (std / sqrt(segments)), as in
scripts/convergence_study.py.
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


def kaplan_yorke_dimension(spectrum):
    """Kaplan-Yorke dimension from a spectrum."""
    exps = np.sort(np.asarray(spectrum))[::-1]
    cum = 0.0
    for j, lam in enumerate(exps):
        if cum + lam < 0.0:
            return j + cum / abs(lam) if j else 0.0
        cum += lam
    return float(len(exps))
