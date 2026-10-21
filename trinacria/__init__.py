"""The Trinacria attractor: the D3 dihedral member, two counter-rotating halves."""

from .dynamics import (
    ALPHA, BETA, DELTA, EPS, GAMMA, LAM, NU, OMEGA, DEFAULT_STATE,
    divergence, group_matrices, jacobian, rk4_step, trajectory, velocity,
)

__all__ = [
    "ALPHA", "BETA", "OMEGA", "GAMMA", "DELTA", "NU", "LAM", "EPS",
    "DEFAULT_STATE", "divergence", "group_matrices", "jacobian", "rk4_step",
    "trajectory", "velocity",
]
