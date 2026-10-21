"""The trinacria package must reproduce the equivariance-verified D3 family flow.

tests/test_d_equivariance.py proves that atlas.family_d.velocity is
D_n-equivariant. The trinacria package re-expresses the n = 3 member in real
(x, y, z) coordinates with the canonical parameters baked in. If the two agree
at random states, the real-coordinate specimen inherits the proven
equivariance; if the hand-derivation of conj(w)^2, i*delta*z*conj(w)^2 or
Im(w^3) is wrong, this fails. The package-level equivariance and equilibrium
checks below are independent of the family module.
"""

from __future__ import annotations

import numpy as np
import pytest

from atlas.family_d import group_matrices as fam_group_matrices
from atlas.family_d import velocity as fam_velocity
from trinacria.dynamics import (
    ALPHA, BETA, DELTA, EPS, GAMMA, LAM, NU, OMEGA, divergence, group_matrices,
    jacobian, velocity,
)
from trinacria.equilibria import classify, group_orbits, rotation_sense
from trinacria.lyapunov import kaplan_yorke_dimension, lyapunov_spectrum

# family_d parameter order: (alpha, beta, omega, gamma, delta, nu, lam, eps)
CANON = np.array([ALPHA, BETA, OMEGA, GAMMA, DELTA, NU, LAM, EPS])


def test_canonical_parameters_are_nominee_8():
    assert CANON.tolist() == [1.4773, 1.4668, 0.7168, 0.8446, 1.0567, 2.3978, 1.2975, -1.5491]


def test_real_form_matches_complex_family():
    rng = np.random.default_rng(2026)
    for _ in range(300):
        x, y, z = rng.uniform(-2, 2, 3)
        dw, dz = fam_velocity(complex(x, y), float(z), CANON, 3)
        vx, vy, vz = velocity(np.array([x, y, z]))
        assert vx == pytest.approx(dw.real, abs=1e-12, rel=0)
        assert vy == pytest.approx(dw.imag, abs=1e-12, rel=0)
        assert vz == pytest.approx(float(np.real(dz)), abs=1e-12, rel=0)


def test_jacobian_matches_finite_difference():
    rng = np.random.default_rng(11)
    h = 1e-6
    for _ in range(40):
        s = rng.uniform(-1.5, 1.5, 3)
        J = jacobian(s)
        for k in range(3):
            ds = np.zeros(3); ds[k] = h
            fd = (velocity(s + ds) - velocity(s - ds)) / (2 * h)
            assert np.allclose(J[:, k], fd, atol=1e-5)


def test_divergence_matches_jacobian_trace():
    rng = np.random.default_rng(3)
    for _ in range(40):
        s = rng.uniform(-2, 2, 3)
        assert divergence(s) == pytest.approx(np.trace(jacobian(s)), abs=1e-10)


def test_group_matrices_are_d3():
    G = group_matrices()
    assert len(G) == 6
    for g in G:
        assert np.linalg.det(g) == pytest.approx(1.0, abs=1e-12)
        assert np.allclose(g @ g.T, np.eye(3), atol=1e-12)
    # closed under multiplication, and the same set as the family's D_3
    for a in G:
        for b in G:
            assert any(np.allclose(a @ b, c, atol=1e-12) for c in G)
    fam = fam_group_matrices(3)
    assert all(any(np.allclose(g, f, atol=1e-12) for f in fam) for g in G)
    # three rotations about z (trace 1 + 2cos) and three half-turns (trace -1)
    traces = sorted(round(float(np.trace(g)), 9) for g in G)
    assert traces == [-1.0, -1.0, -1.0, 0.0, 0.0, 3.0]
    # kappa is the half-turn about x, not a mirror
    assert np.allclose(G[3], np.diag([1.0, -1.0, -1.0]))


def test_package_velocity_is_d3_equivariant():
    rng = np.random.default_rng(7)
    G = group_matrices()
    for _ in range(200):
        s = rng.uniform(-2, 2, 3)
        f = velocity(s)
        for g in G:
            assert np.max(np.abs(velocity(g @ s) - g @ f)) < 1e-12


def test_package_jacobian_is_d3_conjugate():
    # J(gX) = g J(X) g^T follows from equivariance; checks the hand-derived J.
    rng = np.random.default_rng(8)
    G = group_matrices()
    for _ in range(50):
        s = rng.uniform(-2, 2, 3)
        J = jacobian(s)
        for g in G:
            assert np.max(np.abs(jacobian(g @ s) - g @ J @ g.T)) < 1e-11


@pytest.mark.parametrize(
    "name,m",
    [
        ("z-flip", np.diag([1.0, 1.0, -1.0])),
        ("y-mirror", np.diag([1.0, -1.0, 1.0])),
        ("inversion", -np.eye(3)),
        ("C6", np.array([[0.5, -np.sqrt(3) / 2, 0.0], [np.sqrt(3) / 2, 0.5, 0.0], [0.0, 0.0, 1.0]])),
    ],
)
def test_canonical_flow_breaks_larger_groups(name, m):
    rng = np.random.default_rng(9)
    states = rng.uniform(-1.5, 1.5, (40, 3))
    viol = [np.max(np.abs(velocity(m @ s) - m @ velocity(s))) > 1e-6 for s in states]
    assert all(viol), name


def test_on_axis_equilibria_are_exact_roots():
    r = np.sqrt(NU)
    for eq in ([0.0, 0.0, 0.0], [0.0, 0.0, r], [0.0, 0.0, -r]):
        assert np.max(np.abs(velocity(np.array(eq)))) < 1e-12


def test_origin_is_a_saddle_with_double_eigenvalue_minus_beta():
    eigs, label = classify(np.zeros(3))
    assert label == "saddle"
    assert np.allclose(np.sort(eigs.real), [-BETA, -BETA, NU], atol=1e-12)
    assert np.max(np.abs(eigs.imag)) < 1e-12


def test_on_axis_pair_counter_rotates():
    r = np.sqrt(NU)
    up, down = np.array([0.0, 0.0, r]), np.array([0.0, 0.0, -r])
    orbits = group_orbits([np.zeros(3), up, down])
    sizes = sorted(o["size"] for o in orbits)
    assert sizes == [1, 2]
    assert all(o["complete"] and o["spectrum_gap"] < 1e-9 for o in orbits)
    assert rotation_sense(up) == pytest.approx(-2 * OMEGA * r, abs=1e-12)
    assert rotation_sense(down) == pytest.approx(2 * OMEGA * r, abs=1e-12)


def test_lyapunov_segments_average_to_spectrum():
    spec, seg = lyapunov_spectrum(n_steps=2000, dt=0.005, transient=500, segments=4)
    assert seg.shape == (4, 3)
    assert np.allclose(seg.mean(axis=0), spec, atol=1e-12)
    assert kaplan_yorke_dimension([0.3, 0.0, -1.5]) == pytest.approx(2.2)
