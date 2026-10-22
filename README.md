# The Attractor Atlas

*Five novel strange attractors, together with the quality-diversity engine that finds them: **Aurelia** (C₃), **Naiad** (C₂), **Cassiopea** (C₄) and **Mobula** (S₄ rotoreflection), discovered, certified and named in June 2026, and **Trinacria** (D₃, the chiral dihedral group), added in October 2026. They are the first entries in a [periodic table of strange attractors ordered by symmetry](docs/SYMMETRY_PROGRAM.md); Aurelia, the moon jellyfish of the family, leads the tour.*

**[Explore the first four live in your browser](https://sakeeb91.github.io/attractor-atlas/)**: a WebGL viewer with 400,000 orbit points, 3,000 particles riding the flow in real time, orbit controls, live parameter sliders, and a switcher between attractors. Trinacria is not in the viewer yet.

**[Read the story](STORY.md)**: a long-form narrative of the whole arc, from the failed first family through the 600-universe search and the certification of Aurelia, to the quality-diversity engine that found its siblings Naiad and Cassiopea, and on to Mobula — the first member whose symmetry is a rotation-plus-reflection.

**[Read the methods](docs/METHODOLOGY.md)**: the complete methodology — every construction and formula (the equivariant families, term audits, invariant theory), the search and certification protocols, the step-size artifact and its gate, the polyhedral obstruction evidence, and the symmetry-law protocol, with a claim-by-claim reproducibility map.

![The Aurelia attractor](gallery/aurelia_hero.png)

## The system

The Aurelia attractor lives in a three-dimensional autonomous flow with two quadratic
nonlinearities in the plane and a cubic one on the axis:

```
dx/dt = a(z − b)x − cy + b(x² − y²)
dy/dt = cx + a(z − b)y − 2bxy
dz/dt = a(1 + z) − z³ − c(x² + y²)
```

with the canonical parameters

```
a = 1.4,   b = 0.5,   c = 1.9
```

and any generic initial condition near the origin, e.g. `(0.3, 0.2, 0.1)`.

Writing the planar coordinates as a single complex number **w = x + iy**, the system
collapses to two elegant equations:

```
dw/dt = [a(z − b) + ic] w + b w̄²
dz/dt = a(1 + z) − z³ − c|w|²
```

The conjugate-squared term `w̄²` is the heart of the design: under a rotation of the
plane by 120°, `w → e^(2πi/3) w` maps solutions to solutions, so the flow is
**equivariant under the cyclic group C₃**. The strange attractor inherits this
three-fold symmetry — visible as the triskelion in the crown view below. The `z`
equation couples back through the radius `|w|²` alone, preserving the symmetry while
pumping trajectories vertically.

## Why "Aurelia"

*Aurelia aurita* is the moon jellyfish. Seen from the side, the attractor is a
translucent bell with veils trailing beneath it and a luminous column rising through
its center — the unstable axis along which trajectories are reinjected. From the
Latin *aureus*, "golden": orbits spend their slow hours in the gilded crown at the
top of the bell before being flung back into the dark. The name was unclaimed in the
dynamical-systems literature at the time of discovery.

## Gallery

| The crown (top view) | The bell (side view) |
|---|---|
| ![Crown view](gallery/aurelia_crown.png) | ![Bell view](gallery/aurelia_bell.png) |

![Contact sheet](gallery/aurelia_sheet.png)

One full orbit of the camera (`python scripts/render_animation.py`):

![Rotating view](gallery/aurelia_rotation.gif)

## Certification of chaos

All numbers below are reproducible with `python scripts/verify_chaos.py`
(recorded in [`results/verification.json`](results/verification.json)).

| Property | Value |
|---|---|
| Lyapunov spectrum (Benettin/QR, 600k steps) | **λ₁ = +0.2327**, λ₂ = +0.0013 ≈ 0, λ₃ = −1.5119 |
| Kaplan–Yorke dimension | **D_KY ≈ 2.155** (fractal) |
| λ₁ across 5 random initial conditions | 0.195 … 0.252 (always positive) |
| Convergence study (3 dt × 3 lengths × 10 ICs) | λ₁ = +0.230 ± 0.009, D_KY = 2.151 ± 0.006 ([gallery/convergence.png](gallery/convergence.png), [results/convergence.json](results/convergence.json)) |
| Equilibria | exactly one: saddle-focus at (0, 0, 1.5229) |
| Eigenvalues at the equilibrium | 1.432 ± 1.9i (unstable spiral), −5.558 (strong contraction) |
| Boundedness | confirmed over 4,000,000 RK4 steps |
| Time-averaged divergence | −1.290 (dissipative: volumes contract onto the attractor) |
| Attractor extent | x ∈ [−1.40, 1.60], y ∈ [−1.67, 1.27], z ∈ [−0.98, 1.50] |

The signature `(+, 0, −)` spectrum with a fractal dimension strictly between 2 and 3
certifies a **strange attractor**, not a limit cycle or a transient.

### One equilibrium to rule them all

The flow has a *single* fixed point — a saddle-focus on the symmetry axis with
Shilnikov-type geometry. Trajectories spiral outward around the axis (eigenvalues
1.432 ± 1.9i), spill over the rim of the bell, fall along the outer veils, and are
recaptured by the strongly contracting direction (−5.558) at the base of the central
column. Stretch, fold, repeat — forever, never twice the same way.

### Route to chaos and limit cycles

Sweeping the rotation rate `c` with `a, b` fixed
(`python scripts/bifurcation_scan.py`) shows the chaotic sea and its periodic
windows: where the largest Lyapunov exponent falls to ≈ 0, the attractor collapses
to a stable limit cycle; where it is positive, chaos reigns. The canonical
`c = 1.9` sits inside a robust chaotic band.

![Bifurcation scan](gallery/bifurcation.png)

### The homoclinic hunt

Shilnikov's theorem needs two things: the eigenvalue inequality (satisfied by
all three systems, saddle index ν = 0.19–0.26 < 1) and a homoclinic orbit —
the unstable manifold returning exactly to the equilibrium. A shooting hunt
(`python scripts/homoclinic_hunt.py`) launches the manifold from the unstable
eigenplane and sweeps the rotation rate. **Naiad**'s manifold returns to
within **7.5×10⁻⁵** of its saddle-focus at w_rot = 1.833, with a ladder of
seven near-connection windows accumulating below the canonical parameter
(where the distance is still 7.6×10⁻⁴). **Cassiopea**'s returns to within
6.5×10⁻⁴. For those two, homoclinic (Shilnikov) bifurcations almost certainly
thread the family. **Aurelia** is the holdout: along the entire c-band its
manifold keeps a standoff distance of ≈0.13 — the reinjection column re-enters
the spiral plane at finite radius, and any connection must lie off the c-axis.

![Homoclinic hunt](gallery/homoclinic.png)

## Novelty

Searched June 10, 2026, before publication:

* **The functional form** — a 3-D autonomous flow built from the C₃-equivariant
  complex coupling `w̄²` (i.e. the pair `x² − y²`, `−2xy`) joined to cubic `z`
  dynamics through `|w|²` — does not appear in the known attractor catalogues
  (Lorenz, Rössler, Thomas, Aizawa/Langford, Halvorsen, Dadras, Chen, Sprott A–S,
  Rabinovich–Fabrikant, multi-scroll/multi-wing families) or in web/literature
  searches for the equation terms. (Full positioning pass against the
  equivariant-covers, multi-scroll, and normal-form literature:
  [docs/RELATED_WORK.md](docs/RELATED_WORK.md).)
* **Closest relatives, and how this differs:** Field & Golubitsky's celebrated
  *symmetric chaos icons* use equivariant polynomial couplings of the same spirit —
  but they are 2-D discrete **maps**, not continuous 3-D flows. The Aizawa/Langford
  attractor is a 3-D flow with rotational coupling — but it is SO(2)-symmetric
  (continuous rotation), has no equivariant `w̄²` term, and has a qualitatively
  different two-equilibrium structure. Discrete Cₙ rotation symmetry *does* occur
  in published 3-D flows — the proto-Lorenz n-fold covers of Miranda & Stone, the
  cover-system theory of Gilmore & Letellier, Thomas's cyclically symmetric system,
  and the engineered multi-scroll families — but in all of these the symmetry is
  either imposed by lifting known dynamics through an angle-multiplying map
  (leaving the vector field singular on the axis and the dynamics locally identical
  to the image system) or acts without Aurelia's plane-plus-axis anatomy. The
  Aurelia system occupies a sparsely populated intersection: *globally polynomial
  flow* × *discrete C₃ equivariance by construction, not by covering* × *single
  organizing saddle-focus*. Full survey:
  [docs/RELATED_WORK.md](docs/RELATED_WORK.md).
* **The name** "Aurelia attractor" had no prior usage.
* **Measured novelty.** Beyond the literature search, novelty is now a number:
  reducing attractors to shape fingerprints (D2 distance histograms + PCA shape
  ratios, see `atlas/fingerprint.py`) and comparing against all 135 known
  chaotic systems in the [dysts](https://github.com/williamgilpin/dysts)
  catalog, Aurelia's nearest known neighbor is the Aizawa attractor at
  distance 0.185. The catalog's own median nearest-neighbor spacing is 0.115;
  83% of known systems sit closer to something known than Aurelia does. By
  this metric, Aurelia is more isolated in shape space than five-sixths of
  the established bestiary, and its nearest relative is exactly the system
  the qualitative analysis predicted. The metric itself is validated in
  `scripts/validate_fingerprint.py`: documented Lorenz-family variants
  cluster (positive control), the algebraically defined Sprott family does
  not (negative control — the fingerprint reads geometry, not algebra), and
  the within-system noise floor is ≈0.05, putting all three of our systems
  3–5 noise floors from anything known. Under re-simulation Aurelia's
  distance is stable at 0.184, with the InteriorSquirmer system
  statistically tied with Aizawa — at this isolation, ties between distant
  neighbors are expected.
  Map: [gallery/fingerprint_map.png](gallery/fingerprint_map.png) ·
  details: [results/fingerprint_validation.json](results/fingerprint_validation.json).

## Discovery method

The attractor was found by Monte-Carlo search over a seven-parameter
C₃-equivariant family (`python scripts/search_parameters.py`), filtering for
bounded orbits with a strongly positive largest Lyapunov exponent, then judged by
eye for beauty — a workflow whose ancestor is Sprott's automated aesthetic hunt
over 2-D maps in 1993. The winning region of parameter space collapsed — remarkably — to
just **three distinct constants** (a, b, c) = (1.4, 0.5, 1.9) with the largest
Lyapunov exponent essentially unchanged.

## The atlas engine: searching for the next Aurelia

`atlas/` is a quality-diversity search engine that hunts for *new* chaotic
attractors in the C&#8345;-equivariant generalization of Aurelia's design
(rotation orders n = 2 to 5, so two-, four-, and five-fold flowers are in
reach). Instead of optimizing a single score, it runs MAP-Elites: descriptor
space (rotation order × Lyapunov exponent × vertical aspect) is tiled into
cells, and each cell keeps the candidate most *novel* relative to a reference
set of known attractors.

Novelty is a measured number, not a literature claim: every candidate's orbit
is reduced to a shape fingerprint (D2 pairwise-distance histogram + PCA shape
ratios + fill factor), and its novelty is the minimum fingerprint distance to
the [dysts](https://github.com/williamgilpin/dysts) catalog of known chaotic
systems plus everything already discovered in the archive.

Evolved chaotic flows exist in the literature — genetic programming over
Lorenz-like vector fields (Pan & Das 2015), and most recently the 20,000-system
Panda corpus bred from the dysts catalog as training data (Lai, Bao & Gilpin
2025) — but, to our knowledge, this is the first *quality-diversity* search over
a family of dynamical systems, with novelty measured against a catalog of known
attractors rather than chaos strength optimized
([docs/RELATED_WORK.md](docs/RELATED_WORK.md)).

The whole shape space at a glance — the dysts catalog as dots, the three
discoveries as stars on the rim of the known cloud (metric MDS of
fingerprint distances, Kruskal stress-1 = 0.13):

![Shape space of known attractors with the three new systems](gallery/fingerprint_map.png)

```bash
pip install dysts                                # reference catalog (one-time)
python scripts/build_catalog_fingerprints.py     # fingerprint known systems (one-time)
python scripts/run_atlas.py 12 96                # search: 12 generations, batches of 96
python scripts/render_atlas_montage.py           # render the most novel elites
```

The batched evaluator integrates an entire population simultaneously with
vectorized complex RK4 and reproduces Aurelia's Lyapunov exponent as a
sanity check. The archive persists across runs (`results/atlas/archive.json`),
so the search compounds.

## Naiad: the first attractor the atlas found

The first thing the search turned up worth keeping is **Naiad**, the two-fold
(C₂) sibling of Aurelia, named for the water nymph of fountains and springs.
Where Aurelia rotates three-fold, Naiad uses rotation order n = 2, so the
conjugate term is simply `conj(w)`: an anisotropic stretch that makes the
flow symmetric under a half-turn `(x, y) → (−x, −y)`. The attractor is a
flaring **fountain**: a wide turbulent bowl, a bright jet up the central
axis, a narrow stem at the base. Seen from above it is a lens, not a disk,
because two-fold symmetry stretches it along one direction.

![The Naiad attractor](gallery/naiad_hero.png)

| The lens (top) | The jet (side) |
|---|---|
| ![Naiad lens](gallery/naiad_lens.png) | ![Naiad jet](gallery/naiad_jet.png) |

```
dx/dt = (a(z − b) + g)x − w·y
dy/dt = w·x + (a(z − b) − g)y          a=1.2, b=0.68, w=3.7, g=1.7,
dz/dt = mu + nu·z − z³ − lam(x² + y²)   mu=1.2, nu=2.0, lam=2.1
```

Certified the same way as Aurelia (`python scripts/verify_naiad.py`,
recorded in [`results/naiad_verification.json`](results/naiad_verification.json)):

| Property | Value |
|---|---|
| Lyapunov spectrum | **λ₁ = +0.296**, λ₂ = +0.001 ≈ 0, λ₃ = −2.696 |
| Kaplan–Yorke dimension | **D_KY ≈ 2.110** |
| λ₁ across 5 random initial conditions | 0.284 … 0.341 (always positive) |
| Convergence study (3 dt × 3 lengths × 10 ICs) | λ₁ = +0.294 ± 0.009, D_KY = 2.110 ± 0.004 |
| Equilibria | exactly one: saddle-focus at (0, 0, 1.651) |
| Eigenvalues at the equilibrium | 1.166 ± 3.286i (unstable spiral), −6.18 |
| Time-averaged divergence | −2.386 (dissipative) |
| Novelty vs the dysts catalog | nearest known system at distance **0.29–0.34** across re-simulations (the rank-1 identity flips between SprottJerk and ForcedFitzHughNagumo at this isolation) |

That novelty distance of ≈0.3 makes Naiad the *most* isolated member of the
family in shape space (Aurelia sits at 0.18), against a catalog whose median
nearest-neighbor spacing is 0.115. Like Aurelia, it is organized by a single
Shilnikov-type saddle-focus, here with a faster spiral.

## Cassiopea: the four-fold member

**Cassiopea** completes a three-member family (C₂ Naiad, C₃ Aurelia, C₄
Cassiopea), named for *Cassiopea*, the upside-down jellyfish whose bell
carries a four-leaf-clover marking: a real animal with genuine four-fold
symmetry. The equivariant term is `conj(w)³`, the unique cubic surviving a
quarter-turn `w → iw`. Seen from above the attractor is a **four-armed
pinwheel star** with a white spiral core; from the side, a layered bell
with twin curtains.

![The Cassiopea attractor](gallery/cassiopea_hero.png)

| The star (top) | The bell (side) |
|---|---|
| ![Cassiopea star](gallery/cassiopea_star.png) | ![Cassiopea bell](gallery/cassiopea_bell.png) |

```
dx/dt = a(z − b)x − w·y + g(x³ − 3xy²)
dy/dt = w·x + a(z − b)y + g(y³ − 3x²y)     a=2.0, b=0.85, w=1.8, g=0.9,
dz/dt = mu + nu·z − z³ − lam(x² + y²)       mu=1.28, nu=1.8, lam=2.8
```

Certified via `python scripts/verify_cassiopea.py`
([`results/cassiopea_verification.json`](results/cassiopea_verification.json)):

| Property | Value |
|---|---|
| Lyapunov spectrum | **λ₁ = +0.525**, λ₂ ≈ 0, λ₃ = −2.412 |
| Kaplan–Yorke dimension | **D_KY ≈ 2.218** |
| λ₁ across 5 random initial conditions | 0.490 … 0.555 (always positive) |
| Convergence study (3 dt × 3 lengths × 10 ICs) | λ₁ = +0.523 ± 0.008, D_KY = 2.219 ± 0.006 |
| Equilibria | a central saddle-focus at (0, 0, 1.611) plus a C₄ quadruple of distant saddles |
| Eigenvalues at the saddle-focus | 1.522 ± 1.8i (unstable spiral), −5.984 |
| Time-averaged divergence | −1.868 (dissipative) |
| Novelty vs the dysts catalog | nearest known system NoseHoover at distance **0.17–0.22** across re-simulations (Aizawa statistically tied) |

Cassiopea carries the **strongest chaos and deepest fractal structure of the
family** (Aurelia λ₁ = 0.233, Naiad 0.296, Cassiopea 0.525), and it is the
first member with off-axis equilibria: four distant saddles at radius 3,
arranged in a square far below the attractor, with the familiar single
saddle-focus engine on the axis.

## Mobula: the rotoreflection member

The three jellyfish above are all symmetric under a pure **rotation**. **Mobula**
breaks that series: it is the first member symmetric under an *improper*
operation — a quarter-turn **combined with a vertical flip**,
`σ: (w, z) → (i·w, −z)`, and under neither the quarter-turn nor the flip alone.
This is the order-4 rotoreflection group **S₄** (Schoenflies), the point-group
symmetry of a devil ray (genus *Mobula*) breaching and somersaulting clear of
the water. From the side the attractor is a **broad-winged ray of layered
veils**; from above, a **four-fold pinwheel**. Two terms make the symmetry S₄
rather than C₄: `δ·z·conj(w)` and `ε·Im(w²)` — each allowed under the
rotoreflection but forbidden under the pure quarter-turn.

A targeted [literature strike](docs/RELATED_WORK.md) found no prior *constructed,
certified chaotic flow* whose attractor carries a rotoreflection symmetry that
is neither an inversion nor a reflection: the map literature (Reiter et al.) and
the order-2 inversion covers (Letellier–Gilmore) are the nearest prior art.

![The Mobula attractor](gallery/mobula_hero.png)

| The pinwheel (top) | The wing (side) |
|---|---|
| ![Mobula pinwheel](gallery/mobula_pinwheel.png) | ![Mobula wing](gallery/mobula_wing.png) |

```
dx/dt = (α·z² − β + δ·z)x − ω·y + γ(x³ − 3xy²)    α=1.4104, β=0.6701,
dy/dt = (α·z² − β − δ·z)y + ω·x + γ(y³ − 3x²y)     ω=2.4266, γ=0.3, δ=0.5534,
dz/dt = ν·z − z³ − λ·z(x² + y²) + 2ε·xy            ν=3.3, λ=2.1355, ε=1.2153
```

Equivariance is **unit-tested numerically** (`tests/`): the real-coordinate flow
above is checked to equal the complex S₄ family to 10⁻¹², and to *break* C₄ and
the pure flip. Certified via `python scripts/verify_mobula.py`
([`results/mobula_verification.json`](results/mobula_verification.json)):

| Property | Value |
|---|---|
| Lyapunov spectrum | **λ₁ = +0.330**, λ₂ ≈ 0, λ₃ = −1.238 |
| Kaplan–Yorke dimension | **D_KY ≈ 2.266** |
| λ₁ across 5 random initial conditions | 0.334 … 0.384 (always positive) |
| Equilibria | on-axis σ-conjugate pair (0, 0, ±1.817) with *identical* eigenvalues, plus four off-axis saddles forming one order-4 σ-orbit |
| Eigenvalues at (0,0,±1.817) | 3.984 ± 2.21i (unstable spiral), −6.6 |
| S₄ symmetry residual | 0.0136, versus 0.0455 for a pure-flip control (3.3× better) — genuine rotoreflection, not flip symmetry |
| Novelty vs the dysts catalog | nearest known system NuclearQuadrupole at distance 0.273 |

The family's chaos strength is **not monotonic in symmetry order**
(λ₁: C₂ 0.296, C₃ 0.233, C₄ 0.525, S₄ 0.330) — the systematic question the
[symmetry program](docs/SYMMETRY_PROGRAM.md) sets out to map. Because S₄ forbids
a constant in `dz/dt`, Mobula's on-axis equilibria are the *symmetric* triple
`z(ν − z²) = 0`, and the symmetry is visible in the fixed-point structure: σ
exchanges the two off-axis-z equilibria and they carry identical eigenvalues.

## Trinacria: the dihedral member

Every earlier member turns about a single vertical axis. **Trinacria** adds
half-turns about *horizontal* axes. Its symmetry group is the chiral dihedral
group **D₃** (order 6): the 120° rotation `ρ: (w, z) → (e^(2πi/3)·w, z)` about
the vertical axis, together with the half-turn `κ: (w, z) → (w̄, −z)`, that is
`(x, y, z) → (x, −y, −z)`, about the x axis. The six elements are three
rotations about z and three half-turns about horizontal axes at 0° and ±120°.
All six are proper rotations (determinant +1), so the attractor has no mirror
symmetry of any kind.

The half-turn changes how the flow can spin. A constant rotation rate is not
κ-equivariant, so the only admissible spin is odd in z: the term `iω·z·w` turns
the upper half one way and the lower half the other, and κ swaps the two. On
the attractor this shows as a **counter-rotating core**: within radius 0.6 of
the axis (about 55% of the time) the upper half winds counterclockwise and the
lower half clockwise at equal and opposite rates (mean angular velocity ±0.59
for r < 0.1, ±0.47 for 0.1 < r < 0.3, ±0.29 for 0.3 < r < 0.6). Farther out
each half turns weakly the other way, so the net angular momentum of each half
is close to zero.

*Trinacria* Heiberg 1863 is an extinct marine diatom genus (Paleocene to Eocene
fossils from Denmark), named for the triskelion of Sicily: three legs turning
one way, and turning the other way when the disc is flipped over. That reversal
is what κ does to this flow, and a diatom frustule is itself two valves stacked
top and bottom. Web searches for "Trinacria attractor" and for Trinacria as a
chaotic or dynamical system found no prior use, both at selection time and again
after certification (logged in [`results/d3_selection.json`](results/d3_selection.json)
and [`docs/RELATED_WORK.md`](docs/RELATED_WORK.md)).

![The Trinacria attractor](gallery/trinacria_hero.png)

| The triskelion (top) | The profile (side) |
|---|---|
| ![Trinacria triskelion](gallery/trinacria_triskelion.png) | ![Trinacria profile](gallery/trinacria_profile.png) |

```
dx/dt = (α·z² − β)x − ω·z·y + γ(x² − y²) + 2δ·z·xy      α=1.4773, β=1.4668,
dy/dt = (α·z² − β)y + ω·z·x − 2γ·xy + δ·z(x² − y²)       ω=0.7168, γ=0.8446,
dz/dt = ν·z − z³ − λ·z(x² + y²) + ε(3x²y − y³)           δ=1.0567, ν=2.3978,
                                                         λ=1.2975, ε=−1.5491
```

In complex form this is the n = 3 member of a new dihedral family
([`atlas/family_d.py`](atlas/family_d.py)):

```
dw/dt = (α·z² − β + iω·z)·w + γ·w̄² + iδ·z·w̄²
dz/dt = ν·z − z³ − λ·z|w|² + ε·Im(w³)
```

Three coefficient sets would enlarge the group, and none of them is zero here:
`ω = δ = ε = 0` restores the z-flip and the vertical mirrors (D₃ₕ),
`ω = γ = 0` restores the diagonal mirrors and the inversion (D₃d), and
`γ = δ = ε = 0` restores rotation by every angle. Equivariance under all six
elements, and violation of every larger group (C₂ₙ, z-flip, mirrors, inversion,
S₂ₙ, and every copy of O and I that contains Dₙ), is unit-tested for the
family at n = 2 to 6 ([`tests/test_d_equivariance.py`](tests/test_d_equivariance.py);
the O and I copies at n = 2 to 5)
and for the specimen package
([`tests/test_trinacria_matches_family.py`](tests/test_trinacria_matches_family.py)).

The specimen came from a fine-dt screen of 13 nominees of a random probe of the
family box ([`scripts/search_d3.py`](scripts/search_d3.py),
[`results/d3_selection.json`](results/d3_selection.json)). No nominee passed
every pre-registered gate. Trinacria (nominee #8) passed all but the basin
gate, and it was accepted under a relaxed basin gate, disclosed below. Its
parameters are the 4-decimal probe values with no polishing; β = 1.4668 sits
near the box face 1.5.

Certified via `python scripts/verify_trinacria.py`
([`results/trinacria_verification.json`](results/trinacria_verification.json)).
The battery is stricter than Mobula's: it adds a dt-halving arm, a
spectrum-versus-divergence check, a D₃-orbit test of the equilibria, stronger
symmetry controls, and a single-attractor test.

| Property | Value |
|---|---|
| Lyapunov spectrum | **λ₁ = +0.328 ± 0.002**, λ₂ = −0.0006, λ₃ = −1.492 (mean ± SE over 24 orbits, T = 2000 each, dt = 0.00125); a single orbit scatters by ±0.009, so quote **λ₁ = 0.33 ± 0.01** |
| Kaplan–Yorke dimension | **D_KY ≈ 2.219** |
| dt halving (same ICs, same T) | λ₁ = 0.332, 0.332, 0.332, 0.328 at dt = 0.01, 0.005, 0.0025, 0.00125; successive changes 0.2%, 0.2%, 1.3%; no monotone drift; spot check at dt = 0.000625: 0.335 |
| Σλ versus mean divergence | −1.1649 versus −1.1639 (0.09% apart) |
| λ₁ across 5 near-attractor ICs | 0.305 … 0.358 (spread 16%) |
| Boundedness | 4M steps at dt = 0.004, max \|X\| = 1.573; mean divergence −1.164 |
| Centring (64 orbits, 2.75 × 10⁵ time units) | centroid / extent = (0.0023, 0.0004, 0.0017); the x offset is 2.4 SE, inside the 3 SE gate |
| Equilibria | 24, in six complete D₃-orbits of sizes 1, 2, 3, 6, 6, 6; spectra agree within each orbit to 2 × 10⁻¹³ |
| Origin | a saddle with the double real eigenvalue −β = −1.467 and ν = 2.398 along z (not a saddle-focus, unlike the C_n specimens) |
| On-axis pair (0, 0, ±1.548) | κ-conjugate saddle-foci, identical eigenvalues 2.075 ± 1.110i, −4.796, and **opposite rotation sense** (J_xy − J_yx = ±2.220) |
| In-plane triple | saddle-foci at radius 1.737 exactly on the three half-turn axes, eigenvalues −2.958 ± 7.749i, 1.467 |
| D₃ symmetry residual | all five non-identity elements 0.0044–0.0046, against a self floor of 0.0045 (ratio 1.006) |
| Controls (not symmetries) | z-flip 0.070 (15×), y-mirror 0.070 (15×), C₆ 0.113 (25×), inversion 0.076 (17×), 37° rotation 0.083 (18×), 144 elements of O and I containing D₃ ≥ 18× |
| Single attractor, not a κ-conjugate pair | 48/48 near-attractor ICs land on the same symmetric set; orbits launched from κx₀ and ρx₀ fill the same cloud (cross residuals 0.0047–0.0054, at the floor); P(z > 0) = 0.504 ± 0.002, pooled KS(z, −z) = 0.012 |
| Counter-rotation (64-orbit ensemble) | mean angular velocity +0.202 ± 0.001 above the mid-plane, −0.207 ± 0.001 below |
| **Basin (disclosed)** | **24/48** random ICs in [−2, 2]³ stay bounded (56 ± 2% of a 480-IC sample) and **every bounded orbit lands on the symmetric attractor**; the rest escape to infinity. Same protocol: Aurelia 48/48, Naiad 48/48, Cassiopea 43/48, Mobula 41/48 |
| Novelty vs the dysts catalog | nearest Coullet at 0.212, but JerkCircuit (0.241) is within the 0.055 noise floor, so the nearest-neighbour identity is not stable; the shape is distinct (above the 0.115 median spacing) |
| Distance to the other specimens | Cassiopea 0.217, Aurelia 0.290, Naiad 0.626, Mobula 0.644 |

What the numbers do not say:

- **Escape.** Trinacria has the smallest basin of the five. Half of a
  [−2, 2]³ box escapes to infinity. The pre-registered gate (48/48 land) fails
  and is recorded as failing. The specimen is certified under the relaxed gate
  "every bounded orbit lands on the symmetric set" (270/270 in the large
  sample), a choice made explicitly and recorded with its reasons in
  [`results/d3_selection.json`](results/d3_selection.json) (`user_decision`).
- **Slow mixing between the halves.** The orbit stays on one side of the
  mid-plane for up to about 40 time units, so single-window averages
  (centroid, P(z > 0), KS) scatter by more than their gates allow. Those gates
  are taken on the 64-orbit ensemble with across-orbit errors. Two ensemble
  statistics whose exact value is zero sit at 2.4 SE (x centroid) and 2.8 SE
  (above/below angular-velocity mismatch, 2.4% of the mean speed); they are
  reported as they are.
- **Gates changed after seeing results**, each kept in the JSON with its
  original verdict: the basin gate (above); the z-rotation scan control, whose
  first version used the 5° rotation (2.9×, trivially close to the identity)
  and now uses the angles a larger dihedral group would add; and the
  image-launch fingerprint distance, which is reported but not gated because
  the fingerprint ignores rotations.

**Literature.** Dihedral symmetry in chaotic *maps* is classical (Field and
Golubitsky, Chossat and Golubitsky, Reiter). For flows, Letellier and Gilmore
(*J. Phys. A* 40, 5597, 2007) classify the symmetry groups available to 3-D
systems and build equivariant flows by lifting an image system through a
covering map; that construction can produce D₃-equivariant chaotic flows, and
we could not read the full text from this environment. Nearby symmetric flows
that we checked are not D₃: Thomas' cyclically symmetric system is S₆, three
identical coupled cells act as C₃ᵥ or D₃d, and the Qi et al. (2008) four-wing
system is C₂. So the claim is narrow: Trinacria is the first D₃ specimen of
this atlas and, to our knowledge, the first designed and certified polynomial
flow whose attractor carries exactly the chiral group D₃. It is not a claim that
the group is new to chaotic flows. Details and the search log are in
[`docs/RELATED_WORK.md`](docs/RELATED_WORK.md) (Part B, Claim 4).

## The polyhedral wall: an irreducibility obstruction

The [symmetry program](docs/SYMMETRY_PROGRAM.md) set out to climb past the cyclic
and rotoreflection groups to the **Platonic** ones — a strange attractor with the
rotation symmetry of a tetrahedron (T, order 12), an octahedron (O, 24), and, as
the intended signature image, an icosahedron (I, 60). The equivariant families
were constructed and their symmetry proven to machine precision
(`atlas/family_{t,o,i}.py`, `tests/test_{t,o,i}_equivariance.py` — 22 tests
including the degree-6 icosahedral invariant built as `Σ_{g∈I} ((g·v)·X)⁶`).

**No chaotic flow exists in them — and that is the finding.** The cyclic and
rotoreflection groups act *reducibly* on R³ (a plane rotation plus a height
axis), which is exactly what supplies the saddle-focus engine behind all four
jellyfish. The polyhedral rotation groups act *irreducibly*, so by Schur's lemma
the only commuting linear map is `ρ·I`: **there is no linear engine.** An
exhaustive, dt-robust search (polynomial degrees 3–5, gradient and non-gradient
rotational equivariants; uniform-random and fixed-point-stability-targeted;
~10⁵ evaluations) finds no strange attractor — bounded orbits relax to equilibria
or limit cycles. Apparent chaos at a coarse step is an RK4 artifact: a candidate
reading λ₁ ≈ +0.4 at dt = 0.01 collapses to λ₁ ≈ −20 (a stable fixed point) at
dt = 0.0025. Even trigonometric "labyrinth" flows (the Thomas mechanism) give
only weak (λ₁ ≈ 0.04), spatially *localized* — i.e. symmetry-broken — chaos,
because a genuinely symmetric attractor must wrap the origin, where the flow
reduces to its non-chaotic polynomial core. This is strong numerical evidence
(not a proof) for why the literature has chaotic *maps* with polyhedral symmetry
but no chaotic *flow*. What the polyhedral groups host is the rung below chaos:
equilibria and symmetric limit cycles.

## The chaos-vs-symmetry-group law

Giving every group the same search budget and the same fine-dt re-certification
([`scripts/symmetry_law.py`](scripts/symmetry_law.py) → [`results/symmetry_law.csv`](results/symmetry_law.csv))
turns the bestiary into one quantitative claim. The robust comparable is the
**maximum largest-Lyapunov-exponent that survives long-integrator certification**
under equal evaluation budget (the chaotic *fraction* is only loosely comparable
across families with different parameter counts, so it is reported but not
headlined):

| group | action on R³ | \|G\| | max certified λ₁ | named specimen |
|---|---|---:|---:|---|
| C₂ | reducible | 2 | 0.41 | Naiad (λ₁ 0.296) |
| C₃ | reducible | 3 | 0.59 | Aurelia (0.233) |
| C₄ | reducible | 4 | 0.56 | Cassiopea (0.525) |
| C₅ | reducible | 5 | 0.51 | — |
| C₆ | reducible | 6 | 0.27 | — |
| C₇ | reducible | 7 | 0.33 | — |
| C₈ | reducible | 8 | 0.18 | — |
| S₄ | reducible (improper) | 4 | 0.63 | Mobula (0.330) |
| D₃ | reducible (chiral dihedral) | 6 | 0.69 | Trinacria (0.33) |
| **T** | **irreducible** | 12 | **0.00** | — (obstruction) |
| **O** | **irreducible** | 24 | **≈0**† | — (obstruction) |
| **I** | **irreducible** | 60 | **≈0**† | — (obstruction) |

†The lone near-threshold polyhedral survivors are not genuine symmetric attractors:
their λ₁ fails to converge as dt → 0, and their orbits are off-center and
anisotropic (symmetry-broken), not invariant sets of the full group.

The max certified λ₁ is the *family maximum* found by the matched search; it
exceeds the curated named specimen's λ₁ (specimens were chosen for beauty and
novelty, not maximal stretching). Two regimes are visible
([`gallery/symmetry_law.png`](gallery/symmetry_law.png)): the reducible-action
groups host chaos with a **non-monotonic** dependence on order, while the
irreducible-action polyhedral groups fall off an **irreducibility cliff** to zero.
The amount of chaos a symmetry permits depends less on *how much* symmetry there
is than on *what kind*. The D₃ row, added with Trinacria, makes the same point
within one order: C₆ and D₃ both have six elements, and their maxima are 0.27
and 0.69. (The D₃ maximizer was checked to be centred and symmetric; it has one
more parameter than the Cₙ family, so the caveats above apply.) (Null check, after Letellier–Gilmore: in a cover/image
construction the cover's exponent is independent of the covering group; our
families are independently searched, not covers of a shared image, so the
group-dependence is a real property of the equivariant design space.)

## Reproduce everything

```bash
pip install -r requirements.txt

python scripts/verify_chaos.py        # certify Aurelia, write results/verification.json
python scripts/render_gallery.py      # render the Aurelia gallery
python scripts/bifurcation_scan.py    # sweep c, plot the route to chaos
python scripts/search_parameters.py   # re-run the original discovery search
python scripts/verify_naiad.py        # certify Naiad, write results/naiad_verification.json
python scripts/render_naiad.py        # render the Naiad gallery
python scripts/verify_cassiopea.py    # certify Cassiopea
python scripts/render_cassiopea.py    # render the Cassiopea gallery
python scripts/verify_mobula.py       # certify Mobula (S4 rotoreflection)
python scripts/render_mobula.py       # render the Mobula gallery
python scripts/search_d3.py           # screen the 13 D3 nominees (~40 min)
python scripts/verify_trinacria.py    # certify Trinacria (D3 dihedral, ~10 min)
python scripts/render_trinacria.py    # render the Trinacria gallery
python scripts/basin_comparison.py    # shared 48-IC basin test, all five specimens
python scripts/convergence_study.py   # error-barred Lyapunov spectra
python scripts/validate_fingerprint.py # validate the novelty metric, draw the shape map
python scripts/homoclinic_hunt.py     # hunt Shilnikov homoclinic connections
python scripts/symmetry_law.py all    # chaos-vs-symmetry law: C2..C8, S4, D3, T, O, I + figure
python -m pytest tests/               # numerical equivariance + cross-check tests
```

Or use it as a library:

```python
from aurelia import trajectory

pts = trajectory(n_steps=1_000_000)   # (N, 3) points on the attractor
```

## License

MIT — see [LICENSE](LICENSE).
