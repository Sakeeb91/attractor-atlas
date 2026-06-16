# Implications for Network Theory and Computer Science

*What the Attractor Atlas — three certified C₂/C₃/C₄-equivariant strange attractors
plus the quality-diversity engine that finds them — means outside dynamical systems.*

This repository is, on its face, a chaos project: a `Cₙ`-equivariant family of 3-D
polynomial flows (`atlas/family.py`), a MAP-Elites discovery engine
(`atlas/map_elites.py`), a shape-fingerprint novelty metric (`atlas/fingerprint.py`),
and a certification battery (`scripts/verify_*.py`). But the *machinery* it uses and
the *objects* it produces sit on well-trodden bridges into computer science and
network theory. This note separates **what the code already is** (a CS artifact, read
correctly) from **what it implies and enables** (concrete bridges, with the prior art
that grounds each one).

---

## Part A — Computer Science

### A1. The discovery loop is quality-diversity / illumination, not optimization

The engine in `atlas/` is a textbook **MAP-Elites** (Mouret & Clune 2015) driven by
**novelty search** (Lehman & Stanley 2011) — two of the foundational algorithms of
*open-ended* evolutionary computation. The CS-relevant move is the inversion of the
objective:

- `descriptor_cell(n, lle, aspect)` tiles a *behavior* space (symmetry order × largest
  Lyapunov exponent × vertical aspect) into cells; each cell keeps one elite
  (`Archive.consider`). This is **illumination** — "fill a map of diverse solutions" —
  rather than "find the single best," and it is exactly the antidote to the mode
  collapse that single-objective evolution suffers (Pan & Das 2015 evolved hundreds of
  chaotic flows by maximizing λ₁ and got Lorenz twins; QD's cells + novelty exist to
  avoid that).
- The selection pressure is **novelty as distance to a reference set**
  (`fingerprint.novelty`), not a fixed fitness. This is the practical realization of
  "abandoning objectives."

The transferable CS contribution here is the **domain of the novelty signal**: instead
of novelty relative to the current population (the usual recipe), this measures novelty
against an *external scientific catalog* (the 135-system `dysts` set,
`results/catalog_fingerprints.json`) plus the live archive. That is a reusable pattern
for any "discover things genuinely new to the field" problem — drug scaffolds, circuit
topologies, proof tactics, program sketches — where the reference set is the corpus of
what humanity already knows, and the search is explicitly pushed to its frontier.

### A2. Novelty = nearest-neighbor search in a learned-free embedding

`atlas/fingerprint.py` is a small but complete **content-based retrieval** system, and
recognizing it as such connects this project to information retrieval and vector
databases:

- The **D2 shape distribution** (Osada et al. 2002, ACM TOG) — a histogram of pairwise
  distances over a normalized point cloud — is a classic 3-D shape descriptor from
  computer graphics, chosen here precisely because it is invariant to translation,
  scale, and rigid rotation (the things that don't matter), exactly the invariances a
  good retrieval key needs.
- The Hellinger embedding (`hist = sqrt(hist/sum)`) turns the histogram into a vector
  whose Euclidean distance approximates a statistical divergence — the same trick used
  in bag-of-features image search.
- `novelty(fp, reference)` is a brute-force **k-NN (k=1) query** in a 35-dimensional
  feature space. At catalog scale this is fine; the moment the archive grows, this is
  the textbook entry point for ANN indexing (FAISS, HNSW), and the whole novelty loop
  becomes a streaming nearest-neighbor problem — a recognizable systems/algorithms
  engineering task rather than a chaos one.

The validation story (`scripts/validate_fingerprint.py`: Lorenz variants cluster =
positive control; algebraic Sprott family does *not* cluster = negative control; a
~0.05 within-system noise floor) is, in CS terms, a **retrieval-quality evaluation**:
precision on known-positive pairs, a null model, and a measurement-noise floor. It is
how you'd validate any similarity metric before trusting its rankings.

### A3. Equivariance "by construction" is the same idea as geometric deep learning

`atlas/family.py` builds its vector field from exactly the terms that survive a group
action: for `Cₙ` the only symmetry-breaking coupling allowed is the resonant
`conj(w)^(n-1)` (the textbook `Zₙ`-resonant normal-form term, Golubitsky–Stewart–
Schaeffer Vol. II). `docs/SYMMETRY_PROGRAM.md` makes the recipe explicit: enumerate the
monomials that commute with the symmetry, keep those, forbid the rest, and *unit-test
equivariance numerically* before searching.

That "keep only the group-equivariant maps" construction is **identical in spirit to
equivariant / group-convolutional neural networks** (Cohen & Welling 2016) and the
broader geometric-deep-learning program (Bronstein et al. 2021): both restrict a
function family to those that respect a symmetry, so the symmetry is guaranteed rather
than learned. The Atlas is a clean, low-dimensional, fully interpretable instance of
that design principle — a useful teaching/benchmark artifact for it, and a hint that
the discovery engine itself could be made symmetry-aware end to end.

### A4. Automated scientific discovery + reproducibility-as-code

The repo is a concrete "AI-for-Science" pipeline: propose → integrate → measure (λ
spectrum, Kaplan–Yorke dimension, equilibria) → screen for novelty → keep, certify, and
name. Its CS-flavored discipline is worth stating explicitly because it's the part most
reusable:

- **Every claim has a script and a `results/*.json`** (the `verify_*` / `convergence_*`
  / `homoclinic_hunt` scripts each emit a machine-readable artifact). This is
  reproducible-research engineering — claims are executable, not prose.
- The **archive persists across runs** (`Archive.save/load`), so search is incremental
  and compounding — a checkpointed, resumable long-running job, not a one-shot script.
- The novelty objective is *measured against a public benchmark*, which is the same
  hygiene that makes ML leaderboards meaningful.

It also sits in a recognizable lineage of automated discovery — Sprott's 1993 aesthetic
map search, Goertzel's 1995 GA, Pan & Das 2015 genetic programming, and the 2025 "Panda"
corpus of ~20,000 evolved ODEs used as *training data* for a forecasting foundation
model. The distinction this project draws (QD + catalog-grounded novelty + per-system
certification, vs. corpus generation) is precisely an argument about *what kind of
search produces scientifically usable artifacts* — a methods question that generalizes
well beyond attractors.

### A5. The evaluator is SIMD/HPC engineering

`evaluate_batch` integrates an **entire population simultaneously** with vectorized
complex RK4 and a vectorized two-orbit Benettin Lyapunov estimate, culling diverged
orbits in place (`cull()`) so dead candidates stop generating overflow. This is
data-parallel scientific computing: a batch dimension that maps directly onto SIMD/GPU
lanes, branch-free masking instead of per-candidate control flow, and a fixed-work inner
loop. The comment in `fingerprint.py` about NumPy 2.1's histogram fast-path miscounting
edge-landing float64 values is a real **numerical-correctness** war story — the kind of
floating-point subtlety that decides whether a measured invariant is trustworthy.

### A6. Where chaos meets classical CS

The *objects* produced — certified chaotic flows with tunable, group-indexed structure —
plug into several long-standing CS areas:

- **Pseudo-randomness and chaos-based cryptography.** A strongly positive largest
  Lyapunov exponent (λ₁ ≈ 0.23–0.53 across the family) is sensitive dependence —
  exponential decorrelation — which is the property chaotic PRNGs and chaos cryptosystems
  exploit. A *parametric, symmetry-indexed* family of certified chaotic sources with
  known spectra is a tidy substrate for studying keystream quality vs. a structural
  knob (the rotation order `n`), with the usual sharp caveat that empirical chaos is
  necessary but not sufficient for cryptographic security.
- **Predictability limits and computability.** The Lyapunov time `1/λ₁` is a hard horizon
  on forward prediction and a clean illustration of how finite-precision computation
  loses information about a deterministic system — the practical face of questions about
  the computability and complexity of real dynamical systems (e.g. the difficulty of
  *deciding* chaos, and shadowing as the reason simulations remain meaningful at all).
- **New benchmarks for time-series ML / Koopman / reservoir computing.** `dysts` itself
  exists as an *interpretable forecasting benchmark* (Gilpin 2021, NeurIPS D&B). The
  three new systems — with full spectra, error bars, and a *measured* shape distance
  from everything in `dysts` — are ready-made out-of-distribution test cases: their
  isolation in shape space is exactly what you want when probing whether a learned
  forecaster or Koopman/operator-learning model generalizes off-catalog, and the
  `Cₙ` symmetry makes them targeted tests for *symmetry-aware* forecasters.

---

## Part B — Network Theory

The Atlas is not a graph project, so the connections here are bridges — but they are
load-bearing ones, because the repository's organizing idea (**symmetry groups**) is the
same idea that governs dynamics on symmetric *networks*.

### B1. The symmetry-group framing *is* the network cluster-synchronization story

`docs/SYMMETRY_PROGRAM.md` proposes a "periodic table of strange attractors indexed by
symmetry group," and the family is built so that group equivariance forces invariant
structure in phase space. This is mathematically the **same machinery** that network
science uses for **cluster synchronization**: the automorphism group of a network
partitions its nodes into orbits, and those orbits predict which nodes synchronize and
how the synchronous state loses stability (Pecora, Sorrentino, Hagerstrom, Murphy &
Roy, *Nat. Commun.* 2014; Golubitsky & Stewart's groupoid formalism for networks).

The correspondence is concrete:

| Atlas (this repo) | Symmetric networks |
|---|---|
| `Cₙ` equivariance of the vector field | automorphism group of the graph |
| invariant coordinate axis / symmetric attractor support | synchronization manifold / invariant subspace |
| `conj(w)^(n-1)` = lowest term respecting the symmetry | equivariant coupling that preserves cluster partitions |
| Melbourne–Ashwin "average symmetry of an attractor" | symmetry of the synchronous pattern actually observed |

So the project is, unintentionally, a low-dimensional **sandbox for the group theory of
synchronization**: every result here about how a discrete symmetry shapes an attractor
has a direct analogue in how a network's automorphisms shape its synchronous patterns.

### B2. The certified units are node dynamics for coupled-oscillator networks

Each system is a self-sustained oscillator (the `i·ω·w` term sets an internal rotation
rate) that is also chaotic. That is exactly the ingredient network dynamics wants for
studying **chaos synchronization** and **network stability**:

- Couple `N` copies of Aurelia/Naiad/Cassiopea diffusively on a graph and the
  synchronizability is governed by the **Master Stability Function** (Pecora & Carroll
  1998) evaluated along *these* systems' transverse Lyapunov spectra — which this repo
  already computes per system.
- Because the units carry a *known* discrete symmetry, the network's own symmetry and
  the node symmetry interact, giving a controlled setting for symmetry-induced cluster
  states. Most MSF studies use Lorenz/Rössler nodes; a symmetry-indexed family is a new,
  structured set of node models with a tunable knob (`n`, and the family's
  non-monotonic chaos strength: λ₁ for C₂/C₃/C₄ = 0.296/0.233/0.525).

### B3. The equilibrium structure is a heteroclinic network in waiting

Network theory studies not only graphs of nodes but **heteroclinic networks** in phase
space — sets of saddle equilibria/cycles joined by heteroclinic connections, used as
models of *winnerless competition* and sequential neural computation (Ashwin, Field;
Rabinovich et al.). This repo's equilibrium data is a direct entry point:

- Aurelia/Naiad have a *single* organizing saddle-focus (Shilnikov anatomy) — the
  `homoclinic_hunt` script (`results/homoclinic_hunt.json`) already searches for the
  homoclinic connection that would close a one-node loop.
- **Cassiopea is the interesting one**: a central saddle-focus *plus a `C₄` quadruple of
  off-axis saddles* (four saddles at radius 3, in a square — see README). That is the
  skeleton of a *symmetric heteroclinic network*: five vertices, with the symmetry group
  acting on the four outer ones. Mapping the unstable-manifold connections among those
  five equilibria would turn Cassiopea's phase space into an actual directed graph whose
  structure is forced by `C₄` — a concrete, small, symmetric heteroclinic-network
  specimen produced as a by-product of the search.

### B4. Each attractor is a generator of complex networks (time-series → graph)

A standard bridge from dynamics to network science is to **turn a trajectory into a
graph** and study its topology: **recurrence networks** (Marwan, Donner et al.),
**visibility graphs** (Lacasa et al.), and **ordinal/transition networks**. Each
construction maps a chaotic time series to a complex network whose degree distribution,
clustering, and community structure fingerprint the dynamics.

The three certified orbits (`from aurelia import trajectory`, and the siblings) are
ready-made generators for this. Two things make them more than yet another input:

1. **Symmetry should imprint on the network.** A `Cₙ`-symmetric attractor should yield a
   recurrence/transition network with an approximate `Cₙ` automorphism — a clean,
   designed test case for "does network symmetry detect dynamical symmetry?"
2. The fractal (Kaplan–Yorke) dimension and Lyapunov spectrum are *certified*, so any
   correlation found between the derived network's topology (e.g. recurrence-network
   transitivity ↔ dimension, established in the literature) and the known invariants is
   checkable against ground truth here.

### B5. The archive and fingerprint map are themselves a network of dynamical systems

The clearest network object is hiding in plain sight. `fingerprint.distance` defines a
metric over *systems*; the MDS picture in `gallery/fingerprint_map.png` (Kruskal
stress-1 = 0.13) is an embedding of the **similarity graph** whose nodes are dynamical
systems and whose weighted edges are fingerprint distances. Standard network-science
questions apply directly and are partly already answered:

- **Outlier / centrality:** Aurelia at nearest-neighbor distance 0.185, Naiad at ~0.3,
  against a catalog median NN spacing of 0.115 — these are *peripheral, low-degree nodes*
  in a k-NN graph. "83% of known systems sit closer to something than Aurelia does" is a
  network statement about the density of the neighborhood.
- **Community structure:** the validation that Lorenz-family variants cluster while the
  algebraic Sprott family does not is exactly a **community-detection** result on the
  similarity graph — geometry-defined communities exist, algebra-defined ones don't.
- **Graph construction choices matter:** building a k-NN graph vs. an ε-graph over these
  fingerprints, and asking about its connectivity, hubs, and the placement of the three
  new systems on the periphery, is an off-the-shelf network analysis that would sharpen
  the "how isolated is novel?" question the repo already cares about.

---

## What this concretely enables (shortlist)

1. **A symmetry-indexed node library for network dynamics**: package the family as
   drop-in oscillators for MSF / cluster-synchronization studies, with the symmetry knob
   as the independent variable.
2. **Cassiopea's 5-equilibrium graph**: compute the heteroclinic connections among the
   central saddle-focus and the `C₄` outer saddles → a designed symmetric
   heteroclinic-network specimen.
3. **Recurrence/visibility-network analysis** of the three certified orbits, testing
   whether dynamical `Cₙ` symmetry shows up as network automorphism.
4. **Treat the fingerprint archive as a graph**: k-NN/ε-graph over `dysts` + discoveries,
   then community detection and centrality — turning "novelty distance" into network
   isolation, with the three discoveries as labeled peripheral nodes.
5. **Off-catalog ML benchmarks**: ship the three systems (with certified spectra and
   measured catalog distance) as out-of-distribution test cases for forecasting / Koopman
   / reservoir-computing models, including symmetry-aware variants.

## Honest caveats

- The repository contains **no graphs today**; every Part B item is an enabled direction,
  not an implemented result. The strong claim is only that the project's *symmetry-group
  organizing principle* coincides with the mathematics network science already uses for
  synchronization — so the bridge is short, not that it has been crossed here.
- The CS connections in Part A2–A3 and A5 are *descriptions of what the code already is*
  (retrieval, equivariant construction, batched numerics); A1 and A4 are established
  (MAP-Elites / novelty search / reproducible-research), with the novel twist being the
  catalog-grounded novelty signal. A6's cryptography link is real but carries the
  standard warning: chaos is necessary, not sufficient, for security.
- "Implication" ≠ "validated finding." None of these bridges has been measured in this
  repo; this document is a map of where the work *reaches*, with the prior art that would
  anchor each reach.

### Pointers

Code: `atlas/family.py` (Cₙ family + batched evaluator), `atlas/map_elites.py` (QD
archive), `atlas/fingerprint.py` (D2 + PCA novelty), `scripts/run_atlas.py` (search
loop), `scripts/homoclinic_hunt.py` (saddle-connection search). Context:
`docs/SYMMETRY_PROGRAM.md` (the symmetry-group program), `docs/RELATED_WORK.md`
(positioning and citations).
