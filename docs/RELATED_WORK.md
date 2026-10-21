# Related work

This document unions the two literature passes that ran concurrently on this
repository and were reconciled at the `main` → `symmetry-program` merge:

- **Part A — Methods-paper positioning.** The literature-positioning survey for
  the methods paper: what prior work did, how this project differs, and which
  public claims it forces us to reword. Authoritative verdict table at the end of
  the part.
- **Part B — Phase-0 polyhedral, rotoreflection & chaos-law verification.** The
  symmetry program's Phase-0 "literature verification strike": targeted checks of
  the three novelty claims that underpin the polyhedral (T/O/I) specimens, the S₄
  rotoreflection specimen, and the chaos-vs-symmetry-group law.
  Claim 4 (added 2026-10-21) checks dihedral D_n flows before any wording about the
  D3 specimen Trinacria.

The two parts overlap deliberately (both cite Field–Golubitsky, Letellier–Gilmore,
GSS, Ashwin–Melbourne) but from different angles: Part A positions the *method*
(quality-diversity search over equivariant families), Part B positions the
*specimens and the law*. Where a number appears in both, Part A's validated ranges
(noise-floor calibrated in the fingerprint-validation study) take precedence.

---

# Part A — Methods-paper positioning

Literature positioning pass for the methods paper (gap item 1 of
[PAPER_PLAN.md](PAPER_PLAN.md)). Every entry records what the prior work did, how this
project differs, and which of our public claims it forces us to reword. The verdict
table at the end is the authoritative list of claims that survived. Searches were run
June 11, 2026 (web search, arXiv, Crossref); bibliographic data below was verified
against Crossref/arXiv records, not quoted from memory.

*A companion document with the same filename exists on the `symmetry-program` branch,
scoped to that program's three novelty claims (polyhedral flows, rotoreflection
groups, the chaos-vs-group law). When the branches merge, the two documents should be
concatenated as separate top-level parts.*

**Headline verdicts.**

1. The broad claim "no published 3-D flow with discrete Cn equivariance" **does not
   survive**: equivariant covers (Miranda–Stone, Letellier–Gilmore), Thomas's
   cyclically symmetric system, and the rotation-symmetric multi-scroll/multi-wing
   literature are all counterexamples. The narrow claim that survives is stated in
   §2 below.
2. The claim that the *method* — quality-diversity search over an equivariant family
   with novelty measured as fingerprint distance to a catalog of known attractors —
   has no published precedent **survives**, but must be worded against two nearby
   ancestors: genetic programming over Lorenz-like flows (Pan & Das 2015) and the
   evolutionary 20,000-system Panda corpus (Lai, Bao & Gilpin 2025).
3. "Exactly one equilibrium, which is unusual" **does not survive as stated**:
   one-equilibrium chaotic flows are an established genre (jerk systems; the
   one-stable-equilibrium hidden-attractor literature). What is distinctive is the
   *organizing role* of the single saddle-focus together with the global symmetry.
4. The names (Aurelia, Naiad, Cassiopea attractors) re-checked June 11, 2026:
   still unclaimed. **Survives.**

---

## 1. Symmetric chaos: the theory we inherit

**Field, M. & Golubitsky, M., *Symmetry in Chaos: A Search for Pattern in
Mathematics, Art, and Nature*. Oxford University Press, 1992; 2nd ed., SIAM, 2009.**
What they did: built chaotic *maps* of the plane equivariant under cyclic and dihedral
groups (the "symmetric icons") and under crystallographic groups ("quilts"),
explicitly motivated by visual beauty; developed the machinery of equivariant maps,
admissible symmetries, and symmetry detectives. How we differ: their icons are
two-dimensional discrete maps; our systems are three-dimensional continuous
autonomous flows. The aesthetic-plus-equivariance motivation is the same, and the
paper should present the family as a flow-analogue of the icon program. Claims to
reword: none — README/STORY already state the maps-vs-flows distinction correctly.

**Golubitsky, M., Stewart, I. & Schaeffer, D., *Singularities and Groups in
Bifurcation Theory, Vol. II*. Applied Mathematical Sciences 69, Springer, 1988.**
What they did: equivariant bifurcation/normal-form theory. For a planar system with
Z_q rotation symmetry, the lowest-order non-axisymmetric equivariant ("resonant")
term is precisely the conjugate power w̄^(q−1). How we differ: we do not extend the
theory; we *use* the resonant term as a global polynomial ansatz, coupled to an
axial cubic. Claims to reword: any phrasing that presents w̄^(n−1) as a novel
invention ("a small piece of algebraic luck") should carry the citation: the term is
the textbook Z_n-resonant coupling; the design choice is using it globally, as the
*only* symmetry-breaking term, in a 3-D flow tuned for chaos.

**Melbourne, I., Dellnitz, M. & Golubitsky, M., "The structure of symmetric
attractors," *Arch. Rational Mech. Anal.* 123(1):75–98, 1993.
Ashwin, P. & Melbourne, I., "Symmetry groups of attractors," *Arch. Rational Mech.
Anal.* 126(1):59–78, 1994.** What they did: classified which subgroups of a symmetry
group can occur as the (average) symmetry of an attractor; introduced the
instantaneous-vs-average symmetry distinction. Relevance: gives the correct
vocabulary for our observation that each attractor carries the full C_n symmetry on
average; also frames the open question about symmetry-broken siblings and
symmetry-increasing crises. Claims to reword: none; cite in the paper's method
section.

## 2. C_n-equivariant 3-D flows that already exist (the claim-killers)

**Miranda, R. & Stone, E., "The proto-Lorenz system," *Phys. Lett. A*
178(1–2):105–113, 1993. DOI 10.1016/0375-9601(93)90735-I.** What they did:
constructed the quotient ("proto-Lorenz") of the Lorenz system by its Z_2 symmetry,
then built *n-fold covers* of the proto-Lorenz, yielding flows whose attractors have
n-fold rotation symmetry about an axis ("n-eared" Lorenz-like attractors), for
arbitrary n — published thirty-three years before Aurelia. How we differ: the
covers are obtained by *lifting known dynamics* through the angle-multiplying map;
away from the rotation axis a cover is locally diffeomorphic to its image, so the
dynamics (Lyapunov exponents, branched-manifold structure modulo symmetry) is the
Lorenz dynamics re-presented, and the lifted vector fields are not polynomial in
the natural coordinates (they carry inverse powers of the cylindrical radius and
are typically singular or degenerate on the rotation axis). Aurelia/Naiad/Cassiopea
are single globally smooth polynomial vector fields, not covers of any catalogued
system, and their attractors enclose the symmetry axis where their organizing
saddle-focus sits. Claims to reword: every occurrence of "no published 3-D flow
with discrete C_n equivariance" (README §Novelty's relatives framing, STORY §The
name).

**Letellier, C. & Gilmore, R., "Covering dynamical systems: Twofold covers,"
*Phys. Rev. E* 63:016206, 2001. DOI 10.1103/PhysRevE.63.016206. And Gilmore, R. &
Letellier, C., *The Symmetry of Chaos*. Oxford University Press, 2007.** What they
did: developed the systematic theory of cover/image ("equivariant/invariant")
dynamical systems: any 3-D flow can be quotiented by, or lifted to, discrete
rotation symmetries, with explicit algorithms, including three- and four-fold covers
of Rössler- and Lorenz-class systems and analysis of how the rotation-axis position
changes the topology. This is the definitive prior art on discrete rotational
symmetry in 3-D chaotic flows. How we differ: same as for Miranda–Stone — cover
systems re-present the dynamics of a known image system (that is the point of the
construction); they are obtained top-down from a chosen image and symmetry, not
discovered bottom-up by search, and are generically non-polynomial/singular at the
axis. Our claim must therefore be about the *family and how it was found*, not about
equivariance per se. Claims to reword: same as above; additionally, the paper's
related-work section must cite this book prominently (a referee from this community
is the most likely reviewer).

**Thomas, R., "Deterministic chaos seen in terms of feedback circuits: Analysis,
synthesis, 'labyrinth chaos'," *Int. J. Bifurcation Chaos* 9(10):1889–1905, 1999.
DOI 10.1142/S0218127499001383.** What he did: the cyclically symmetric system
ẋ = sin y − bx (and cyclic permutations) — invariant under the coordinate
permutation (x,y,z)→(y,z,x), which is geometrically a 2π/3 rotation about the space
diagonal (1,1,1). A published 3-D flow with a C_3 rotational symmetry, in every
standard catalog. How we differ: Thomas's symmetry permutes the three coordinates
(rotation about the diagonal, no invariant coordinate plane); ours acts in a
distinguished plane with an invariant axis, which is what produces the
jellyfish-bell anatomy (planar petals + axial column). Thomas's system also has
many equilibria (a lattice of them in the labyrinth limit), not a single
saddle-focus. Claims to reword: "no published 3-D flow with discrete C_n
equivariance" again — Thomas alone falsifies it.

**Lü, J. & Chen, G., "Generating multiscroll chaotic attractors: Theories, methods
and applications," *Int. J. Bifurcation Chaos* 16(4):775–858, 2006. DOI
10.1142/S0218127406015179.** What they did: the survey of the multi-scroll
literature: families of 3-D systems engineered (piecewise-linear switching,
staircase/hysteresis nonlinearities, sine functions) to produce any prescribed
number and arrangement of scrolls, including rotationally arranged ("symmetrical
multi-scroll") attractors. How we differ: multi-scroll systems are *designed* — the
scroll count and placement is the input, typically via non-smooth nonlinearities —
whereas our systems are smooth low-degree polynomial flows found by search inside a
fixed equivariant ansatz, with the attractor's structure an output. Claims to
reword: the README "multi-scroll/multi-wing families" passing mention stays; the
positioning sentence should make the smooth-polynomial-vs-engineered distinction
explicit.

**Yang, Y., Huang, L., Xiang, J., Bao, H. & Li, H., "Design of multi-wing 3D chaotic
systems with only stable equilibria or no equilibrium point using rotation
symmetry," *AEU – Int. J. Electron. Commun.* 135:153710, 2021. DOI
10.1016/j.aeue.2021.153710. Also Wang, X. & Chen, G., "Symmetrical multi-petal
chaotic attractors in a 3D autonomous system with only one stable equilibrium,"
*Proc. 4th Int. Workshop on Chaos-Fractals Theories and Applications*, 82–85, 2011.
DOI 10.1109/IWCFTA.2011.11.** What they did: modern instances of the cover
technique — rotation symmetry used as a *construction tool* to multiply wings/petals
of hidden attractors, including around stable equilibria. How we differ: as above
(top-down covers vs bottom-up search in a fixed smooth family); these papers also
show the technique is alive and standard, which strengthens the case that
equivariance alone cannot be our claim. Claims to reword: same as the Lü–Chen entry.

## 3. The "new chaotic system" genre and its standards

**Sprott, J. C., "Some simple chaotic flows," *Phys. Rev. E* 50(2):R647–R650, 1994.
DOI 10.1103/PhysRevE.50.R647.** What he did: computer search over millions of
randomly generated quadratic ODE systems, yielding the 19 minimal chaotic flows
(Sprott A–S); the ancestor of every automated hunt for chaotic ODEs, ours included.
How we differ: Sprott searched for *algebraic minimality* with no structural prior;
we search inside a symmetry-constrained ansatz, score *shape novelty against a
catalog*, and keep a diversity archive rather than a list of hits. Claims to reword:
none — but the paper must name this as the ancestor (the plan already says so).

**Sprott, J. C., "A proposed standard for the publication of new chaotic systems,"
*Int. J. Bifurcation Chaos* 21(9):2391–2394, 2011. DOI 10.1142/S021812741103009X.**
What he did: editorial criteria meant to stem the flood of cosmetically "new"
chaotic systems: a new system should be simpler than anything known with the same
behavior, or exhibit behavior of genuine scientific interest, not merely differ in
coefficients. Consequence for us: this is exactly why the paper is a **methods
paper** — the three attractors are exhibits of the search method, and the paper
must self-apply Sprott's standard when presenting them (we claim measured shape
novelty and a parametric equivariant family, not "yet another chaotic system").
Claims to reword: none in README/STORY; governs the paper's framing.

**One-equilibrium chaotic flows.** **Sprott, J. C., "Simplest dissipative chaotic
flow," *Phys. Lett. A* 228(4–5):271–274, 1997. DOI 10.1016/S0375-9601(97)00088-1**
(the minimal jerk flow, with a single equilibrium); **Wang, X. & Chen, G., "A
chaotic system with only one stable equilibrium," *Commun. Nonlinear Sci. Numer.
Simul.* 17(3):1264–1272, 2012. DOI 10.1016/j.cnsns.2011.07.017**; **Molaie, M.,
Jafari, S., Sprott, J. C. & Golpayegani, S. M. R. H., "Simple chaotic flows with one
stable equilibrium," *Int. J. Bifurcation Chaos* 23(11):1350188, 2013. DOI
10.1142/S0218127413501885** (twenty-three of them, in the hidden-attractor
program). What they did: established that chaotic flows with exactly one
equilibrium — even a stable one — are a recognized genre. How we differ: our point
is not the equilibrium *count* but the geometry: a single saddle-focus with a 2-D
unstable spiral manifold organizing the entire symmetric recirculation (Shilnikov
anatomy), with the symmetry axis through it. Claims to reword: STORY's "precisely
one equilibrium, which is unusual" — drop or qualify "unusual."

**Shilnikov, L. P., "A case of the existence of a countable number of periodic
motions," *Sov. Math. Dokl.* 6:163–166, 1965.** The saddle-focus homoclinic
theorem behind the "Shilnikov-type" language; gap item 4 (the homoclinic hunt)
exists to support this connection numerically. Claims to reword: none if the hunt
is reported honestly (near-homoclinic distances, not a proof).

## 4. Automated discovery of chaotic systems

**Sprott, J. C., "Automatic generation of strange attractors," *Computers &
Graphics* 17(3):325–332, 1993. DOI 10.1016/0097-8493(93)90082-K.** What he did:
fully automated generation of chaotic 2-D quadratic *maps*, filtered by Lyapunov
exponent and correlation dimension, explicitly in service of aesthetics — the
direct ancestor of "search for beautiful attractors," twenty-three years early.
How we differ: maps vs flows; symmetry ansatz; catalog-grounded novelty; QD
archive. Claims to reword: STORY's "whether the looking could be industrialized"
narrative should nod to this ancestry (one clause suffices).

**Goertzel, B., "Rapid generation of strange attractors with the eugenic genetic
algorithm," *Computers & Graphics* 19(1):151–156, 1995. DOI
10.1016/0097-8493(94)00130-Q.** GA acceleration of Sprott's map search. Same
differences.

**Pan, I. & Das, S., "When Darwin meets Lorenz: Evolving new chaotic attractors
through genetic programming," *Chaos, Solitons & Fractals* 76:141–155, 2015. DOI
10.1016/j.chaos.2015.03.017.** What they did: multi-gene genetic programming over
modifications of the Lorenz vector field, maximizing the largest Lyapunov exponent;
reported over a hundred evolved chaotic attractor variants. The closest published
ancestor for *evolving new chaotic flows*. How we differ: (i) single-objective LLE
maximization vs quality-diversity (their own figures show many Lorenz-like twins;
QD's cells + novelty exist precisely to avoid that collapse); (ii) no quantitative
novelty measure against known systems; (iii) free-form GP trees vs a structured
equivariant family with interpretable parameters; (iv) no certification battery
(spectra, D_KY, equilibria, bifurcations) per discovered system. Claims to reword:
any "first evolutionary search for new chaotic flows" — we make no such claim; the
paper's related work must cite this.

**Lai, J., Bao, A. & Gilpin, W., "Panda: A pretrained forecast model for chaotic
dynamics," arXiv:2505.13755, 2025.** What they did: generated ~2×10⁴ novel chaotic
ODEs by evolutionary search seeded from the dysts catalog, as *training data* for a
forecasting foundation model (chaoticity-filtered, not individually studied). The
largest automated chaotic-system generation effort to date, from the same group as
dysts. How we differ: (i) goal — training-set augmentation vs discovery and
characterization of individually certified, presentable systems; (ii) selection —
chaoticity filters vs novelty-against-the-catalog as the explicit objective, with a
MAP-Elites archive organized by symmetry order and dynamical descriptors; (iii)
structure — unconstrained recombinations vs an interpretable equivariant
normal-form ansatz; (iv) output — a corpus vs named systems with full spectra,
equilibrium analysis, bifurcation scans, and error bars. Claims to reword: the
engine must not be described as the first or largest automated generator of chaotic
systems; its claim is the QD + catalog-grounded-novelty + symmetry combination and
the certification standard.

## 5. Catalogs, benchmarks, fingerprints

**Gilpin, W., "Chaos as an interpretable benchmark for forecasting and data-driven
modelling," *NeurIPS Datasets and Benchmarks*, 2021. arXiv:2110.05266.** What he
did: dysts — the machine-readable catalog (131 systems at publication; 135 in the
version we fingerprint) with standardized integration and annotated invariants;
also embedded the catalog in feature space to map its diversity. Relevance: dysts
is our novelty reference set, and Gilpin's feature-space embedding is the closest
precedent for our fingerprint map figure (different features: forecasting-oriented
statistics vs geometry-only shape descriptors). Claims to reword: none; cite both
roles.

**Osada, R., Funkhouser, T., Chazelle, B. & Dobkin, D., "Shape distributions,"
*ACM Trans. Graphics* 21(4):807–832, 2002.** The D2 pairwise-distance histogram
our fingerprint uses (plus PCA eigenvalue ratios and a fill factor). Claims to
reword: none; provenance citation.

**Langford, W. F., "Numerical studies of torus bifurcations," in *Numerical Methods
for Bifurcation Problems*, ISNM 70, Birkhäuser, 285–295, 1984. DOI
10.1007/978-3-0348-6256-1_19.** The system that circulates in catalogs (including
dysts) as the "Aizawa attractor" — Aurelia's nearest fingerprint neighbor — traces
to Langford's torus-bifurcation studies. Cite when naming the neighbor.

## 6. Quality-diversity search

**Lehman, J. & Stanley, K. O., "Abandoning objectives: Evolution through the search
for novelty alone," *Evolutionary Computation* 19(2):189–223, 2011. DOI
10.1162/evco_a_00025. Mouret, J.-B. & Clune, J., "Illuminating search spaces by
mapping elites," arXiv:1504.04909, 2015. Pugh, J. K., Soros, L. B. & Stanley, K. O.,
"Quality diversity: A new frontier for evolutionary computation," *Front. Robot. AI*
3:40, 2016. DOI 10.3389/frobt.2016.00040.** What they did: novelty search,
MAP-Elites, and the QD framing, in evolutionary robotics and beyond. Relevance: the
engine is a textbook MAP-Elites with novelty-as-quality; the behavior descriptors
(symmetry order × λ₁ × aspect) and the novelty signal (fingerprint distance to a
*scientific catalog* rather than to the population) are the domain transfer.
Searches on June 11, 2026 for QD/novelty-search applied to discovering chaotic
attractors or dynamical systems found no prior work; the nearest neighbor is the
intrinsically-motivated pattern-discovery line in self-organizing systems (e.g.,
Reinke et al., "Intrinsically motivated discovery of diverse patterns in
self-organizing systems," arXiv:1908.06663), which targets cellular automata, not
ODE attractors. Claims to reword: state the method claim as "to our knowledge, the
first quality-diversity search over a family of dynamical systems with novelty
measured against a catalog of known attractors" — and keep the "to our knowledge."

## 7. Claim-by-claim verdict

| Claim (where) | Verdict | Action |
|---|---|---|
| "No published 3-D flow with discrete C_n equivariance [through this coupling]" (STORY §The name; README §Novelty's framing of relatives as only F&G + Aizawa) | **Fails** in broad form: proto-Lorenz covers, Letellier–Gilmore covers, Thomas 1999, rotation-symmetric multi-scroll | Reword to: smooth globally polynomial flow, equivariant by construction (resonant coupling), *not a cover of a known system*, single organizing saddle-focus, found by search — and cite the counterexample literature |
| F&G icons are 2-D discrete maps, not flows (README, STORY) | Survives | Keep |
| Aizawa/Langford is SO(2)-symmetric with different equilibrium structure; nearest catalog neighbor at 0.185 (README) | Survives (measured) | Keep; cite Langford 1984 as the source of the "Aizawa" system |
| "Exactly one equilibrium, which is unusual" (STORY §Anatomy) | **Fails as stated**: one-equilibrium chaotic flows are an established genre | Drop "unusual"; emphasize the organizing saddle-focus + symmetry-axis geometry instead |
| Names Aurelia/Naiad/Cassiopea unclaimed (README, STORY) | Survives (re-checked 2026-06-11) | Keep |
| Novelty metric: nearest-catalog-distance with D2+PCA fingerprints; Aurelia 0.185 / Naiad 0.34 / Cassiopea 0.224 vs median NN spacing 0.115 (README) | Survives (it is a measurement; method validated in gap 3) | Keep; cite Osada 2002 + Gilpin 2021 |
| The QD engine is novel (README §atlas engine) | Survives **narrowly**: no QD/novelty-search precedent found for attractor discovery, but evolutionary generation of chaotic flows exists (Pan & Das 2015) and at scale (Panda 2025) | Add positioning sentence + "to our knowledge"; never claim "first automated/evolutionary discovery of chaotic systems" |
| "Found by search ... judged by eye for beauty" (README §Discovery method) | Survives, with ancestry | Nod to Sprott 1993 automated aesthetic search |
| w̄^(n−1) is "the unique term" surviving the rotation (README, STORY, module docstrings) | Survives as mathematics, but is textbook equivariant theory | Present as the known Z_n-resonant normal-form coupling (cite GSS Vol. II), not as a discovery |

## What the paper must do (summary for §2 Related Work)

1. Open with Sprott 1994/2011: the genre and its standard; declare the methods-paper
   framing.
2. Symmetric chaos: F&G (maps, aesthetics), GSS (the resonant coupling), Melbourne/
   Ashwin (attractor symmetry vocabulary).
3. Equivariant 3-D flows exist: Miranda–Stone, Letellier–Gilmore (covers), Thomas
   (cyclic), multi-scroll (engineered). State precisely what is new here: the family
   is *not a cover*, is globally polynomial, parameterized by symmetry order n, and
   was *discovered* by search rather than constructed.
4. Automated discovery: Sprott 1993 → Goertzel 1995 → Pan & Das 2015 → Panda 2025;
   position the engine as QD + catalog-grounded novelty + certification.
5. QD lineage: Lehman–Stanley, Mouret–Clune, Pugh et al.; the transfer is the
   novelty-vs-scientific-catalog signal and symmetry-indexed cells.


---

# Part B — Phase-0 polyhedral, rotoreflection & chaos-law verification

Verification of the three novelty ("virginity") claims underpinning the symmetry program
(docs/SYMMETRY_PROGRAM.md). Conducted 2026-06-11 via targeted web search (multiple engines),
direct reading of primary-source PDFs where accessible, and follow-up fetches of abstracts.
Scope caveat: this is a careful but finite search, not a systematic review; verdicts are
"no prior art found", not "no prior art exists". Every claim below distinguishes MAPS
(discrete iterations) from FLOWS (autonomous ODEs), which is where most of the apparent
prior art dissolves.

---

## Claim 1 — Polyhedral chaotic flows

**Question.** Has anyone exhibited a chaotic autonomous 3-D ODE flow whose strange attractor
has tetrahedral (T), octahedral (O), or icosahedral (I) point-group symmetry?

**VERDICT: NO PRIOR EXHIBIT FOUND — the land is clear.** Everything that looks like prior
art is either (a) an iterated map, (b) a flow whose attractor is not chaotic, (c) a chaotic
flow in dimension 4, or (d) abstract admissibility theory with no constructed example.

### (a) The map literature (polyhedral symmetry, but discrete iterations)

- G. F. Brisson, K. M. Gartz, B. J. McCune, K. P. O'Brien, C. A. Reiter,
  "Symmetric attractors in three-dimensional space", *Chaos, Solitons & Fractals* **7**(7),
  1033–1051 (1996). Iterated equivariant polynomial maps with cube symmetry.
- C. A. Reiter, "Chaotic attractors with the symmetry of the tetrahedron",
  *Computers & Graphics* **21**(6), 841–848 (1997). Iterated maps.
- C. A. Reiter, "Chaotic attractors with the symmetry of the dodecahedron" (preprint,
  Lafayette College, NSF DMS-9805507, https://webbox.lafayette.edu/~reiterc/dodec/dodec02_pp.pdf).
  **Read directly for this verification**: the construction is f(x) = Σ_{σ∈G} σ⁻¹P(σx) for a
  random polynomial P — a group-averaged iterated MAP, Monte-Carlo searched over coefficients,
  with a Lyapunov-exponent aesthetic gate (0.01–0.6). Icosahedral rotation group I (order 60)
  and full group Ī (order 120). Explicitly a map; no differential equations anywhere.
- K. C. Jones, C. A. Reiter, "Chaotic attractors with cyclic symmetry revisited",
  *Computers & Graphics* **24**(2), 271–282 (2000). Maps.
- M. Field, M. Golubitsky, *Symmetry in Chaos: A Search for Pattern in Mathematics, Art and
  Nature*, OUP 1992 / 2nd ed. SIAM 2009. Symmetric icons and quilts — planar iterated maps.

### (b) Flows with related symmetry whose attractors are not chaotic

- J. Guckenheimer, P. Holmes, "Structurally stable heteroclinic cycles",
  *Math. Proc. Camb. Phil. Soc.* **103**, 189–192 (1988). The classic R³ ODE
  ẋᵢ = (λ + a₁xᵢ² + a₂xᵢ₊₁² + a₃xᵢ₊₂²)xᵢ, equivariant under the order-24 group generated by
  cyclic coordinate permutation and rotation-type sign changes (a rotation subgroup of the
  cube's symmetry group). The attractor is a structurally stable heteroclinic cycle — not chaotic.
- I. S. Labouriau, A. C. Murza, "Hopf bifurcation with tetrahedral and octahedral symmetry",
  arXiv:1407.2866 (2014). Periodic solutions via the equivariant Hopf theorem — not chaos.

### (c) Chaotic equivariant flows in the wrong dimension

- J. Guckenheimer, P. Worfolk, "Instant chaos", *Nonlinearity* **5**, 1211–1222 (1992).
  Verified via their survey (J. Guckenheimer, P. Worfolk, "Dynamical systems: some
  computational problems", arXiv:chao-dyn/9304010, §6): a stable equilibrium bifurcates
  directly to a small-amplitude chaotic attractor — but the vector field lives on **R⁴**
  (ẋᵢ = (λ + a₁xᵢ² + …)xᵢ ± e·xⱼxₖxₗ, equivariant under signed cyclic permutations of four
  coordinates, the orientation-preserving subgroup of the Field–Swift group), with the
  attractor on an invariant 3-sphere. Not a 3-D flow, and not a T/O/I point group of
  physical 3-space.

### (d) Admissibility theory (abstract; no constructed specimens)

- I. Melbourne, M. Dellnitz, M. Golubitsky, "The structure of symmetric attractors",
  *Arch. Rational Mech. Anal.* **123**, 75–98 (1993).
- P. Ashwin, I. Melbourne, "Symmetry groups of attractors",
  *Arch. Rational Mech. Anal.* **126**, 59–78 (1994). Admissible subgroups Σ ⊆ Γ classified;
  results illustrated for finite subgroups of O(3).
- M. Field, I. Melbourne, M. Nicol, "Symmetric attractors for diffeomorphisms and flows",
  *Proc. London Math. Soc.* (3) **72**, 657–696 (1996). **The key theoretical landmark**: for
  flows in dimension ≥ 3 with no reflection hyperplanes there are essentially no restrictions
  on which symmetry groups attractors can realize. T, O, I are admissible for flows in R³ —
  the theory says the specimens can exist; nobody appears to have built them.

**Consequence for the program.** Phase 2 (polyhedral specimens) keeps its headline claim:
*first explicitly constructed, certified chaotic 3-D flows whose attractors carry the
polyhedral rotation groups*. The map literature becomes the "previously, only maps" framing;
Field–Melbourne–Nicol admissibility is the theory the constructions instantiate.

---

## Claim 2 — Rotoreflection (improper) symmetry by design

**Question.** Has anyone designed a chaotic flow whose attractor is symmetric under a
rotoreflection (S₄/S₆ in Schoenflies notation — rotation composed with reflection) without
being symmetric under the rotation or the reflection separately?

**VERDICT: NO PRIOR ART FOUND FOR S₄/S₆; the order-2 improper group (inversion, S₂) IS
occupied.** The S₄ specimen of Phase 1 stands, with the claim stated precisely.

### What is occupied

- C. Letellier, R. Gilmore, "Covering dynamical systems: Twofold covers",
  *Phys. Rev. E* **63**, 016206 (2001). **Read directly for this verification.** They state
  there are exactly three order-2 groups in R³ — reflections σ_z, rotations R_z(π),
  inversions P: (x,y,z)→(−x,−y,−z) — and construct/catalog equivariant systems for R_z(π)
  (Lorenz, Burke–Shaw) and P (induced Lorenz, Kremliovsky). Inversion is an improper
  operation (S₂), so *improper symmetry per se* is occupied at order 2. They explicitly
  exclude reflection groups because reflection-symmetric strange attractors must be
  disconnected.
- R. Gilmore, C. Letellier, *The Symmetry of Chaos*, Oxford University Press (2007).
  Systematic n-fold covers with cyclic rotation symmetry C_n, plus the order-2 groups above.

### The precision boundary

The group S₄ = {E, S₄, C₂, S₄³} (order 4, cyclic, generated by a quarter-turn-plus-flip) is
**not** reachable by the occupied constructions: it contains no inversion, no reflection,
and its improper elements are genuine order-4 rotoreflections. Its only proper subgroup is
C₂. A designed flow whose attractor symmetry group is exactly S₄ — symmetric under the
quarter-turn-plus-flip but not under quarter-turn or flip alone — appears in no search we
ran ("improper symmetry attractor", "rotoreflection equivariant flow", "S4 point group
chaotic", and variants). Note also that rotoreflections fix no plane, so the
Letellier–Gilmore disconnectedness obstruction for reflections does not apply: a connected
S₄-symmetric attractor is permitted.

### Caveats and near-fields

- C. Letellier, R. Gilmore, "Symmetry groups for 3D dynamical systems",
  *J. Phys. A: Math. Theor.* **40**, 5597–5620 (2007). The general taxonomy paper. We could
  not access the full text during this strike; it plausibly *enumerates* improper groups as
  possibilities. The defensible claim is therefore about the **constructed, certified
  specimen**, not about the group being unconsidered.
- J. S. W. Lamb, J. A. G. Roberts, "Time-reversal symmetry in dynamical systems: a survey",
  *Physica D* **112**, 1–39 (1998). "Reversible-equivariant" systems combine spatial actions
  with time reversal — a conceptually distinct notion (our σ acts on space only). Cited to
  disambiguate, not as prior art.

**Consequence for the program.** Phase 1 proceeds with the novelty claim phrased as: *the
first designed chaotic flow whose attractor symmetry group is the order-4 rotoreflection
group S₄ ⊂ O(3) — an improper symmetry that is not an inversion and contains no reflection*.
Hedge against the 2007 taxonomy paper in any publication.

---

## Claim 3 — A chaos-vs-symmetry-group law

**Question.** Is there any systematic study of how maximal Lyapunov exponent (or any chaos
measure) varies across symmetry groups in matched families of equivariant flows?

**VERDICT: NOTHING FOUND. This appears genuinely unmapped** — searches surface only
adjacent-but-different work:

- Multi-scroll/multi-wing engineering: e.g., "Design of multi-wing chaotic systems with
  higher largest Lyapunov exponent", *Chaos, Solitons & Fractals* (2022). Scroll/wing *count*
  is a geometric feature, not a symmetry group; no matched-budget cross-group comparison.
- A. Biswas, S. S. Chaurasia, P. Parmananda, S. Sinha, "Asymmetry induced suppression of
  chaos", *Scientific Reports* **10** (2020). Coupled-oscillator ensembles driven to steady
  states by dissimilar external systems — about network heterogeneity, not point groups.
- C. Bick et al., "Chaos in symmetric phase oscillator networks", *Phys. Rev. Lett.* **107**,
  244101 (2011). Existence of chaos in S_N-symmetric networks; no group-by-group law.

One structural insight from the covers literature sharpens our Phase 3 claim: in the
Letellier–Gilmore cover/image construction, the cover is *locally identical* to its image,
so its Lyapunov exponents are independent of the covering group — a C₂ cover and a C₄ cover
of the same image have the same λ₁ by construction. Any nontrivial dependence of attainable
λ₁ on the group in *our* families (which are independently searched, not covers of a common
image) is therefore a statement about the equivariant design space itself, not an artifact —
and the cover construction provides the natural null hypothesis to state in the paper.

**Consequence for the program.** Phase 3 is the robust headline regardless of Phases 1–2:
the λ₁-vs-group map (C₂ 0.296, C₃ 0.233, C₄ 0.525 already in hand) appears to be the first
systematic measurement of its kind.

---

## Decision (per the Phase 0 decision rules)

All three claims survive verification. Therefore:

1. **Phase 1 (S₄ specimen)** proceeds with its novelty claim, phrased per the precision
   boundary above.
2. **Phase 2 (polyhedral specimens)** keeps the flagship headline.
3. **Phase 3 (the law)** remains the robust spine of the paper.

## Search log

Queries run (2026-06-11, WebSearch + targeted WebFetch of PDFs/abstracts; Perplexity
unavailable — quota): "tetrahedral symmetry chaotic attractor flow ODE icosahedral
octahedral"; Brisson tetrahedron title search; "rotoreflection improper symmetry strange
attractor equivariant flow S4 point group"; "does symmetry suppress chaos" + Lyapunov
variants; Letellier–Gilmore symmetry-groups/covers searches; Ashwin–Melbourne admissibility;
Field "Dynamics and Symmetry" + icosahedral equivariant; Guckenheimer–Worfolk "instant
chaos" (×2); multi-scroll polyhedral; "number of scrolls"/symmetry-order vs largest Lyapunov
exponent; Cn-equivariant λ₁ comparisons. Primary sources read in full or in relevant part:
Reiter dodecahedron preprint (PDF), Letellier–Gilmore PRE 63 016206 (PDF, pp. 1–3),
Guckenheimer–Worfolk survey (arXiv chao-dyn/9304010, §6 text extraction).

---

## Claim 4: Dihedral (chiral D_n) chaotic flows (Trinacria check, 2026-10-21)

**Question.** Has anyone exhibited a chaotic autonomous 3-D flow whose attractor carries the
chiral dihedral point group D3 (rotation by 120 degrees about z plus half-turns about three
horizontal axes, all det +1, no mirror and no inversion)? More broadly, is D_n occupied for
flows?

**VERDICT: no explicit D3 exhibit found, but the general construction is not ours to
claim.** We found no 3-D flow with a certified chaotic attractor carrying exactly chiral D3.
Dihedral chaos is classical for maps, cover constructions make D_n-equivariant flows
available in principle, and D2 is plausibly occupied by cross-product quadratic systems. The
claim is therefore about one constructed and certified specimen, not about the group.

### (a) Dihedral symmetry in maps (classical)

- P. Chossat, M. Golubitsky, "Symmetry-increasing bifurcation of chaotic attractors",
  *Physica D* **32**, 423-436 (1988). Planar maps with D_m symmetry; the Chossat-Golubitsky
  map is the standard example.
- M. Field, M. Golubitsky, *Symmetry in Chaos*, OUP 1992 / SIAM 2009. Symmetric icons:
  planar iterated maps with dihedral and cyclic symmetry.
- I. Melbourne, M. Dellnitz, M. Golubitsky, "The structure of symmetric attractors",
  *Arch. Rational Mech. Anal.* **123**, 75-98 (1993). Admissible attractor symmetries for
  D_m-equivariant planar maps (for example, no D_k with 2 < k < m).
- C. A. Reiter and coauthors (Lafayette College): cyclic, dihedral and polyhedral
  equivariant polynomial MAPS in two and three dimensions, Monte-Carlo searched (see Claim 1).

None of these is a flow, so none occupies D3 for ODEs.

### (b) Cover constructions for flows (the main hedge)

- C. Letellier, R. Gilmore, "Symmetry groups for 3D dynamical systems",
  *J. Phys. A: Math. Theor.* **40**, 5597-5620 (2007), and R. Gilmore, C. Letellier,
  *The Symmetry of Chaos*, OUP 2007. Bibliographic record confirmed; the full text was again
  not readable from this sandbox. They classify the symmetry groups available to 3-D
  systems and build equivariant flows by lifting an image system through a covering map.
  Lifting through a map invariant under a D_n action would produce D_n-equivariant chaotic
  flows, so such flows are constructible from any chaotic image system and may already appear
  as examples there. Any publication must cite this and must not claim the group.

### (c) Neighbouring symmetric 3-D flows that are not D3

- Cyclically symmetric (circulant) flows: R. Thomas, *Int. J. Bifurcation Chaos* **9**,
  1889 (1999); the Halvorsen system; K. Rajagopal, F. Nazarimehr, S. Jafari, J. C. Sprott and
  coauthors on circulant chaotic flows with coexisting attractors (*IJBC*, 2019 and 2025).
  The cyclic permutation (x, y, z) -> (y, z, x) is a 120-degree rotation about (1, 1, 1), so
  these carry C3; with odd nonlinearity (Thomas: sin) they also carry the inversion, giving
  S6. We checked Thomas' system numerically: the half-turn (x, y, z) -> (-y, -x, -z) about
  an axis perpendicular to (1, 1, 1) fails at 40/40 random states (it reverses the cyclic
  order), so Thomas is S6, not D3 or D3d.
- Three identical coupled cells (Golubitsky-Stewart coupled-cell theory) carry the abstract
  group S3, which is isomorphic to D3. Acting by coordinate permutations on R^3, though, the
  transpositions are mirrors (det -1), so the point group is C3v, not the chiral D3; with an
  odd vector field it grows to D3d, which contains D3 but is larger. We did not search this
  area systematically.
- Four-wing quadratic systems. G. Qi, G. Chen, M. A. van Wyk, B. J. van Wyk, Y. Zhang,
  "A four-wing chaotic attractor generated from a new 3-D quadratic autonomous system",
  *Chaos Solitons Fractals* **38**, 705-721 (2008): x' = a(y - x) + eyz, y' = cx + dy - xz,
  z' = -bz + xy. By direct substitution (checked numerically at 40 random states) it is
  equivariant under (x, y, z) -> (-x, -y, z) only: the half-turns about x and y both fail
  because a(y - x) and cx + dy mix x and y. So the Qi system is C2, not D2. This corrects
  the earlier working caveat that it might occupy D2. Systems whose only nonlinear terms are
  the cross products yz, xz, xy with a diagonal linear part (for example x' = ax - yz,
  y' = -by + xz, z' = -cz + xy) are D2-equivariant under all three coordinate half-turns,
  so D2 is plausibly occupied by such systems; we did not identify the first named one.
- The Guckenheimer-Holmes cycle (Claim 1 (b)) has signed cyclic permutations only and no
  transposition-type half-turn, so it carries no D3; it is not chaotic in any case.

### Allowed wording

- Allowed: "Trinacria is the first D3 specimen of this atlas"; "to our knowledge, the first
  designed and certified polynomial flow whose attractor carries exactly the chiral dihedral
  group D3", always followed by the Letellier-Gilmore hedge above.
- Forbidden: "first dihedral attractor", "first D_n-equivariant chaotic flow", or any claim
  that the group itself is new.
- Shape statements stay catalog-relative (results/trinacria_verification.json,
  novelty_vs_dysts and novelty_vs_atlas).

### Name check

"Trinacria attractor", "Trinacria chaotic system", "Trinacria dynamical system": no
dynamical-systems use found (2026-10-21, repeating the selection-stage search recorded in
results/d3_selection.json). Trinacria Heiberg 1863 is an extinct marine diatom genus
(Hemiaulaceae) with Paleocene-Eocene fossils from Denmark.

### Search log (Claim 4)

WebSearch, 2026-10-21 (WebFetch fails on DNS here, so verdicts rest on snippets and on our
own substitution checks): Letellier-Gilmore 2007 J. Phys. A 40 5597; "dihedral symmetry D3
equivariant three-dimensional chaotic flow attractor Lyapunov"; Qi-Chen 2008 four-wing
symmetry; Chossat-Golubitsky 1988 symmetry-increasing bifurcation dihedral; "dihedral
symmetric chaotic attractor ODE three dimensions D_n equivariant vector field"; Reiter
dodecahedron and cyclic-symmetry papers; Field-Golubitsky Symmetry in Chaos; Sprott and
coauthors on circulant flows; cross-product four-wing system x' = ax - yz; three coupled cells
with S3 symmetry; the three Trinacria name queries. An independent reviewer ran six further
searches with the same outcome (no 3-D D3 flow with a Lyapunov analysis surfaced).

