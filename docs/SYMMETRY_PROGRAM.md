# The Symmetry Program: a periodic table of strange attractors

A research program for this repository (github.com/Sakeeb91/attractor-atlas). Read README.md and
STORY.md first for project state; docs/PAPER_PLAN.md holds the earlier methods-paper plan, which
this program SUPERSEDES AND ABSORBS: the paper's science question becomes the one below.

## The program in one sentence

Complete the map of which finite symmetry groups of 3-D space host minimal strange attractors:
explicit, certified, named, rendered specimens per group, plus an empirical law of how chaos
strength depends on the group. The existing family (Naiad C2, Aurelia C3, Cassiopea C4, all
certified, all in the live viewer) provides the first three entries.

## Why this and not something else (strategic context, already decided)

- "New chaotic system" papers are a saturated genre (Sprott's 2011 IJBC standard exists because
  of it). Plain Cn rotational symmetry is occupied territory: Letellier & Gilmore systematically
  build n-fold symmetric "covers" of known attractors (*The Symmetry of Chaos*, 2007).
- The ABSTRACT theory of which groups can be symmetry groups of attractors exists
  (Ashwin & Melbourne, 1990s). What does not exist, as far as this session could tell:
  the constructed zoo. Explicit minimal flows per point group, certified chaotic, compared
  quantitatively, rendered beautifully.
- Already-opened question we own: chaos strength is NOT monotonic in symmetry order
  (lambda1: C2 0.296, C3 0.233, C4 0.525). Nobody appears to have mapped this systematically.

## Phase 0: the literature verification strike (DO THIS FIRST, ~1-2 hours)

"Virgin" is a hypothesis until checked. Verify exactly three claims with targeted searches
(WebSearch, alphaxiv-paper-lookup skill, Google Scholar queries via WebFetch):

1. **Polyhedral flows**: has anyone exhibited a chaotic FLOW (autonomous 3-D ODE) whose attractor
   has tetrahedral/octahedral/icosahedral symmetry? Known near-misses to check and cite:
   Brisson et al. "Chaotic attractors with the symmetry of a tetrahedron" (Computers & Graphics
   1997?, likely MAPS not flows); Reiter & Carter symmetric-chaos map art; Field & Golubitsky
   *Symmetry in Chaos* (maps); Ashwin/Melbourne admissibility theory (abstract);
   Golubitsky-Stewart-Schaeffer *Singularities and Groups in Bifurcation Theory* Vol II
   (equivariant vector field tables, not chaos); Melbourne work on attractors with spherical
   symmetry. Search terms: "icosahedral chaotic attractor flow", "octahedral symmetry strange
   attractor ODE", "polyhedral equivariant chaotic flow", "attractor admissible symmetry group".
2. **Rotoreflection (improper) symmetry**: any designed attractor symmetric under a
   rotation-plus-reflection (S4/S6 point groups, Schoenflies notation) but not the pure rotation?
   Search: "improper symmetry attractor", "rotoreflection equivariant flow", "S4 point group
   chaotic", "reversible equivariant chaotic attractor".
3. **Chaos-vs-symmetry-group law**: any systematic study of maximal Lyapunov exponent across
   symmetry groups in equivariant families? Search: "symmetry order Lyapunov exponent",
   "does symmetry suppress chaos", "equivariant chaos systematic".

Decision rules: if (1) has prior art, the polyhedral specimens become "additional examples" and
the program's headline shifts to the comparative law (3), which is robust. If (2) has prior art,
skip the S4 specimen as a novelty claim but keep it as a family member. If (3) has prior art
(unlikely), the program reduces to the bestiary + paper. Write findings with full citations to
docs/RELATED_WORK.md (create or extend) and STATE THE VERDICT for each claim. Commit before
building anything.

## Phase 1: the S4 rotoreflection attractor (first new specimen, smallest step)

The symmetry: sigma: (w, z) -> (i*w, -z) where w = x+iy. Quarter-turn PLUS vertical flip; the
attractor is its own upside-down quarter-turn but NOT symmetric under quarter-turn alone or flip
alone. sigma^2 = (w,z)->(-w,z) so the group contains the half-turn C2; group is S4 (order 4).

Equivariance conditions for dw/dt = f(w, wbar, z), dz/dt = h(w, wbar, z):
  f(i*w, -i*wbar, -z) = i * f(w, wbar, z)        (note: conj(i*w) = -i*wbar)
  h(i*w, -i*wbar, -z) = - h(w, wbar, z)
Term audit (verify each before trusting; write a numerical equivariance unit test):
  ALLOWED in f: w (with z-EVEN coefficient, e.g. (alpha*z^2 - beta + i*omega)*w);
                wbar^3 (since (-i)^3 = i);
                z*wbar (the SIGNATURE term: allowed under S4, FORBIDDEN under pure C4)
  FORBIDDEN in f: z*w (this is why the existing Aurelia-family term a(z-b)w breaks the flip)
  ALLOWED in h: z, z^3, z*|w|^2, Re(w^2), Im(w^2) (signature terms: odd under sigma),
                z*Re(w^4), z*Im(w^4)
  FORBIDDEN in h: constants, |w|^2 alone, Re(w^4)/Im(w^4) alone (all even under sigma)
Candidate family (7 dials, mirrors the existing family's size):
  dw/dt = (alpha*z^2 - beta + i*omega) w + gamma*wbar^3 + delta*z*wbar
  dz/dt = nu*z - z^3 - lam*z*|w|^2 + eps*Im(w^2)
Procedure: add this family to atlas/ (new module or a family registry), FIRST write the
equivariance test (apply sigma to random states, check F(sigma s) = sigma' F(s) to 1e-12),
then Monte-Carlo/MAP-Elites search for bounded chaos exactly as scripts/run_atlas.py does
(reuse the evaluator pattern in atlas/family.py: batched complex RK4, two-orbit lambda1, gates).
Aesthetic target per project convention: lambda1 > ~0.2, D_KY in 2.1-2.3, layered veils.
Then full certification battery (mirror scripts/verify_naiad.py: long QR spectrum, IC robustness,
equilibria w/ classification, 4M-step boundedness, divergence, novelty vs
results/catalog_fingerprints.json), renderer with a NEW palette family (existing: Aurelia gold,
Naiad cyan, Cassiopea emerald; pick something distinct, e.g. violet/rose), gallery, README
section, STORY chapter, viewer entry (index.html SYSTEMS registry: switcher and handlers are
table-driven, adding an entry is enough; .seg grid-template-columns needs repeat(4,1fr)).
Naming convention: marine lineage (Aurelia=moon jelly, Naiad=water nymph, Cassiopea=upside-down
jelly). For S4 consider sea creatures with up-down ambiguity; verify the chosen name is unclaimed
("<name> attractor" search).

## Phase 2: polyhedral specimens (the flagship)

Target: first certified strange attractor of a designed flow with tetrahedral (T, order 12),
then octahedral (O, 24), then icosahedral (I, 60) ROTATION symmetry. Construction recipe
(verify equivariance numerically; consult Golubitsky-Stewart-Schaeffer Vol II tables if needed):
  - T-equivariant building blocks: linear x; gradient of invariants: grad(xyz) = (yz, zx, xy);
    grad(x^2+y^2+z^2) = 2x; cubic (x^3, y^3, z^3) relates to O; radial terms |x|^2 * x.
  - General dissipative ansatz:
      dx/dt = R x + a * grad(I_G)(x) + b * (x cross grad(I_G)(x)) - (c + d*|x|^2) x
    where I_G is the group's lowest non-trivial invariant (T: xyz; O: x^4+y^4+z^4;
    I: the degree-6 icosahedral invariant, look up its explicit form) and R is a G-commuting
    linear map (for T/O/I irreducible on R^3, R = scalar; rotation must come from the cross term).
  - CAUTION: for T/O/I the linear part commuting with G is only c*Identity (irreducibility), so
    spin enters through nonlinear terms; expect different dynamics than the Cn family. If the
    simple ansatz refuses to go chaotic, widen: two invariants, higher equivariants, or couple a
    second radial variable.
  - Search each family with the QD engine (extend descriptor cells: group label replaces n).
  - Icosahedral render should be the program's signature image.

## Phase 3: the chaos-vs-symmetry law (the finding)

Matched-budget searches per group: Cn for n=2..8, the S4/S6 improper groups, T, O, I (and Dn
dihedral if a Dn-equivariant family is constructed: needs an additional anti-symmetry in the
plane, e.g. invariance under w -> wbar). For each group: same evaluator settings, same number of
candidate evaluations (e.g. 20k), record the maximum certified lambda1 and D_KY achieved and the
fraction of parameter space that is chaotic. Output: results/symmetry_law.csv + one figure
(lambda1_max vs group, ordered by group order) + a section of analysis. This converts the
bestiary into a quantitative claim regardless of what Phase 0 finds.

## Phase 4 (optional): symmetry-increasing crisis in the Aurelia family

Hunt parameter regions where the C3 family hosts three conjugate symmetry-broken attractors that
merge into one symmetric attractor as a parameter crosses a crisis (compare attractor supports
from symmetric-related ICs; detect merging via fingerprint distance between the three). Known in
maps (Chossat-Golubitsky); thin for flows. One script + figure if found.

## Acceptance criteria and conventions

- One logical change = one atomic commit, pushed to github.com/Sakeeb91/attractor-atlas.
- Every certification claim has a script and a results/*.json; every render visually inspected
  (view the PNG; if 3-D matplotlib renders come out pure black, check numpy >= 2.3 FIRST,
  dysts installs downgrade it; histogram calls need explicit bin edges; np.ptp(arr) not arr.ptp()).
- Equivariance of any new family is unit-tested numerically before any search runs.
- Novelty is measured (fingerprint distance vs results/catalog_fingerprints.json) AND hedged
  (catalog is 135 systems, not the literature).
- Names: marine lineage, verified unclaimed.
- Final deliverable of the program: updated paper skeleton (docs/PAPER_PLAN.md section list)
  where Section 4 becomes the symmetry law + bestiary; STORY.md gains chapters as specimens land.
- Report at each phase boundary: verdicts (Phase 0), certified numbers (Phases 1-2), the law
  figure (Phase 3).

## Outcomes (as executed, June 2026)

**Phase 0 (literature strike): done.** All three novelty claims survived a targeted search;
verdicts in `docs/RELATED_WORK.md` (Part B). No prior polyhedral chaotic *flow* found in the
literature -- only maps (Field-Golubitsky, Reiter) and admissibility theory (Melbourne-Dellnitz-
Golubitsky). That absence turns out to be meaningful (see Phase 2).

**Phase 1 (S4 rotoreflection specimen): done.** Mobula, certified and named (lambda1 = 0.330,
D_KY = 2.266), in the viewer and the bestiary. The S4 group acts *reducibly* on R^3 (a 2-D plane
rotation plus a 1-D height), so it inherits the C_n family's saddle-focus engine and goes chaotic
readily.

**Phase 2 (polyhedral specimens): the irreducibility obstruction (a negative result that is the
finding).** The T, O, and I equivariant families are built and their symmetry is proven
numerically (`atlas/family_{t,o,i}.py`, `tests/test_{t,o,i}_equivariance.py`). But unlike C_n/S4,
the polyhedral rotation groups act *irreducibly* on R^3, so by Schur the only commuting linear map
is rho*I: **there is no saddle-focus linear engine.** An exhaustive, dt-robust search found NO
strange attractor in these families:

  - Polynomial ansatze (degree-3 complete, plus degree-4/5 gradient and non-gradient rotational
    equivariants; both uniform-random and fixed-point-stability-targeted; ~10^5 evaluations):
    bounded orbits relax to equilibria or limit cycles. The "all fixed points unstable + bounded"
    configurations (forced by Poincare-Hopf to host a cycle or a strange set) yield limit cycles.
  - Coarse-dt Benettin estimates *fabricate* chaos around stable equilibria: a candidate with
    lambda1 ~ +0.4 at dt=0.01 collapses to lambda1 ~ -15 to -25 at dt=0.0025 (a stable fixed
    point). Every apparent hit reversed sign under dt refinement -- the standard re-certification
    caught all of them.
  - Trigonometric "labyrinth" T-equivariant flows (Thomas-style, the only chaotic-flow route in
    the literature with related symmetry) give at most weak (lambda1 ~ 0.04) chaos that is
    spatially *localized* to a single labyrinth cell -- i.e. symmetry-BROKEN, not a single
    T-invariant attractor. Near the origin (the only T-fixed point a symmetric attractor could
    center on) the trig flow reduces to its non-chaotic polynomial Taylor expansion.

This is strong numerical evidence -- not a proof -- that low-degree equivariant 3-D flows under an
*irreducible* polyhedral action resist chaos, and is a candidate explanation for the literature's
absence of any polyhedral chaotic flow. The achievable polyhedral invariant sets are equilibria
and limit cycles, which *can* carry the full rotation group; chaos cannot (here).

**Phase 3 (the chaos-vs-symmetry-group law): done -- and it now carries the Phase 2 finding.**
`scripts/symmetry_law.py` runs a matched-budget search over C2..C8, S4, T, O, I, re-certifies the
top-K of each with the long integrator, and records the maximum certified lambda1 per group:
`results/symmetry_law.csv`, figure `gallery/symmetry_law.png`. The law shows two regimes: the
*reducible-action* groups (C_n, S4) host chaos with a non-monotonic dependence on order, while the
*irreducible-action* polyhedral groups (T, O, I) sit at ~zero certified chaos -- the irreducibility
cliff. The robust comparable is max certified lambda1 under equal evaluation budget (the chaotic
*fraction* is only loosely comparable across families with different parameter counts).

## Outcome added October 2026: the dihedral family and Trinacria

**D_n family: built.** The dihedral family anticipated in Phase 3 exists as
`atlas/family_d.py`. The extra generator is the half-turn kappa:(w, z) -> (conj(w), -z),
not the planar mirror w -> conj(w) suggested above: the mirror has determinant -1, and the
half-turn keeps the group chiral. kappa forbids a uniform spin, so the rotation rate is odd in
z. `tests/test_d_equivariance.py` proves equivariance for n = 2 to 6 and excludes every larger
group (C_2n, z-flip, mirrors, S_2n, inversion, and the copies of O and I containing D_n).

**Trinacria (D3): certified, with a disclosed basin.** lambda1 = 0.33 +- 0.01 (ensemble mean
0.328 +- 0.002 at dt 0.00125, converged under dt halving), D_KY 2.219, group residuals at the
self floor and controls at least 15 times worse, a single symmetric attractor rather than a
kappa-conjugate pair. About half of [-2,2]^3 escapes to infinity; every bounded orbit lands on
the attractor. See the README section and `results/trinacria_verification.json`.

**Law: D3 row added.** At the matched budget D3 reaches a maximum certified lambda1 of 0.689
(D_KY 2.49), the largest of any group so far, against 0.272 for C6, the cyclic group of the
same order (`results/law/law_D3.json`, run with PYTHONHASHSEED=0).
