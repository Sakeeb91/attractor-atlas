# Hand-off: Symmetry Program Phases 2 & 3

> **Status (2026-10-22): historical, except sections 1 and 2.** Phases 2 and 3 were executed
> in June 2026, and the D3 law row and the Trinacria specimen (chiral dihedral D3) were added
> in October 2026. Read this file as follows:
> - **Section 0 is obsolete.** The branch divergence it describes was resolved by the merge
>   in commit `bad965d`; the `symmetry-program` branch and the separate paper worktree path no
>   longer need any action.
> - **Sections 1 and 2 are still the living workflow conventions.** New specimen work should
>   mirror them (equivariance test first, atomic commits, certify before naming, view every
>   render, the per-specimen file layout).
> - **Sections 3 and 4 are superseded** by [docs/METHODOLOGY.md](METHODOLOGY.md) sections 10
>   (the polyhedral irreducibility obstruction: no certified T/O/I specimen exists) and 11 (the
>   chaos-vs-symmetry-group law), and by the "Outcomes" section of
>   [docs/SYMMETRY_PROGRAM.md](SYMMETRY_PROGRAM.md). Do not start building the polyhedral
>   specimens described there.

Continue the symmetry program in `github.com/Sakeeb91/attractor-atlas` (local: "Attractor Origin").
Read, in order: `docs/SYMMETRY_PROGRAM.md` (the master plan — Phases 2/3 are §"Phase 2" and
§"Phase 3"), `docs/RELATED_WORK.md` (Phase-0 verdicts), `README.md` and `STORY.md` (project state),
then this file (the operational recipe and the math already worked out for you).

Phases 0 (literature strike) and 1 (Mobula, the S4 rotoreflection specimen) are **done, certified,
and pushed** on branch `symmetry-program`. Your job: **Phase 2** (polyhedral T/O/I specimens — the
flagship) and **Phase 3** (the chaos-vs-symmetry-group law — the finding).

---

## 0. Branch and merge situation — READ FIRST

Two workstreams ran concurrently and the branches have **diverged**:

- `origin/main` is **13 commits ahead** on the *paper* side: a literature-positioning survey,
  Lyapunov convergence study (gap 2), fingerprint validation with a **validated novelty noise-floor**
  and softened/validated novelty numbers (gaps 3, `f6c44b8`/`f9e324f`), a Shilnikov homoclinic hunt
  (gap 4), and a **compiled paper draft** in `paper/main.tex`. The paper builds in a separate
  worktree at `/Users/sakeeb/Code repositories/attractor-atlas-paper` (checked out on `main`,
  `tectonic main.tex`).
- `symmetry-program` is **9 commits ahead** on the *program* side: Phase 0 + Mobula (Phase 1).
- **Both branches define `docs/RELATED_WORK.md` with different content** (main = paper positioning;
  symmetry-program = Phase-0 polyhedral verification). This is a known collision.

**Recommended first move:** merge `origin/main` into `symmetry-program` *before* building Phase 2, so
you inherit the paper draft (Phase 3 must update its §4), the validated novelty noise-floor (use it
for all novelty/law claims), and the convergence/fingerprint infrastructure. Resolve the
`docs/RELATED_WORK.md` conflict by **unioning the two as separate top-level parts** (Part A: paper
positioning; Part B: Phase-0 polyhedral/rotoreflection verification). After merging, re-check the
Mobula README novelty phrasing against main's *validated ranges* (main revised some catalog
distances; my Mobula novelty 0.273 was computed fresh and should be consistent, but verify).
Run `git status -sb` before every commit — the main checkout has been flipped between branches
before. `pytest tests/` must stay green through the merge.

If you prefer not to merge yet, you can build Phase 2/3 on `symmetry-program` and reconcile at the
end — but the divergence only grows. Merging early is cleaner.

---

## 1. Non-negotiable workflow (how Phase 1 was done; mirror it exactly)

1. **Equivariance test FIRST (TDD).** Before searching any new family, write a numerical
   equivariance unit test and watch it pass to ~1e-12. Pattern: `tests/test_s4_equivariance.py`.
   For a group G with generators {g_i}, assert `F(g·s) == g·F(s)` for random states/params and every
   generator (cheap to also loop the full generated group). Also assert the flow is **NOT**
   equivariant under any *larger* group you must exclude (for S4 we excluded pure C4 and pure flip;
   for T you must exclude the reflections of T_h, i.e. confirm `xyz` not `|x|`-type terms; for the
   law's honesty, confirm a T-family is not accidentally O-symmetric, etc.).
2. **One logical change = one atomic commit, pushed.** Nine commits for Mobula; see
   `git log --oneline cf879fa..9836a17`. Commit messages state the certified numbers.
3. **Beauty is a hard acceptance criterion.** Reject mathematically-valid-but-ugly candidates. Target
   λ1 > ~0.2, D_KY in 2.1–2.3, "layered veils." **Always view the rendered PNG** (Read the image);
   if a 3-D matplotlib render comes out pure black, check `numpy >= 2.3` FIRST.
4. **Novelty is measured AND hedged** (fingerprint distance vs `results/catalog_fingerprints.json`,
   135 systems, not the literature). Use main's validated noise-floor once merged.
5. **Certify before naming.** Run the full battery (`scripts/verify_mobula.py` is the template) and
   write `results/<name>_verification.json` before you commit a name.
6. **Names: marine lineage, verified unclaimed** (search `"<name> attractor"`). Used so far:
   Aurelia (moon jelly, C3), Naiad (water nymph, C2), Cassiopea (upside-down jelly, C4),
   Mobula (devil ray, S4), Trinacria (fossil marine diatom with a triskelion-shaped frustule,
   chiral dihedral D3). For polyhedra, think creatures/forms with the right symmetry, e.g.
   radiolaria/`Circogonia` (icosahedral!), `Lithocubus`/cubic radiolarians (octahedral),
   `Dorataspis` or a tetrahedral radiolarian (tetrahedral). Radiolaria are the obvious lineage:
   single-celled marine organisms whose silica skeletons realize exactly these point groups.
   *Circogonia icosahedra* is a real icosahedral radiolarian — a gift for the flagship I-specimen.
7. **Report at each phase boundary** (Phase 2: certified numbers per specimen; Phase 3: the law figure).

### Environment landmines (cost real time — do not re-learn)
- `python3` is homebrew 3.14; numpy is pinned to **2.3.5** (a dysts install can downgrade it to 2.1.x,
  which (a) renders 3-D matplotlib scatter pure black and (b) has a histogram fast-path bug — fixed in
  `atlas/fingerprint.py` via explicit bin edges). If renders go black, check numpy version first.
- Use `np.ptp(arr)`, not `arr.ptp()`. pytest 9.0.2, scipy 1.17.0 present.
- **browser-harness coordinate gotcha (discovered in Phase 1):** `click_at_xy` takes **page pixels**,
  but `capture_screenshot` images are ~2.5× downscaled. Read the element rect via `js("...
  getBoundingClientRect() ...")` and click those coordinates; don't click screenshot-pixel coords.
  `js()` eval is intermittently flaky in Comet — retry / use `ensure_real_tab()`. Start with
  `new_tab()`, never `goto_url()`. Serve the viewer over http (`python3 -m http.server`) — `file://`
  breaks the ES-module importmap.
- Front-end work (the viewer) requires the skill trio impeccable + taste-skill + ui-ux-pro-max
  (global CLAUDE.md). The viewer's `SYSTEMS` registry and `.seg` switcher are **table-driven**: add an
  entry + a `<button id="sys-<key>">` + widen `.seg` `grid-template-columns`. For a 5th system the
  segmented control will need `repeat(5,1fr)` and likely a smaller panel/typography (it's already at
  4 columns / 312px panel / 11px buttons after Mobula).

---

## 2. The codebase patterns to mirror (file-by-file, using Mobula as the worked example)

Per specimen `<name>` (a marine name), create a package mirroring `mobula/`:
- `<name>/dynamics.py` — constants, `velocity(state, *params)`, `jacobian`, `divergence`,
  `rk4_step`, `trajectory`. Docstring derives the real form from the complex/vector form.
- `<name>/equilibria.py` — `find_equilibria` (analytic on-axis roots + multi-start `fsolve`),
  `classify` (eigs + stability label).
- `<name>/lyapunov.py` — `lyapunov_spectrum` (Benettin QR), `kaplan_yorke_dimension`.
- `<name>/__init__.py` — re-export.
- `tests/test_<name>_equivariance.py` — the equivariance gate (write FIRST).
- `tests/test_<name>_matches_family.py` — cross-check the real-coordinate package equals the
  search-family evaluator at canonical params, + Jacobian-vs-finite-difference + divergence-vs-trace.
  (Mobula's analogue caught a hand-derivation slip risk; keep it.)
- `scripts/verify_<name>.py` — the battery; writes `results/<name>_verification.json`. For polyhedra,
  add a **symmetry residual** like Mobula's `symmetry_residuals()`: nearest-neighbour RMS of
  `g(orbit)` vs the orbit for each generator g, with a control that is NOT a symmetry, to prove the
  symmetry is genuine.
- `scripts/render_<name>.py` — gallery with a NEW palette (used: Aurelia gold, Naiad cyan, Cassiopea
  emerald, Mobula violet-rose; pick distinct — e.g. icosahedral = iridescent/opal or deep-ocean
  bioluminescent blue-green; the I render is the program's signature image).
- Viewer entry in `index.html` (`SYSTEMS` registry + button + grid width).
- README section + STORY chapter + `docs/PAPER_PLAN.md` §4 update.

The search family (the batched evaluator) is separate from the specimen package. For Cn it's
`atlas/family.py` (parameterized by `n`); for S4 it's `atlas/family_s4.py`. For each polyhedral group
write `atlas/family_<g>.py` with `PARAM_LO/HI`, `N_PARAMS`, `param_names()`, a single-state
`velocity(...)` (used by the equivariance test), a batched `_deriv/_rk4/evaluate_batch`, and
`sample_params/mutate_params`. The S4 module is the cleanest template — copy its structure.

The QD search itself: `scripts/run_atlas.py` loops `ORDERS=(2,3,4,5)` over `atlas/family.py`. The S4
search was run as standalone scratch scripts (`scratch/search_s4.py` → `refine_s4.py` →
`render_finalists.py`), which is fine — scratch/ is gitignored. **Lesson from Phase 1:** the batched
two-orbit Benettin estimate is *optimistic*; several S4 "alive" candidates **escaped under the long
fine integrator**. Always (a) widen the search box until the optimum is interior, then (b) re-certify
the top-K finalists with the real `trajectory()` integrator and drop escapees before naming.

---

## 3. Phase 2 — the polyhedral specimens (T → O → I, in that order)

**The flagship claim** (verified clear in `docs/RELATED_WORK.md`): the first *constructed, certified
chaotic 3-D flow* whose strange attractor carries a polyhedral **rotation** group. Do T (order 12)
first to validate the whole pipeline, then O (24), then I (60). The icosahedral render is the
signature image of the program.

### THE CENTRAL DIFFICULTY (internalize before coding)
For T, O, I the action on R³ is **irreducible**, so the only G-commuting linear map is a **scalar**.
There is **no saddle-focus engine** like the Cn family's `[a(z−b)+iω]w` linear core — the origin's
linearization is just `ρ·I`. **All chaos must come from the nonlinear terms; rotation enters only
through a cross-product term.** Expect this to be harder to make chaotic than the Cn/S4 families. If
the simple cubic ansatz refuses to go chaotic after a real search, **widen**: add a second invariant's
gradient, add higher-degree equivariants (e.g. `r²·grad(I_G)`), or couple a second radial variable.
Budget for this; it is the program's main technical risk.

### The construction recipe (general)
Dissipative ansatz, with `r² = x²+y²+z²`:
```
dX/dt = (ρ − d·r²) X  +  a·grad(I_G)(X)  +  b·( X × grad(I_G)(X) )
```
`I_G` = the group's lowest non-trivial invariant; `grad(I_G)` and `X × grad(I_G)` are both
G-equivariant (gradient-of-invariant, and cross-of-equivariants under SO(3)). 4 dials (ρ, d, a, b);
widen as needed.

### Tetrahedral T (order 12) — DERIVED AND VERIFIED FOR YOU
Generators (test against both): `σ3:(x,y,z)→(y,z,x)` (3-fold about (1,1,1)),
`σ2:(x,y,z)→(x,−y,−z)` (2-fold about x). Lowest invariant `I_T = xyz` (degree 3; invariant under T,
**broken by O** — this is the T-signature). `grad(xyz) = (yz, zx, xy)`.
`X × grad(xyz) = (x(y²−z²), y(z²−x²), z(x²−y²))`. The family:
```
ẋ = (ρ − d·r²) x + a·yz + b·x(y²−z²)
ẏ = (ρ − d·r²) y + a·zx + b·y(z²−x²)
ż = (ρ − d·r²) z + a·xy + b·z(x²−y²)
```
I verified both terms are T-equivariant under σ3 and σ2 by hand; still **write the test and watch it
pass** (and assert it is NOT O-equivariant: the `a·yz` term breaks O). Note the `(yz,zx,xy)` quadratic
is the classic gradient nonlinearity; the integrable core `ẋ=a·yz,…` is broken by the linear/radial
and rotational terms. Suggested search box: ρ∈[0.5,3], d∈[0.5,3], a∈[0.5,4], b∈[−3,3] — tune.

### Octahedral O (order 24) — DERIVED AND VERIFIED FOR YOU
Generators: `c4z:(x,y,z)→(−y,x,z)` (4-fold about z), `σ3:(x,y,z)→(y,z,x)` (O ⊃ T). Lowest invariant
beyond r² is `I_O = x⁴+y⁴+z⁴` (degree 4; `xyz` is NOT O-invariant). `grad = (x³,y³,z³)` (drop the 4
into `a`). `X × (x³,y³,z³) = (yz(z²−y²), zx(x²−z²), xy(y²−x²))`. The family:
```
ẋ = (ρ − d·r²) x + a·x³ + b·yz(z²−y²)
ẏ = (ρ − d·r²) y + a·y³ + b·zx(x²−z²)
ż = (ρ − d·r²) z + a·z³ + b·xy(y²−x²)
```
Verified O-equivariant under c4z and σ3 by hand; write the test anyway.

### Icosahedral I (order 60) — THE HARD ONE (do last; signature render)
Lowest invariant beyond r² is **degree 6**, `I_6`. You need (a) the 60 rotation matrices and (b) `I_6`:
- **Generators:** the Reiter "dodecahedron" preprint
  (`webbox.lafayette.edu/~reiterc/dodec/dodec02_pp.pdf`, already cited in RELATED_WORK.md) gives two
  generating rotation matrices `T1, T2` in terms of τ = (1+√5)/2 (golden ratio). Also see
  Golubitsky–Stewart–Schaeffer *Singularities and Groups in Bifurcation Theory* Vol II for the
  icosahedral invariant tables, and Field–Golubitsky *Symmetry in Chaos*.
- **`I_6` two ways:** (a) look up the closed-form degree-6 icosahedral invariant in the standard
  orientation (2-fold axes on coordinate axes); or (b) **construct it numerically, robustly:** generate
  the 60 matrices from T1,T2, pick a generic unit vector v, set `I_6(X) = Σ_{g∈I} ((g·v)·X)^6` — this
  is automatically I-invariant. Fit/expand it to a polynomial (or differentiate it directly for the
  flow). Verify by checking `I_6(g·X)=I_6(X)` for all g. Then `grad(I_6)` is the degree-5 equivariant
  for the ansatz. The I family will likely need the widening strategies above to go chaotic.
- Equivariance test: generate the full group, assert `F(g·X)=g·F(X)` for the two generators (suffices).

For each family: equivariance test → standalone QD search (copy `scratch/search_s4.py`, swap the
evaluator) → widen box until interior optimum → certify top-K finalists with the long integrator,
drop escapees → pick the most beautiful → full battery + symmetry residual → name (radiolarian
lineage) → render (new palette) → package → viewer → README/STORY/PAPER_PLAN. One atomic commit each.

---

## 4. Phase 3 — the chaos-vs-symmetry-group law (the robust finding)

This is the spine of the paper regardless of Phase 2's outcome. The story so far is already
**non-monotonic**: λ1_max so far C2 0.296, C3 0.233, C4 0.525, S4 0.330. Make it systematic.

**Protocol:**
- Groups: Cn for n = 2..8 (use `atlas/family.py`, already n-parameterized), S4 (`atlas/family_s4.py`),
  T/O/I (Phase 2 families). Optionally Dn (dihedral) — needs an extra in-plane antisymmetry; the Cn
  family is *not* Dn because the `iω` rotation breaks `w→conj(w)` (only Dn when ω=0). Construct a Dn
  family only if you want the data point; it's optional in the master plan.
- **Matched budget:** identical evaluator settings and the **same number of candidate evaluations**
  per group (e.g. 20k). Record per group: order |G|, n_eval, **max λ1 that survives long-integrator
  certification** (NOT the optimistic batch estimate — re-certify the top-K), the D_KY at that max,
  and the **chaotic fraction** (bounded & λ1 > threshold).
- **Fairness caveat to state honestly:** families have different param counts (Cn 7, S4 8, T 5, O/I 4),
  so chaotic-fraction is only loosely comparable across families; the robust comparable is
  **max certified λ1 under equal evaluation budget**. Say so.
- **Null hypothesis to cite** (from RELATED_WORK.md): in the Letellier–Gilmore cover/image
  construction the cover's LLE is *independent of the covering group* (cover is locally identical to
  its image). Our families are independently searched, not covers of a shared image, so any
  group-dependence we measure is a real property of the equivariant design space, not an artifact.
- **Outputs:** `results/symmetry_law.csv` (columns: group, order, n_eval, max_lambda1_certified,
  D_KY_at_max, chaotic_fraction, notes), one figure `gallery/symmetry_law.png` (λ1_max vs group,
  ordered by group order), and an analysis section added to README and to the paper's revised §4.7.
  Use main's validated novelty noise-floor for any novelty statements.

---

## 5. Deliverable bookkeeping & acceptance criteria

- Final program deliverable: the paper's **§4 = bestiary + symmetry law** (`docs/PAPER_PLAN.md`
  already reframed; `paper/main.tex` lives on `main` — update it after merging). STORY.md gains a
  chapter per specimen (Mobula chapter is the template). README gains a section per specimen.
- Every certification claim → a script + a `results/*.json`. Every render → viewed (PNG). Every new
  family → equivariance unit-tested numerically before any search. Novelty measured + hedged. Names
  marine + verified unclaimed. One logical change = one atomic commit, pushed.
- `pytest tests/` green at all times.
- Report at each phase boundary: Phase 2 → certified numbers per polyhedral specimen + the I render;
  Phase 3 → `symmetry_law.csv` + the λ1-vs-group figure + the analysis.

## 6. Key files / pointers
- Master plan: `docs/SYMMETRY_PROGRAM.md`. Phase-0 verdicts: `docs/RELATED_WORK.md`.
- Worked example (mirror this): `atlas/family_s4.py`, `mobula/`, `tests/test_s4_equivariance.py`,
  `tests/test_mobula_matches_family.py`, `scripts/verify_mobula.py`, `scripts/render_mobula.py`,
  `results/mobula_verification.json`, and the `index.html` SYSTEMS `mobula` entry.
- Search infra: `atlas/family.py` (Cn), `atlas/fingerprint.py` (D2+PCA fingerprints, explicit bins),
  `atlas/map_elites.py`, `results/catalog_fingerprints.json` (135 dysts systems).
- Paper (on `main`, in worktree `…/attractor-atlas-paper`): `paper/main.tex`, plus
  `results/convergence.json`, fingerprint-validation results, Shilnikov results from the paper
  workstream — use these.
