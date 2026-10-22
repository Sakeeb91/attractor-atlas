# The Jellyfish in the Equations

*How a hunt for beauty in a seven-dimensional parameter space turned up a new strange attractor, organized by a single unstable point, that looks like a moon jellyfish and answers to three numbers.*

![The Aurelia attractor](gallery/aurelia_hero.png)

There is exactly one place in the Aurelia system where nothing happens.

It sits on the vertical axis, at the point (0, 0, 1.523), and if you could balance there, the equations would hold you motionless forever. You cannot balance there. The point repels almost everything that approaches it, flinging nearby trajectories outward in a widening spiral, over the rim of a bell-shaped shell, down along translucent veils, and back up through a narrow central column to spiral out all over again. No trajectory ever repeats. No trajectory ever leaves. The result, traced over millions of orbits, is the luminous object above: a strange attractor that did not exist in the literature before June 10, 2026, and that looks, from the side, remarkably like a moon jellyfish.

This is the story of how it was found, why it has the shape it has, and what it took to certify that the thing in the picture is genuine chaos rather than an elaborate accident.

## An old kind of object

A dynamical system is just a rule for how a point moves: given where you are, the rule tells you where you go next. Weather, planetary orbits, the voltage in a neuron, all can be written this way. For three centuries after Newton, mathematicians expected such rules to produce tame behavior in the long run. Trajectories should settle down: to a standstill, or to a repeating loop called a limit cycle.

In 1963 the meteorologist Edward Lorenz found a third possibility hiding in a toy model of convection. His trajectories settled onto something, but the something was neither a point nor a loop. It was an infinitely layered butterfly-shaped set on which motion remained forever erratic, forever bounded, and exquisitely sensitive: two trajectories starting a hair apart would diverge exponentially until they bore no resemblance to each other. The set attracted everything around it, yet motion on it never repeated. Such objects came to be called strange attractors, and they are the geometric face of chaos.

Since Lorenz, a bestiary has accumulated. Rössler's single spiral fold. Thomas's cyclically symmetric weave. The Aizawa attractor's rotating column. Sprott's alphabetized menagerie of minimal chaotic flows. Each is a specific system of ordinary differential equations, usually three of them, with specific constants, that happens to fold space onto a fractal set. Finding a new one is not hard if you do not care what it looks like. Finding one worth looking at is a different problem.

## Designing with symmetry

The hunt began with an aesthetic constraint, not a mathematical one: the attractor had to be beautiful. That constraint immediately suggested symmetry.

In the early 1990s, the mathematicians Michael Field and Martin Golubitsky showed how to build chaotic maps whose strange attractors carry the symmetries of polygons and wallpaper patterns. Their book *Symmetry in Chaos* is full of images that look like stained-glass icons, and the trick behind them is a single idea called equivariance. A system is equivariant under a symmetry if the symmetry maps solutions to solutions: rotate a trajectory by 120 degrees, and you get another valid trajectory. An attractor of such a system either has that symmetry itself or comes in a symmetric family of copies. Field and Golubitsky's icons, however, are two-dimensional discrete maps: snapshots, not flows. The motion jumps from point to point.

The question that started this project was simple: what happens if you demand the same kind of symmetry from a continuous three-dimensional flow, the kind of system Lorenz wrote down?

The cleanest way to force three-fold rotational symmetry is to work in a complex variable. Write the horizontal coordinates as a single complex number, w = x + iy, so that rotating the plane by 120 degrees means multiplying w by e^(2πi/3). Most polynomial terms break this symmetry. But the conjugate-squared term w̄² respects it perfectly: rotating w by 120 degrees rotates w̄² by exactly the same angle, because the conjugate reverses the rotation and squaring doubles it, and minus two-thirds of a turn is the same as plus one-third. It is the unique quadratic term with this property — a fact equivariant bifurcation theory catalogued long ago, where the conjugate power w̄^(q−1) is the standard resonant coupling for q-fold symmetry. The design choice here is to let that one term carry all the symmetry-breaking in a flow tuned for chaos.

So the candidate family took the form

```
dw/dt = [α(z − β) + iω] w + γ w̄²
dz/dt = μ + νz − z³ − λ|w|²
```

Read it as a machine with three interlocking parts. The term iωw spins the horizontal plane like a turntable. The term α(z − β)w makes that spin unstable or stable depending on the altitude z: above the level z = β trajectories spiral outward, below it they are drawn back in. And the vertical equation is a thermostat: the cubic −z³ caps the altitude while the −λ|w|² term pushes z down whenever trajectories wander far from the axis. Energy pumps in along the axis, spills outward, gets reined back. The γw̄² term, the symmetry carrier, warps the circular spiral into a three-lobed one.

Written back in real coordinates, the system is

```
dx/dt = α(z − β)x − ωy + γ(x² − y²)
dy/dt = ωx + α(z − β)y − 2γxy
dz/dt = μ + νz − z³ − λ(x² + y²)
```

Seven constants. Seven dials, any setting of which is a different universe with different long-run behavior. Almost all of those universes are boring.

## Six hundred universes

There was a first attempt, and it failed in an instructive way. Before the equivariant family above, the search began with a different design: a cyclic system blending sine nonlinearities with bilinear cross-terms, a cousin of the Thomas attractor. It produced chaos readily enough, with largest Lyapunov exponents around 0.02 to 0.04, but the attractors were thin, wiry bands, like a few loops of bent wire. Mathematically chaotic; visually dead. The lesson: weak chaos looks weak. An attractor only develops layered, veil-like structure when the stretching is vigorous, and the stretching rate is exactly what the largest Lyapunov exponent measures. Beauty, it turns out, correlates with a number.

So the second campaign made that number the objective. A Monte-Carlo search drew random settings of the seven dials, six hundred of them, and for each one integrated the equations and estimated the largest Lyapunov exponent: take two trajectories separated by one part in a hundred million, watch how fast the gap grows, renormalize before it saturates, average the growth rate. Settings whose trajectories flew off to infinity were discarded. Settings that converged to fixed points or limit cycles, exponent zero or negative, were discarded. The filter demanded a strongly positive exponent and an orbit that actually filled a respectable volume of space.

Out of six hundred random universes, exactly two passed.

The stronger of the two, with an exponent of 0.27, sat at the unlovely coordinates α = 1.444, β = 0.481, ω = 1.955, γ = 0.485, μ = 1.433, ν = 1.413, λ = 1.905. Rendered, it was immediately the keeper: a domed shell of nested veils with a bright column up the middle, and, seen from directly above, a triskelion, three petals chasing each other around the axis, the C₃ symmetry made visible.

| The crown | The bell |
|---|---|
| ![Crown view](gallery/aurelia_crown.png) | ![Bell view](gallery/aurelia_bell.png) |

Then came the small surprise that gave the system its final form. Those seven awkward constants were not seven independent quantities at all. Nudging each toward a rounder value, with the Lyapunov exponent as a guardrail, they collapsed: α, μ, and ν all wanted to be 1.4; γ and β both wanted to be 0.5; ω and λ both wanted to be 1.9. The chaos did not merely survive the rounding, it strengthened slightly: under matched long measurements, the original constants score an exponent of 0.222 and the rounded ones 0.229. The seven-dial machine was secretly a three-dial machine:

```
dx/dt = a(z − b)x − cy + b(x² − y²)
dy/dt = cx + a(z − b)y − 2bxy        with  a = 1.4,  b = 0.5,  c = 1.9
dz/dt = a(1 + z) − z³ − c(x² + y²)
```

Whether this collapse points at something structural or is a happy coincidence of this corner of parameter space is one of the open questions at the end of this story.

## The anatomy of the jellyfish

Why this shape? The answer lives at that single motionless point.

Set the right-hand sides to zero and solve. On the axis, where x = y = 0, the planar equations vanish identically and the vertical one becomes a(1 + z) = z³, a cubic with exactly one real root: z = 1.523. Off the axis, a numerical hunt from four hundred random starting points finds nothing. The Aurelia system has precisely one equilibrium; Lorenz, by comparison, has three, and Aizawa-type systems two. One-equilibrium chaotic flows have a small literature of their own — minimal jerk systems are the classic examples — so the spare count is not itself the novelty. What matters is the work this lone point does: everything in the picture is organized by it.

Linearize the flow there and the local character emerges. Two eigenvalues form the complex pair 1.432 ± 1.9i: an unstable spiral, unwinding outward at rate 1.432 while rotating at rate 1.9 (the rotation dial c, surfacing again). The third eigenvalue is −5.558: a powerfully attracting direction. The point is a saddle-focus, and the interplay of its outward spiral and inward crush is the engine of the whole attractor. Trajectories drift up the axis into the equilibrium's neighborhood, get caught in the unstable spiral, and are slung outward and over the rim. They fall along the outside of the bell, where the thermostat in the z-equation catches them, and the strongly contracting direction funnels them back to the base of the central column to ride up and be flung again. Stretch, fold, reinject, forever.

This recirculating geometry is the signature of what dynamicists call Shilnikov-type chaos, after the Russian mathematician who proved in 1965 that a saddle-focus whose trajectories nearly reconnect with the equilibrium forces the existence of infinitely many unstable periodic orbits, chaos in its purest combinatorial form. The Aurelia attractor wears that mechanism on its sleeve: the veils are the stretching, the bright column is the reinjection, and the gilded crown at the top, where the colors in the renders run to gold, is the slow neighborhood of the saddle-focus where trajectories linger before their next ejection.

## Certifying the strange

A picture can lie. A trajectory plotted for a finite time can look chaotic while secretly being a long transient on its way to a limit cycle, or a periodic orbit of enormous period. Certification means numbers, and the numbers form a checklist.

**The Lyapunov spectrum.** A three-dimensional flow has three Lyapunov exponents, measuring exponential stretching along three independent directions. The clean way to compute all three is to carry a small orthogonal frame along the trajectory, letting the linearized flow shear it while repeatedly re-orthogonalizing with QR decomposition; the logarithms of the stretching factors accumulate into the exponents. Over 600,000 steps, the Aurelia spectrum converges to

> λ₁ = +0.2327, λ₂ = +0.0013, λ₃ = −1.5119

This is the textbook signature (+, 0, −) of a strange attractor in three dimensions. The positive exponent is sensitivity: nearby trajectories separate by a factor of e every 4.3 time units. The zero exponent (numerically +0.001) is obligatory for any flow, reflecting the trivial freedom of sliding along the trajectory itself. The large negative exponent is the contraction that glues trajectories onto the attractor. And the sum, strongly negative, means volumes shrink: the system is dissipative, averaging a contraction rate of −1.29, so the attractor occupies zero volume. Rerunning the estimate from five random starting points gives exponents from 0.195 to 0.252, always positive; the chaos is not an artifact of one lucky initial condition.

**The fractal dimension.** Kaplan and Yorke conjectured that the exponents encode the attractor's dimension: count how many exponents you can sum before the total goes negative, then interpolate. For Aurelia, D = 2 + (λ₁ + λ₂)/|λ₃| = 2.155. More than a surface, less than a solid: the veils in the renders are not one sheet but infinitely many, packed into thickness zero, a structure the dimension 2.155 quantifies precisely.

**Boundedness.** Four million integration steps, no escape, every coordinate confined to a box roughly three units on a side. Combined with dissipativity, this rules out the embarrassing alternative that the pretty pictures depict a slow explosion.

All of it is reproducible from the repository with one command per claim, and recorded in [`results/aurelia_verification.json`](results/aurelia_verification.json).

## The edge of chaos

Chaos is not the system's only mood. Hold a and b fixed and sweep the rotation rate c from 1.2 to 2.6, computing the largest Lyapunov exponent at each of 281 settings, and the attractor's life story unfolds sideways:

![Bifurcation scan](gallery/bifurcation.png)

Below c ≈ 1.40, the exponent hugs zero: the system settles into limit cycles, closed loops of perfectly periodic motion, the fate Lorenz's predecessors expected of all differential equations. At c ≈ 1.40 chaos ignites, burning through a band to about 1.62 before collapsing back to periodicity. A second, broader chaotic continent runs from about 1.81 to 2.33, and the canonical c = 1.9 sits comfortably inside it. Within each band, narrow periodic windows slice through, places where the chaos momentarily crystallizes into order before dissolving again, the same alternation famously found in the logistic map. Past 2.33 the system retires into limit cycles for good, save for one last flicker of chaos near c = 2.40.

This interleaving matters for the certification, too: it shows the chaotic parameter set is robust, a band rather than a knife-edge, which is why small perturbations of the sliders in the [live viewer](https://sakeeb91.github.io/attractor-atlas/) bend the jellyfish without killing it.

## The name

*Aurelia aurita* is the moon jellyfish, the translucent, domed drifter found in every ocean. The resemblance in the side view needs no argument: the bell, the trailing veils, the luminous core like the jelly's central mass. The Latin root *aureus* means golden, which fits the crown where slow trajectories pool. The name was unclaimed; searches of the dynamical-systems literature in June 2026 found no Aurelia attractor, and no published three-dimensional flow combining discrete C₃ equivariance through the w̄² coupling with cubic vertical dynamics. Discrete rotational symmetry by itself, though, is not new in three-dimensional chaos, and honesty requires the family tree. The Lorenz attractor has n-fold symmetric covers — the proto-Lorenz constructions of Miranda and Stone — and Gilmore and Letellier built an entire theory of such covering systems; the multi-scroll circuit literature engineers rotationally symmetric attractors at will; Thomas's cyclically symmetric system carries a three-fold rotation about the space diagonal. What those constructions share is that the symmetry is imposed on top of known dynamics, by lifting an existing attractor through an angle-multiplying map (which leaves the vector field singular on the axis and the dynamics locally identical to the original) or by piecewise design. Aurelia's claim is narrower: a single smooth polynomial flow, equivariant from birth through the unique resonant coupling, a cover of nothing, with the whole attractor organized by one saddle-focus. Its closest relatives in spirit — Field and Golubitsky's symmetric icons (discrete, two-dimensional) and the Aizawa-Langford flow (continuously rotationally symmetric, two equilibria) — each lack what the other has. Aurelia lives in the gap between them: a continuous flow, a discrete symmetry, a single organizing saddle-focus. A sparsely populated corner, not virgin land; the full survey is in [docs/RELATED_WORK.md](docs/RELATED_WORK.md).

Honesty requires the caveats. The novelty claim is a search claim: not found, rather than proven absent, in a literature too large for any search to exhaust. And the chaos certification is numerical. The Lyapunov exponents and the dimension are measurements with error bars, not theorems. A rigorous proof of chaos for this system, of the kind Warwick Tucker achieved for the Lorenz attractor in 2002 using computer-assisted interval arithmetic, does not exist and would be a genuine contribution.

## The machine that hunts

Aurelia was found by hand, in the sense that a human looked at renders and said yes, that one. The obvious next question was whether the looking could be industrialized. Not the taste, which stays human, but the searching. Automating the search itself is an old idea — Sprott harvested thousands of chaotic maps by computer in the early 1990s, selecting by the same eye, and evolutionary algorithms have since bred chaotic flows by the tens of thousands as machine-learning training data. The part that had no precedent was searching for *difference*: measuring how far each candidate sits from everything already known, and keeping the far ones.

Two ingredients turned the one-off hunt into an engine. The first was a way to *measure* novelty instead of asserting it. Every attractor, known or candidate, gets reduced to a fingerprint: the histogram of pairwise distances between points on its normalized orbit (a classic shape signature, blind to rotation and scale), together with a few gross-anatomy numbers like how planar or how space-filling the cloud is. Against this, a reference library: the dysts catalog, William Gilpin's machine-readable collection of 135 chaotic systems from a century of literature, every Lorenz and Rössler and Sprott in one place. Novelty became a distance: how far is this candidate's fingerprint from the nearest thing anyone has ever seen?

The metric promptly audited the original discovery. Asked for Aurelia's nearest known neighbor among all 135 systems, it answered: the Aizawa attractor, at distance 0.185. That is exactly the relative the qualitative analysis had named, found independently by arithmetic. (A later validation run re-simulated everything and found the distance rock-stable at 0.184 — with a second system, the InteriorSquirmer, statistically tied with Aizawa. Out at the rim of shape space, ties between faraway neighbors are what isolation looks like.) And the distance is large: the catalog's own median nearest-neighbor spacing is 0.115, so by this ruler Aurelia sits farther from everything known than five-sixths of the known systems sit from each other.

The second ingredient was a search algorithm that does not converge. Optimizers find one best thing; what was wanted here was an atlas of *different* things. The tool is MAP-Elites, a quality-diversity method: behavior space is tiled into cells (by symmetry order, chaos strength, and body shape), and each cell remembers the most novel specimen that ever landed in it. Mutants of the archive's elites and fresh random draws compete in batches, a whole population integrated simultaneously. The family Aurelia came from generalizes naturally: the conjugate term w̄^(n-1) imposes n-fold symmetry for any n, so the engine searched rotation orders two through five at once. In its first run it evaluated 4,608 candidate universes in 127 seconds and kept 62.

## Naiad, the fountain

The first thing worth keeping came from the two-fold family, and it taught a lesson about novelty on the way in.

The archive's *most novel* two-fold specimen, distance 0.44 from everything known, turned out under full certification to be barely chaotic at all: its second Lyapunov exponent was negative and its dimension a hair above 2, a nearly-flat sheet wearing an exotic shape. Novelty search rewards weird, and weird includes degenerate. The certification battery is what disciplines the wishlist. The specimen that survived it, with a strongly positive exponent and the layered look of real chaos, scored about 0.3 — still more isolated in shape space than Aurelia itself. At that distance the *identity* of the nearest relative turns out to be rank-unstable: the SprottJerk system in the original run, a forced FitzHugh–Nagumo oscillator under re-simulation. What is robust is the isolation, not the name of the second-place finisher.

![The Naiad attractor](gallery/naiad_hero.png)

Rounded to clean constants, it became **Naiad**, named for the water nymphs of fountains and springs. At n = 2 the symmetry term is simply w̄, which in real coordinates is an anisotropic stretch: it feeds energy along one horizontal axis and drains it along the other, so the flow is symmetric only under a half-turn. The attractor is a fountain seen from the side, a wide turbulent bowl with a brilliant jet up the middle and a narrow stem below. Seen from above it is not a disk but a lens, the two-fold symmetry made visible as elongation.

Its papers are in order: Lyapunov spectrum (+0.296, 0, −2.696), dimension 2.110, bounded over four million steps, dissipative. And the anatomy repeats: exactly one equilibrium, a saddle-focus on the axis at (0, 0, 1.651), spiraling out at rate 1.166 while rotating nearly three times faster than Aurelia's engine, with a contracting direction of −6.18 to catch what falls. Off-axis equilibria are not merely absent; the algebra forbids them, since they would require the rotation rate to be smaller than the stretch, and Naiad's parameters put it well above. One motionless point again, holding a fountain open forever.

## Cassiopea, the star

The four-fold family was searched the same night, and it produced the strongest chaos of the whole project.

The winning candidate's exponent measured 0.490 raw. Rounding its seven constants to clean values did what rounding did for Aurelia, nudging the chaos *up*, to a certified λ₁ = 0.525: more than twice Aurelia's stretching rate, with the family's deepest fractal dimension, 2.218. At n = 4 the symmetry term is w̄³, the unique cubic that survives a quarter-turn of the plane. From above, the attractor is a four-armed pinwheel star with a white spiral core, like a star sapphire cut square. From the side, the familiar family silhouette: a layered bell, twin curtains, a bright jet.

| The star from above | The bell from the side |
|---|---|
| ![Cassiopea star](gallery/cassiopea_star.png) | ![Cassiopea bell](gallery/cassiopea_bell.png) |

The name is **Cassiopea**, and for once the biology cooperates completely. *Cassiopea* is the upside-down jellyfish, a real animal whose bell carries a four-leaf-clover marking: genuine four-fold symmetry, in a jellyfish, echoing a constellation. The lineage holds: moon jelly, water nymph, upside-down jelly. All three names were unclaimed.

Cassiopea also broke the family's pattern in an instructive way. The equilibrium hunt first reported three fixed points: the usual saddle-focus engine on the axis at (0, 0, 1.611), plus two saddles far below the attractor. That count was impossible. A C₄-equivariant flow cannot have exactly two off-axis equilibria; rotating an equilibrium by ninety degrees must give another one, so they come in quadruples. The symmetry was auditing the numerics: the search box had been too small for equilibria sitting at radius three, and widening it produced the missing pair. Five equilibria in all, a central engine and a distant square of saddles, the first member of the family with any structure off the axis. The square sits far from the attractor and the dynamics never visit it, but its existence marks the four-fold family as algebraically richer terrain.

By the fingerprint ruler, Cassiopea sits 0.17–0.22 from the nearest thing known across re-simulations (the Nose-Hoover oscillator, with Aizawa statistically tied), beyond the catalog's median spacing of 0.115.

So the family stood at three, ordered by symmetry: Naiad (C₂, λ₁ = 0.296), Aurelia (C₃, 0.233), Cassiopea (C₄, 0.525). Same skeleton in every case, a single saddle-focus pumping a bell-shaped recirculation; different symmetry group, different face. They are all in the [live viewer](https://sakeeb91.github.io/attractor-atlas/), one switch apart.

## Mobula, the ray

The three jellyfish share a hidden assumption. Each is symmetric under a *pure rotation* — a third of a turn, a half, a quarter, about the vertical axis. But rotations are only half of the symmetries of three-dimensional space. The other half are *improper*: a rotation followed by a reflection, the operations that turn a left hand into a right one. The next question almost asked itself. Is there a strange attractor whose symmetry is not a rotation at all, but a rotation welded to a flip — an object that is its own upside-down mirror image, and is symmetric under neither the turning nor the flipping alone?

Before building, a literature strike, because "nobody has done this" is a hypothesis until checked. The verdict, written up with citations: the artists had made *maps* with the symmetry of tetrahedra and dodecahedra (Reiter and colleagues, iterating polynomials pixel by pixel), and Letellier and Gilmore had built flows with *inversion* symmetry — the point reflection through the origin. But a continuous flow whose strange attractor carries a genuine rotoreflection — improper, yet not the inversion and not a plain mirror — appeared in no search. The land was clear.

The target symmetry is `σ: (w, z) → (i·w, −z)`: spin the plane a quarter turn *and* flip top to bottom, in one indivisible move. Do it twice and you get the half-turn `(w, z) → (−w, z)`; do it four times and you return to the start. The group is the order-four rotoreflection **S₄**, in the chemist's notation — the symmetry of a tetragonal disphenoid, and of a devil ray turning over as it leaps. Working out which terms a vector field may contain without breaking σ is a short exercise in bookkeeping, and it turns up two terms that the four-fold jellyfish Cassiopea is explicitly *forbidden* to have: a `z·w̄` coupling in the plane and an `Im(w²)` coupling in the height. They are the signature of the improper symmetry — each one, on its own, allowed under the rotation-plus-flip but destroyed by the quarter-turn alone. They are exactly what drags the symmetry down from C₄ to S₄.

This time the equivariance was not trusted to a derivation. Before any search ran, a test was written that throws two hundred random states and random parameters at the flow and checks that turning-and-flipping the input turns-and-flips the output, to one part in a trillion — and, just as important, that the *pure* quarter-turn and the *pure* flip both fail. The test was watched to pass before the first candidate was evaluated. Only then did the quality-diversity engine go looking, in a parameter box widened until the best specimens sat in its interior rather than against its walls.

What it found is a ray. From the side the attractor is a broad-winged silhouette of layered veils, two great pectoral curtains flaring from a narrow waist — and the top half and the bottom half are not reflections of each other but *quarter-turned* reflections, the rotoreflection made visible. From above it is a four-fold pinwheel wound around a vortex. The name is **Mobula**, the genus of the devil rays, which breach the surface and somersault clear of the water before falling back — animals whose whole signature motion is a spin married to a flip. The marine lineage holds, and the name was unclaimed.

![The Mobula attractor](gallery/mobula_hero.png)

| The pinwheel from above | The ray from the side |
|---|---|
| ![Mobula pinwheel](gallery/mobula_pinwheel.png) | ![Mobula wing](gallery/mobula_wing.png) |

Its papers are in order, and they carry a new kind of evidence. Lyapunov spectrum (+0.330, 0, −1.238), dimension 2.266, bounded over four million steps, dissipative; the largest exponent holds between 0.33 and 0.38 across five scattered starts. The symmetry shows up in the skeleton: because S₄ forbids a constant in the height equation, the on-axis fixed points are the perfectly symmetric triple `z(ν − z²) = 0` — the origin and a pair at `(0, 0, ±1.817)` that σ swaps and which therefore carry *identical* eigenvalues, an unstable spiral out at 3.98 against a contraction of −6.6. Four more saddles stand off the axis, arranged not in Cassiopea's square but in a single order-four orbit of the rotoreflection. And to prove the symmetry is real and not an illusion of the eye, the certification measures how far the attractor moves when σ is applied to it: 0.0136 of its own size, against 0.0455 when a *pure* flip is applied instead — the rotoreflection maps the cloud onto itself more than three times better than the flip does. By the fingerprint ruler Mobula sits 0.273 from the nearest known system, more novel than Aurelia herself.

So the family stands at four, and the table has gained a column. Three pure rotations and one rotoreflection: Naiad (C₂, λ₁ = 0.296), Aurelia (C₃, 0.233), Cassiopea (C₄, 0.525), Mobula (S₄, 0.330). The chaos does not climb with the order of the symmetry; it wanders. That wandering is the thread the [symmetry program](docs/SYMMETRY_PROGRAM.md) now means to follow — toward the polyhedral groups, and toward a law of how much chaos a symmetry will permit.

## The polyhedral wall

The plan, written down in confidence, was to keep climbing the symmetry ladder. After the cyclic groups and the rotoreflection would come the Platonic ones: a strange attractor with the symmetry of a tetrahedron, then an octahedron, then — the intended signature image of the whole program — an icosahedron, sixty rotations deep, the symmetry of a virus shell and a soccer ball and a radiolarian's glass skeleton. The construction recipe was even written out in advance: take the lowest invariant of each group, `xyz` for the tetrahedron, `x⁴+y⁴+z⁴` for the octahedron, a degree-six polynomial for the icosahedron; push the flow with the gradient of that invariant and spin it with a cross product. The families were built, and a test proved each one carried exactly the right symmetry, down to machine precision. Then the search for chaos was turned on, and nothing came back.

Not "nothing beautiful." Nothing at all. Across a hundred thousand parameter sets, every bounded orbit in the tetrahedral family wound down to a fixed point or settled into a closed loop. The same for the octahedron. The reason, once seen, is not bad luck; it is a wall, and it had been hiding in plain sight in the one-line caution of the program's own plan. The cyclic and rotoreflection groups act on space *reducibly* — they spin a plane and leave a height axis alone — and that split is exactly what let the earlier systems work: a contraction–rotation engine in the plane, a fold along the axis, the saddle-focus that organizes every one of the four jellyfish. The Platonic groups are different. They act *irreducibly*: there is no plane they preserve, no axis they leave alone, no proper subspace at all. And a theorem of Schur says that anything linear which commutes with such an action must be a plain multiple of the identity — the same number in every direction. The engine that made the jellyfish chaotic is, for the polyhedra, mathematically forbidden. Near the origin every direction simply expands at one rate, like a balloon. All of the chaos would have to be manufactured by the nonlinear terms alone.

It can't be — at least not at the degrees a clean construction allows. There is a quieter reason underneath the algebra. When the radial pull settles an orbit onto a sphere, the motion left over lives on a two-dimensional surface, and a famous result forbids two-dimensional flows from being chaotic at all: they have nowhere to tangle. To get chaos the orbit has to keep climbing in and out of its shell, and the symmetric, balloon-like radial force won't let it. Widening the recipe — higher-degree terms, a second rotational engine, forcing every fixed point unstable so that something has to keep moving — only ever produced a limit cycle: a symmetric closed loop, the polyhedron's symmetry made into a single orbiting thread, but not chaos.

There is one known way to make a symmetric flow chaotic that doesn't need the linear engine: replace the polynomials with sines, the "labyrinth" trick behind Thomas's cyclically symmetric attractor, where an orbit wanders a periodic maze of cells. Built into a tetrahedral flow it does, finally, produce a flicker of genuine chaos — but a weak one (a largest exponent near 0.04, a tenth of Mobula's), and a *localized* one, trapped in a single cell of the maze off to one side of the origin. A localized blob is not symmetric: the group's other eleven copies of it sit elsewhere, and the orbit visits only one. The deep trap is that a genuinely symmetric attractor wants to wrap around the origin, the one point every rotation fixes — and at the origin the sines straighten back into their own Taylor series, the polynomial flow that has no chaos to give. The symmetric place is the calm place.

This is, almost certainly, why nobody has published a chaotic *flow* with polyhedral symmetry, only chaotic *maps* (Field and Golubitsky's symmetric icons, Reiter's tetrahedra and dodecahedra): maps iterate in discrete jumps and never have to integrate a vector field past Schur's wall. The signature icosahedral image the program set out to render does not exist as a strange attractor — not, at any rate, anywhere a careful search could reach. What the polyhedral groups will host is the rung below chaos: equilibria and symmetric limit cycles, which *can* carry all sixty rotations. The flagship became a negative result, and a sharper question: the table of strange attractors ordered by symmetry has an edge, and the edge is the line between groups that leave a direction alone and groups that don't.

## The law of the cliff

That edge is exactly what the program's final figure measures. Give every group the same search budget — the same number of candidate flows, the same evaluator, the same fine-grained re-certification that catches the integration-step mirages — and record the strongest chaos that survives. Across the cyclic and rotoreflection groups the answer is a positive, wandering number (the four jellyfish and their siblings), confirming once more that chaos strength does not march in step with symmetry order. Across the three Platonic groups the answer falls to the floor. Plotted against the group, ordered by size, the picture is not a trend but a cliff: chaos is permitted right up to the boundary of reducibility, and forbidden past it. The wandering line and the cliff are the same finding seen twice — the amount of chaos a symmetry will permit depends less on how *much* symmetry there is than on *what kind* it is. The numbers, group by group, are in [`results/symmetry_law.csv`](results/symmetry_law.csv) and the figure in [`gallery/symmetry_law.png`](gallery/symmetry_law.png).

## Open questions

The day that produced this family left more questions than it answered.

- **Does the Shilnikov criterion hold?** In all three systems the saddle-focus eigenvalues satisfy the magnitude inequality in Shilnikov's theorem, but the theorem also requires a homoclinic orbit, a trajectory leaving the equilibrium and returning to it exactly. The hunt has now been run (`scripts/homoclinic_hunt.py`): shooting the unstable manifold and sweeping the rotation rate, Naiad's manifold returns to within 7×10⁻⁵ of its equilibrium, with a ladder of seven near-connection windows accumulating below the canonical parameter, and Cassiopea's to within 7×10⁻⁴. For those two, homoclinic connections almost certainly thread the family. Aurelia is the holdout: along its entire c-band the manifold keeps a standoff distance of 0.13, so its connection, if one exists, requires moving a second parameter. Proving any connection rigorously remains open.
- **Why do seven constants collapse to three?** Aurelia's rounding may be a coincidence of its basin in parameter space, or the three-constant subfamily may be distinguished, perhaps the symmetric system maximizes stretching under some constraint. The pattern repeated for Cassiopea, where rounding raised the exponent again. Twice is suggestive; nobody knows.
- **Why does the four-fold family stretch hardest?** Chaos strength is not monotonic in symmetry order (C₂ 0.296, C₃ 0.233, C₄ 0.525). Whether the cubic equivariant term is intrinsically better at folding, or Cassiopea just sits in a luckier corner of its parameter space, is unexplored.
- **The five-fold gap.** The engine searched n = 5 too, and its five-fold elites were all weakly chaotic. Whether a strong five-pointed sibling exists, a starfish to complete the tide pool, likely needs a longer, targeted run.
- **Are there symmetry-broken siblings?** Equivariant systems can host attractors that break the symmetry and come in conjugate families, merging into one symmetric attractor at a symmetry-restoring crisis. Whether these parameter spaces contain such regimes, three small jellyfish orbiting where one now swims, is unexplored.
- **Does Cassiopea's saddle square matter?** Its four distant off-axis saddles never touch the attractor, but they may organize the boundary of its basin of attraction. Nobody has looked.
- **The other hit.** The original Monte-Carlo search produced a second three-fold chaotic universe, at λ ≈ 0.19, never examined beyond a thumbnail. It is still there, in the search script's output, waiting.

The equations fit on three lines apiece. The objects they contain have infinitely many layers, fractional dimensions, and the silhouettes of creatures that have drifted through the oceans for five hundred million years. You can fly around all four [in your browser](https://sakeeb91.github.io/attractor-atlas/), one switch apart, or trace them yourself from a few lines of Python in this repository. Either way, somewhere in the middle of each one, a single motionless point is holding the whole thing together by pushing everything away.

---

*All quantities cited are computed by the scripts in this repository: the Lyapunov spectra and dimensions by `scripts/verify_chaos.py`, `scripts/verify_naiad.py`, `scripts/verify_cassiopea.py`, and `scripts/verify_mobula.py`; the parameter sweep by `scripts/bifurcation_scan.py`; the original discovery search by `scripts/search_parameters.py`; the quality-diversity engine by `scripts/run_atlas.py` over `atlas/` (and its S₄ family in `atlas/family_s4.py`); the numerical equivariance and cross-check proofs by `pytest tests/`; the literature strike by `docs/RELATED_WORK.md`; and the renders by the `scripts/render_*.py` family.*
