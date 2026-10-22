# Methodology: how the Attractor Atlas was built, searched, certified, and measured

This document is the complete methods record of the attractor-atlas program: every
construction, every formula, every search protocol, every certification procedure, and the
exact chain of evidence behind each conclusion — including the negative result (the
polyhedral irreducibility obstruction) and the artifact that nearly contaminated it. Every
number cited here is produced by a script in this repository and stored in `results/`; the
[reproducibility map](#13-reproducibility-map) at the end links each claim to its code.

Companion documents: [README.md](../README.md) (overview),
[STORY.md](../STORY.md) (narrative), [SYMMETRY_PROGRAM.md](SYMMETRY_PROGRAM.md) (the
research plan and outcomes), [RELATED_WORK.md](RELATED_WORK.md) (literature positioning,
two parts), and the paper draft `paper/main.tex`.

---

## Contents

1. [Scope and claims](#1-scope-and-claims)
2. [Mathematical background and notation](#2-mathematical-background-and-notation)
3. [Phase 0 — literature verification](#3-phase-0--literature-verification)
4. [Symmetry-first design of the candidate families](#4-symmetry-first-design-of-the-candidate-families)
5. [The equivariance-test-first protocol](#5-the-equivariance-test-first-protocol)
6. [Search methodology](#6-search-methodology)
7. [The certification battery](#7-the-certification-battery)
8. [The five certified specimens](#8-the-five-certified-specimens)
9. [The step-size artifact and the dt-robustness gate](#9-the-step-size-artifact-and-the-dt-robustness-gate)
10. [The polyhedral irreducibility obstruction](#10-the-polyhedral-irreducibility-obstruction)
11. [The chaos-vs-symmetry-group law](#11-the-chaos-vs-symmetry-group-law)
12. [Numerical environment and known pitfalls](#12-numerical-environment-and-known-pitfalls)
13. [Reproducibility map](#13-reproducibility-map)
14. [Limitations](#14-limitations)
15. [References](#15-references)

---

## 1. Scope and claims

The program asks one question: **which finite symmetry groups of three-dimensional space
host minimal strange attractors, and how does the strength of the chaos depend on the
group?** It answers with:

1. **Five certified specimens**: explicit polynomial flows on $\mathbb{R}^3$, each
   equivariant under a chosen point group, each certified chaotic by the battery of §7:
   **Naiad** ($C_2$), **Aurelia** ($C_3$), **Cassiopea** ($C_4$), **Mobula** ($S_4$
   rotoreflection — to our knowledge the first *designed* attractor whose symmetry group is
   an improper rotation but not the corresponding proper one), and **Trinacria** (the
   chiral dihedral group $D_3$, added later and certified under the extended battery of
   §7.9; its basin gate was relaxed and the escape fraction is disclosed, §8).
2. **A negative result that is itself the finding**: the polyhedral rotation groups
   $T$ (order 12), $O$ (24), $I$ (60) act irreducibly on $\mathbb{R}^3$; in an exhaustive,
   step-size-robust search of their low-degree equivariant flows, **no strange attractor
   exists** — bounded orbits relax to equilibria or limit cycles (§10). This is strong
   numerical evidence, not a theorem.
3. **A quantitative law**: under matched search budgets, the maximum certified largest
   Lyapunov exponent is 0.18 to 0.69 across the reducible-action groups ($C_2..C_8$, $S_4$,
   $D_3$),
   non-monotonically in group order, and ≈0 for the irreducible polyhedral groups — the
   *irreducibility cliff* (§11).

Hedges are part of the claims: novelty is measured against a 135-system catalog (not the
literature); the obstruction is numerical evidence over stated ansätze (not a proof); the
chaotic-fraction column of the law is only loosely comparable across families (§11.3).

---

## 2. Mathematical background and notation

**Flows.** A system $\dot{X} = F(X)$, $X \in \mathbb{R}^3$, $F$ polynomial (or, in §10.5,
trigonometric). A **strange attractor** is a bounded attracting invariant set with
sensitive dependence — operationally here: largest Lyapunov exponent $\lambda_1 > 0$,
orbits bounded for millions of steps, fractal Kaplan–Yorke dimension $2 < D_{KY} < 3$,
negative mean divergence.

**Lyapunov exponents.** For the tangent flow $\dot{v} = J(X(t))\,v$ with
$J = \partial F/\partial X$, the exponents are
$\lambda_i = \lim_{t\to\infty} \tfrac{1}{t}\ln \sigma_i(t)$ where $\sigma_i$ are the
singular values of the fundamental matrix. For a flow on an attractor that is not a fixed
point, one exponent is 0 (the flow direction). Numerical methods: §6.2 (two-orbit
estimate) and §7.1 (full QR spectrum).

**Kaplan–Yorke dimension.** With exponents sorted $\lambda_1 \ge \lambda_2 \ge \lambda_3$
and $j$ the largest index with $\sum_{i\le j}\lambda_i \ge 0$:

$$D_{KY} \;=\; j + \frac{\sum_{i\le j}\lambda_i}{|\lambda_{j+1}|}.$$

**Equivariance.** A flow is equivariant under a group $G$ acting linearly on
$\mathbb{R}^3$ iff

$$F(g\,X) \;=\; g\,F(X) \qquad \text{for all } g \in G .$$

Then $g$ maps solutions to solutions, and an attractor is either $G$-symmetric as a set or
comes in a conjugate family $\{g\,A\}$.

**The groups used** (Schoenflies notation, all as subgroups of $O(3)$):

| group | order | action on $\mathbb{R}^3$ | generators used |
|---|---|---|---|
| $C_n$ | $n$ | rotation by $2\pi/n$ about $z$ — **reducible** (plane ⊕ axis) | $w \mapsto e^{2\pi i/n} w$, $z \mapsto z$ |
| $S_4$ | 4 | rotoreflection: quarter-turn **and** flip — reducible, improper | $\sigma: (w,z) \mapsto (i\,w,\,-z)$ |
| $T$ | 12 | tetrahedral rotations ($\cong A_4$) — **irreducible** | $\sigma_3:(x,y,z)\!\to\!(y,z,x)$, $\sigma_2:(x,y,z)\!\to\!(x,-y,-z)$ |
| $O$ | 24 | octahedral rotations ($\cong S_4^{\text{perm}}$) — irreducible | $c_{4z}:(x,y,z)\!\to\!(-y,x,z)$, $\sigma_3$ |
| $I$ | 60 | icosahedral rotations ($\cong A_5$) — irreducible | $g_5 = R\big((0,1,\varphi),\,2\pi/5\big)$, $g_3=\sigma_3$; $\varphi = \tfrac{1+\sqrt5}{2}$ |

Throughout, $w = x + iy$ identifies the $(x,y)$-plane with $\mathbb{C}$, and
$r^2 = x^2+y^2+z^2$.

**Why reducibility matters (the central structural fact).** For $C_n$ and $S_4$ the action
splits $\mathbb{R}^3 = \mathbb{C}_w \oplus \mathbb{R}_z$, so an equivariant *linear* part
can be $\big(\alpha + i\omega\big) w$ on the plane and an independent rate on the axis: a
**saddle-focus** (complex contraction/expansion + rotation in the plane, opposite-sign
rate on the axis) is available at the linear level. For $T, O, I$ the 3-dimensional
representation is absolutely irreducible, so by Schur's lemma the only linear map
commuting with the group is $\rho\,\mathbb{1}$ — every direction expands or contracts at
the *same* rate, no linear rotation, **no saddle-focus engine**. All of §10 unfolds from
this sentence.

---

## 3. Phase 0 — literature verification

Before construction, three novelty hypotheses were checked against the literature
(targeted web search, primary-source PDFs; full citations and search log in
[RELATED_WORK.md](RELATED_WORK.md) Part B):

1. **Polyhedral chaotic flows** — no prior exhibit found. All apparent prior art is
   (a) iterated *maps* (Brisson et al. 1996; Reiter 1997 tetrahedron, dodecahedron
   preprint; Field–Golubitsky *Symmetry in Chaos*), (b) non-chaotic flows, (c) chaos in
   dimension 4, or (d) abstract admissibility theory (Melbourne–Dellnitz–Golubitsky) with
   no constructed example.
2. **Rotoreflection-symmetric designed attractors** — none found; $S_4$ specimen proceeds
   as a novelty claim.
3. **A systematic chaos-vs-symmetry-group comparison** — none found; the law is open
   territory. Relevant null: in the Letellier–Gilmore cover/image construction the
   cover's Lyapunov exponent is *independent* of the covering group (§11.4).

A methods-paper positioning pass (Part A of RELATED_WORK.md) additionally calibrated which
claims about the *search method* survive (quality-diversity + catalog-grounded novelty:
survives, with ancestors Pan & Das 2015 and the 2025 Panda corpus to cite).

---

## 4. Symmetry-first design of the candidate families

The design principle inverts the usual order: **choose the group first, then write down
the most general low-degree equivariant vector field with dissipation dials, then search
the family numerically.** Equivariant terms are obtained two ways: (i) complex resonance
conditions on the reducible actions; (ii) invariant theory (gradients of invariants and
cross products) on the irreducible ones.

### 4.1 The $C_n$ family (Naiad, Aurelia, Cassiopea live here)

Under $w \mapsto \varepsilon w$ with $\varepsilon = e^{2\pi i/n}$, a monomial
$w^p\bar{w}^q$ transforms as $\varepsilon^{\,p-q}$; equivariance of $\dot w$ requires
$p - q \equiv 1 \pmod n$. The lowest-order term that is *not* $|w|^2$-radial is
$\bar{w}^{\,n-1}$ (since $-(n-1)\equiv 1$): the classical $\mathbb{Z}_n$-resonant coupling
of equivariant bifurcation theory. The search family (`atlas/family.py`, 7 parameters):

$$
\begin{aligned}
\dot w &= \big[\alpha\,(z-\beta) + i\,\omega\big]\,w \;+\; \gamma\,\bar{w}^{\,n-1},\\
\dot z &= \mu + \nu z - z^3 - \lambda\,|w|^2 .
\end{aligned}
$$

The linear core $[\alpha(z-\beta) + i\omega]w$ is the saddle-focus engine: the plane
rotates at rate $\omega$ and expands/contracts at the $z$-modulated rate $\alpha(z-\beta)$,
while the axis equation feeds back through $-\lambda|w|^2$. The canonical specimens are
particular members:

- **Aurelia** ($n=3$): $\dot w = (a(z-b) + ic)w + b\,\bar w^2$,
  $\dot z = a(1+z) - z^3 - c(x^2{+}y^2)$, with $(a,b,c) = (1.4,\, 0.5,\, 1.9)$. In real
  coordinates: $\dot x = a(z-b)x - cy + b(x^2-y^2)$, $\dot y = cx + a(z-b)y - 2bxy$.
- **Naiad** ($n=2$): coupling $g\,\bar w$ — an anisotropic stretch $(+g$ on $x$, $-g$ on
  $y)$; $(a,b,\omega_{\text{rot}},g,\mu,\nu,\lambda) = (1.2, 0.68, 3.7, 1.7, 1.2, 2.0, 2.1)$.
- **Cassiopea** ($n=4$): coupling $g\,\bar w^3$,
  $(a,b,\omega_{\text{rot}},g,\mu,\nu,\lambda) = (2.0, 0.85, 1.8, 0.9, 1.28, 1.8, 2.8)$.

### 4.2 The $S_4$ rotoreflection family (Mobula)

The generator is $\sigma: (w, z) \mapsto (i\,w,\, -z)$ — a quarter-turn *plus* a vertical
flip. $\sigma^2: (w,z) \mapsto (-w, z)$ is the half-turn, $\sigma^4 = \mathbb{1}$; the
group $\{\mathbb 1, \sigma, \sigma^2, \sigma^3\}$ is cyclic of order 4. Since
$\overline{i w} = -i\bar w$, equivariance of $\dot w = f(w,\bar w,z)$,
$\dot z = h(w,\bar w,z)$ demands

$$
f(i w,\, -i\bar w,\, -z) = i\,f(w,\bar w,z), \qquad
h(i w,\, -i\bar w,\, -z) = -\,h(w,\bar w,z).
$$

The term audit (each verified numerically, §5):

| term | in | under $\sigma$ | verdict |
|---|---|---|---|
| $(\alpha z^2 - \beta + i\omega)\,w$ | $f$ | $z^2$ even, $w \to iw$ | **allowed** |
| $\bar w^3$ | $f$ | $(-i)^3 = i$ | **allowed** |
| $z\,\bar w$ | $f$ | $(-z)(-i\bar w) = i\,(z\bar w)$ | **allowed — the $S_4$ signature** (forbidden under pure $C_4$) |
| $z\,w$ | $f$ | $(-z)(iw) = -i(zw)$ | forbidden (this is why the $C_n$ core $a(z-b)w$ breaks the flip) |
| $\nu z,\; z^3,\; z|w|^2$ | $h$ | odd | allowed |
| $\operatorname{Im}(w^2)$ | $h$ | $\operatorname{Im}(-w^2) = -\operatorname{Im}(w^2)$ | **allowed — second signature** (forbidden under $C_4$) |
| constant, $|w|^2$ | $h$ | even | forbidden |

The family (`atlas/family_s4.py`, 8 parameters):

$$
\begin{aligned}
\dot w &= (\alpha z^2 - \beta + i\omega)\,w + \gamma\,\bar w^3 + \delta\, z\,\bar w,\\
\dot z &= \nu z - z^3 - \lambda\, z\,|w|^2 + \epsilon \operatorname{Im}(w^2).
\end{aligned}
$$

Setting $\delta = \epsilon = 0$ restores pure $C_4$; with either present, $C_4$ is broken
— so the attractor's symmetry is pinned at exactly $S_4$. Note a structural consequence
used in certification: $S_4$ forbids a constant in $\dot z$, so the on-axis equilibria are
the symmetric triple $z(\nu - z^2) = 0$.

In real coordinates (with $\bar w^3 = (x^3-3xy^2) - i(3x^2y - y^3)$,
$\operatorname{Im}(w^2) = 2xy$):

$$
\begin{aligned}
\dot x &= (\alpha z^2 - \beta + \delta z)\,x - \omega y + \gamma(x^3 - 3xy^2),\\
\dot y &= (\alpha z^2 - \beta - \delta z)\,y + \omega x + \gamma(y^3 - 3x^2y),\\
\dot z &= \nu z - z^3 - \lambda z (x^2{+}y^2) + 2\epsilon\, x y .
\end{aligned}
$$

The signature $\delta z\bar w$ is a $z$-modulated anisotropic stretch whose sign reverses
with $z$ — the mechanism that ties the quarter-turn to the flip.

**Mobula canonical parameters:**
$(\alpha,\beta,\omega,\gamma,\delta,\nu,\lambda,\epsilon) =
(1.41039,\ 0.670108,\ 2.426585,\ 0.3,\ 0.553378,\ 3.3,\ 2.135518,\ 1.215326)$.

### 4.3 The polyhedral families: invariant-theory construction

For $G \le SO(3)$, two systematic sources of $G$-equivariant vector fields are

1. **gradients of $G$-invariant polynomials**: if $I(gX) = I(X)$ then
   $\nabla I(gX) = g\,\nabla I(X)$;
2. **cross products of equivariants**: if $u, v$ are equivariant and $g \in SO(3)$ then
   $g u \times g v = g (u \times v)$.

The general dissipative ansatz used for every polyhedral group:

$$
\dot X \;=\; (\rho - d\,r^2)\,X \;+\; a\,\nabla I_G(X) \;+\; b\,\big(X \times \nabla I_G(X)\big)
\;\;[\,+\ \text{further equivariants}\,],
$$

where $I_G$ is the group's lowest non-trivial invariant beyond $r^2$.

**Tetrahedral $T$** (`atlas/family_t.py`). Lowest invariant $I_T = xyz$ (degree 3; it is
*anti*-invariant under the octahedral $c_{4z}$, since $(-y)(x)(z) = -xyz$ — this is
exactly what separates $T$ from $O$). The complete degree-$\le 3$ $T$-equivariant field
has **five** independent terms:

$$
\dot X = (\rho - d r^2) X + a\,\underbrace{(yz,\, zx,\, xy)}_{\nabla(xyz)}
+ b\,\underbrace{\big(x(y^2{-}z^2),\, y(z^2{-}x^2),\, z(x^2{-}y^2)\big)}_{X \times \nabla(xyz)}
+ c\,\underbrace{(x^3, y^3, z^3)}_{\tfrac14\nabla(x^4+y^4+z^4)} .
$$

The $a$-term is the $T$-signature (breaks $O$, and being even under $X \to -X$ also breaks
$T_h$, so the flow is chiral and carries exactly the proper group $T$). The $b$, $c$,
radial terms are $O$-equivariant; setting $a = b = 0$ restores full $O(3)$ symmetry.

**Octahedral $O$** (`atlas/family_o.py`). Lowest invariant beyond $r^2$ is
$J = x^4 + y^4 + z^4$ ($xyz$ is not $O$-invariant). With $K = \tfrac14\nabla J = (x^3,y^3,z^3)$:

$$
\dot X = (\rho - d r^2) X + a\,(x^3,y^3,z^3)
+ b\,\big(yz(z^2{-}y^2),\; zx(x^2{-}z^2),\; xy(y^2{-}x^2)\big),
\qquad X \times K \text{ being the } b\text{-term}.
$$

The cross term is even under inversion, breaking $O_h$ down to the proper rotation group
$O$ (chirality).

**Icosahedral $I$** (`atlas/family_i.py`). The icosahedral group has **no invariant of
degree $< 6$** beyond $r^2$ (the classical fact behind fullerene/quasicrystal symmetry),
so the construction is numerical and exact by design:

1. Build the 60 rotation matrices as the closure of
   $g_5 = R\big((0,1,\varphi),\, 2\pi/5\big)$ (Rodrigues rotation about a vertex axis;
   $\varphi = \tfrac{1+\sqrt5}{2}$) and $g_3: (x,y,z) \mapsto (y,z,x)$. (Caution
   discovered en route: a 2-fold + 5-fold pair about *coordinate-aligned* axes closes to
   the order-10 dihedral $D_5$, not order 60; the $\{g_5, g_3\}$ pair is verified to give
   exactly 60 elements, all proper orthogonal.)
2. Fix a generic unit vector $v$ ($= (0.31, 0.52, 0.79)$ normalized) and define

$$
I_6(X) \;=\; \sum_{g \in I} \big( (g\,v) \cdot X \big)^6,
\qquad
\nabla I_6(X) \;=\; 6 \sum_{g\in I} \big( (g v)\cdot X \big)^5 \,(g v),
$$

which is $I$-invariant **by construction** (the sum runs over the whole group: replacing
$X \to hX$ permutes the orbit $\{gv\}$). Invariance is verified numerically to
$\sim 10^{-16}$, and anisotropy is asserted ($I_6 \ne c\,r^6$ — otherwise the family would
secretly be $SO(3)$-symmetric). The gradient is rescaled by a fixed constant so its mean
norm is 1 on the unit sphere (keeps the $(\rho,d,a,b)$ search box comparable across
T/O/I). The ansatz is the same 4-dial form with $\nabla I_6$ (degree 5) as the equivariant.

**An exact integrability identity used in §10.** For any invariant $I$:

$$
\frac{d}{dt}\,r^2 \Big|_{b\text{-term}} = 2\,X \cdot (X \times \nabla I) = 0,
\qquad
\frac{d}{dt}\,I \Big|_{b\text{-term}} = \nabla I \cdot (X \times \nabla I) = 0 .
$$

The cross ("rotational") term conserves **both** $r^2$ and $I$: alone, it moves along the
closed curves $\{r^2 = \text{const}\} \cap \{I = \text{const}\}$ — it is an integrable
engine and cannot produce chaos by itself.

### 4.4 The $D_n$ dihedral family (Trinacria)

The chiral dihedral group $D_n$ (order $2n$) is generated by the rotation
$\rho_n: (w, z) \mapsto (e^{2\pi i/n} w,\, z)$ and the half-turn about the $x$ axis
$\kappa: (w, z) \mapsto (\bar w,\, -z)$, i.e. $(x, y, z) \mapsto (x, -y, -z)$. Its $2n$
elements are $n$ rotations about $z$ and $n$ half-turns about horizontal axes; all have
determinant $+1$, so the group contains no mirror and no inversion. Equivariance of
$\dot w = f(w, \bar w, z)$, $\dot z = h(w, \bar w, z)$ demands the $C_n$ conditions plus

$$
f(\bar w,\, w,\, -z) = \overline{f(w, \bar w, z)}, \qquad
h(\bar w,\, w,\, -z) = -\,h(w, \bar w, z).
$$

Term audit under $\kappa$ (each verified numerically, §5, `tests/test_d_equivariance.py`):

| term | in | under $\kappa$ | verdict |
|---|---|---|---|
| $(\alpha z^2 - \beta)\,w$ | $f$ | real coefficient, even in $z$ | allowed |
| $i\omega_0\, w$ (uniform spin) | $f$ | $i\omega_0 \bar w \ne \overline{i\omega_0 w}$ | **forbidden**: no uniform spin |
| $i\omega\, z\, w$ | $f$ | $i\omega(-z)\bar w = \overline{i\omega z w}$ | **allowed**: spin odd in $z$ |
| $\gamma\, \bar w^{\,n-1}$ | $f$ | real coefficient | allowed |
| $i\delta\, z\, \bar w^{\,n-1}$ | $f$ | odd $z$ times $i$ | allowed |
| $\nu z,\; z^3,\; z|w|^2$ | $h$ | odd | allowed |
| constant, $|w|^2$ | $h$ | even | forbidden |
| $\epsilon \operatorname{Im}(w^n)$ | $h$ | $\operatorname{Im}(\bar w^n) = -\operatorname{Im}(w^n)$ | allowed |

The family (`atlas/family_d.py`, 8 parameters, same box conventions as §4.2):

$$
\begin{aligned}
\dot w &= (\alpha z^2 - \beta + i\omega z)\,w + \gamma\,\bar w^{\,n-1} + i\delta\, z\,\bar w^{\,n-1},\\
\dot z &= \nu z - z^3 - \lambda\, z\,|w|^2 + \epsilon \operatorname{Im}(w^n).
\end{aligned}
$$

Three coefficient sets enlarge the group, and the gate test shows each set is minimal
(keeping any one member nonzero leaves the larger group broken):
$\omega = \delta = \epsilon = 0$ restores the $z$-flip and the vertical mirrors ($D_{nh}$);
$\omega = \gamma = 0$ restores the rotoreflection $S_{2n}$ and the diagonal mirrors
($D_{nd}$, which for odd $n$ contains the inversion); $\gamma = \delta = \epsilon = 0$
restores rotation by every angle. The test also shows that $C_{2n}$, the $z$-flip, the
vertical and diagonal mirrors, $S_{2n}$ and the inversion each fail at 40/40 random states
($n = 2$ to 6), and that every copy of $O$ and $I$ containing $D_n$ is broken ($n = 2$ to 5).

For $n = 3$, with $c = \alpha z^2 - \beta$:

$$
\begin{aligned}
\dot x &= c\,x - \omega z y + \gamma(x^2 - y^2) + 2\delta z x y,\\
\dot y &= c\,y + \omega z x - 2\gamma x y + \delta z (x^2 - y^2),\\
\dot z &= \nu z - z^3 - \lambda z (x^2{+}y^2) + \epsilon (3x^2 y - y^3).
\end{aligned}
$$

The $\omega$, $\gamma$ and $\delta$ terms are trace-free, so
$\nabla\!\cdot\!F = 2(\alpha z^2 - \beta) + \nu - 3z^2 - \lambda(x^2{+}y^2)$, the same closed
form as Mobula.

**Equilibrium skeleton.** On the axis $\dot z = z(\nu - z^2)$, so the on-axis equilibria
are the origin and the $\kappa$-conjugate pair $(0, 0, \pm\sqrt\nu)$. At $z = 0$ the spin
$\omega z$ vanishes and the planar linear part at the origin is $-\beta\,\mathbb 1$: the
origin is a saddle with a *double real* eigenvalue $-\beta$ and $\nu$ along $z$, not a
saddle-focus. The Shilnikov narrative of the $C_n$ specimens (§8) therefore does not
transfer. The on-axis pair has identical spectra and opposite rotation sense
($J_{xy} - J_{yx} = -2\omega z$ changes sign with $z$). Equilibria in the plane $z = 0$
lie on the half-turn axes and form orbits of size 3; generic equilibria form orbits of
size 6. The origin is the only point fixed by all of $D_3$, so a $D_3$-invariant attractor
must be centred there; an off-centre cloud is symmetry-broken.

**Trinacria canonical parameters** (nominee #8 of the fine-dt screen,
`results/d3_selection.json`, 4-decimal probe values with no polishing):
$(\alpha,\beta,\omega,\gamma,\delta,\nu,\lambda,\epsilon) =
(1.4773,\ 1.4668,\ 0.7168,\ 0.8446,\ 1.0567,\ 2.3978,\ 1.2975,\ -1.5491)$.
$\beta$ sits near the box face 1.5; since no optimization was performed, the face is not a
selection effect on the certified numbers.

---

## 5. The equivariance-test-first protocol

No family is searched before its symmetry is *proven numerically* (test-driven
development applied to mathematics). Each `tests/test_*_equivariance.py` enforces, with
machine-precision tolerances (`1e-12` absolute, looser only where degree-5 polynomials
amplify roundoff):

1. **Group sanity** — generators are proper rotations; the generated group has exactly the
   right order (12 / 24 / 60), computed by closure under multiplication.
2. **Full-group equivariance** — $F(gX) = g\,F(X)$ for *every* element $g$ of the
   generated group (not just the generators), at dozens of random states and random
   parameter vectors drawn from the search box.
3. **Exclusion of larger groups** — the flow must *not* be equivariant under the next
   group up the lattice, tested as ≈100% violation across random states:
   $S_4$ family ✗ pure $C_4$, ✗ pure flip; $T$ family ✗ octahedral $c_{4z}$, ✗ inversion
   (excludes $T_h$); $O$ and $I$ families ✗ inversion, ✗ random $SO(3)$ rotations.
4. **Load-bearing tests** — zeroing precisely the signature coefficients must *restore*
   the larger symmetry ($\delta{=}\epsilon{=}0 \Rightarrow C_4$ for $S_4$;
   $a{=}b{=}0 \Rightarrow O(3)$ for the polyhedral families). This pins the symmetry
   boundary to specific terms and guards against silently dropping one.
5. **Cross-checks for the named specimen** — the real-coordinate package must equal the
   complex-form family evaluator at canonical parameters to $10^{-12}$
   (`tests/test_mobula_matches_family.py`); the hand-derived Jacobian must match central
   finite differences ($h = 10^{-6}$, tolerance $10^{-5}$); the closed-form divergence
   must equal $\operatorname{tr} J$ to $10^{-10}$.

The suite (34 tests at program end) must stay green through every commit.

---

## 6. Search methodology

### 6.1 Batched evaluation

Each family module exposes `evaluate_batch(params, ...)`: an entire batch of $B$ parameter
vectors is integrated *simultaneously* with vectorized classical RK4,

$$
X_{k+1} = X_k + \tfrac{\Delta t}{6}\,(k_1 + 2k_2 + 2k_3 + k_4),
$$

with per-family step/duration defaults (Cn/S4: $\Delta t = 0.01$, transient 8 000 steps,
$\lambda$-window 25 000; polyhedral, post-§9: $\Delta t = 0.005$, transient 16 000,
window 50 000). Diverging members are **culled** every 200 steps (escape radius 50–60,
non-finite check) and frozen at 0 so they stop generating overflow; the `alive` mask
records survivors. Initial conditions are drawn near (but not at) the origin with the
symmetry deliberately broken, so a symmetric transient does not bias toward the
group-fixed point.

### 6.2 The two-orbit Benettin estimate of $\lambda_1$

A shadow orbit is seeded at distance $d_0 = 10^{-8}$ and renormalized every step:

$$
\lambda_1 \;\approx\; \frac{1}{N\,\Delta t}\sum_{k=1}^{N}\ln\frac{d_k}{d_0},
\qquad
X^{\text{shadow}}_{k} \leftarrow X_k + (X^{\text{shadow}}_k - X_k)\,\frac{d_0}{d_k}.
$$

This is the engine of the search (cheap, vectorizes across the batch). It is also
**systematically optimistic at coarse $\Delta t$** — the artifact dissected in §9. It is
never the certified number; it only nominates candidates.

### 6.3 Quality-diversity search (MAP-Elites)

The original $C_n$ atlas run (`scripts/run_atlas.py`, orders $n = 2..5$) is a MAP-Elites
loop: per generation, propose a batch (Gaussian mutations of randomly chosen elites,
$\sigma = 8\%$ of each parameter's box span, mixed with fresh uniform draws); evaluate;
fingerprint survivors (§6.4); score **novelty** as distance to the dysts catalog *plus*
the current archive; offer each survivor to a behavioral archive keyed by rotation order
and coarse shape descriptors, keeping the best per cell. The S4/T searches used the same
propose–evaluate–rank loop as standalone scripts with a top-$K$ survivor pool instead of
a grid archive.

### 6.4 Shape fingerprints and the novelty metric

`atlas/fingerprint.py` builds a 35-dimensional descriptor of a decimated orbit cloud,
invariant to translation, scale, and (statistically) rotation:

- normalize the cloud to zero mean and unit RMS radius;
- **D2 shape distribution** (Osada et al. 2002): histogram of pairwise distances among
  400 random points, 32 explicit bins on $[0,3]$, Hellinger-embedded
  ($h \mapsto \sqrt{h/\sum h}$);
- **PCA eigenvalue ratios** $e_2/e_1$, $e_3/e_1$ (planarity / volumetricity);
- **fill ratio**: occupied fraction of a $12^3$ grid over the bounding box.

Distance between fingerprints: $\;\|\Delta\text{hist}\|_2 + 0.5\,\|\Delta\text{feat}\|_2$.
**Novelty** of a candidate = minimum distance to the 135-system dysts catalog
(`results/catalog_fingerprints.json`).

The metric was validated in a dedicated study (`scripts/validate_fingerprint.py`):
re-simulation noise floor **0.0549** (15k points; 0.0835 for the catalog's 1500-point
fingerprints); median catalog nearest-neighbour spacing **0.115**; rule of thumb —
distances below $\sim 2\times$ the noise floor are not interpretable, and NN identities go
rank-unstable when several systems sit within one floor of each other. All published
novelty numbers are quoted against these validated ranges.

### 6.5 The beauty gate

Aesthetics is a hard acceptance criterion, encoded then human-confirmed. The ranking
heuristic (T-search form):

$$
\text{beauty} = \exp\!\Big(-\frac{(\lambda_1 - \lambda^\ast)^2}{2s^2}\Big)
\cdot (0.3 + \text{iso}) \cdot (0.3 + e_3/e_1) \cdot \big(0.5 + \min(\text{nov}, 0.5)\big)
$$

with target $\lambda^\ast \approx 0.3{-}0.45$, width $s \approx 0.25{-}0.30$; *iso* =
min/max bounding-box span (isotropy — a genuinely polyhedral attractor must fill all
three axes equally); $e_3/e_1$ penalizes flat sheets. Every shortlisted candidate is then
**rendered and looked at** (multi-view montages); "mathematically valid but ugly" is
grounds for rejection. Target profile from the program conventions: $\lambda_1 > \sim 0.2$,
$D_{KY} \in [2.1, 2.3]$, layered veils.

### 6.6 Search-box discipline

If the running optimum pins against a box face, the box is widened until the optimum is
interior (applied in the S4 and T searches; in T it exposed that the apparent optima were
chasing corners — an early symptom of §9). Batch nominees are *always* re-certified with
the long integrator before any naming decision; in Phase 1, several S4 "alive" candidates
escaped under the fine integrator and were dropped.

---

## 7. The certification battery

Run by `scripts/verify_<name>.py`; writes `results/<name>_verification.json`. A specimen
must pass **all** of:

### 7.1 Long Lyapunov spectrum (tangent-space QR / Benettin)

The full frame $Q \in \mathbb{R}^{3\times3}$ is advected by the linearized flow with a
tangent RK4 (state and frame integrated together):

$$
\dot Q = J(X)\,Q, \qquad
Q_{k+1} R_{k+1} = \text{qr}\big(Q_k + \tfrac{\Delta t}{6}(j_1 + 2j_2 + 2j_3 + j_4)\big),
\qquad
\lambda_i = \frac{1}{N \Delta t} \sum_k \ln |(R_k)_{ii}| .
$$

Settings for certification: $N = 600\,000$ steps at $\Delta t = 0.005$ after a 20 000-step
transient (Mobula). Acceptance: $\lambda_1 > 0$ robustly, $\lambda_2 \approx 0$ (flow
direction, $|\lambda_2| \lesssim 10^{-2}$), $\lambda_3 < 0$, $\sum_i \lambda_i < 0$.
A separate convergence study (`scripts/convergence_study.py`,
`results/convergence.json`) established error bars vs window length for all specimens.

### 7.2 Kaplan–Yorke dimension

Formula of §2 applied to the certified spectrum; acceptance $2 < D_{KY} < 3$ (fractal,
genuinely three-dimensional flow attractor).

### 7.3 Initial-condition robustness

$\lambda_1$ re-estimated from 5 scattered random initial conditions (150 000-step
windows); all must be positive and mutually consistent.

### 7.4 Boundedness and dissipation

A 4 000 000-step orbit at $\Delta t = 0.004$ must stay finite and inside a fixed ball
(extent recorded); the time-averaged divergence
$\overline{\nabla \!\cdot\! F} = \overline{\operatorname{tr} J}$ (closed form per family,
e.g. Mobula: $2(\alpha z^2 - \beta) + \nu - 3z^2 - \lambda(x^2{+}y^2)$, the $\delta$- and
cubic terms cancelling pairwise) must be negative.

### 7.5 Equilibria and their classification

On-axis equilibria found analytically (e.g. $S_4$: the symmetric triple
$z(\nu - z^2) = 0$, i.e. $(0,0,0)$ and $(0,0,\pm\sqrt{\nu})$); off-axis equilibria by
multi-start root finding (600 seeds, `scipy.optimize.fsolve`, $\|F\| < 10^{-10}$
acceptance, deduplication at $10^{-6}$). Each equilibrium is classified by the
eigenvalues of $J$ (stable/unstable node/focus, saddle, saddle-focus). The symmetry must
be visible in this skeleton: group-conjugate equilibria must carry identical spectra
(e.g. Mobula's $\sigma$-exchanged pair at $(0,0,\pm 1.817)$, and its four off-axis saddles
forming a single order-4 $\sigma$-orbit).

### 7.6 The symmetry-residual statistic (with control)

Proof that the *attractor as a set* carries the group, not just the equations. For the
centered orbit cloud $\{q_i\}$ and a group generator $g$:

$$
\rho_g \;=\; \frac{1}{\text{extent}}
\sqrt{\;\frac{1}{n}\sum_{i=1}^{n}\ \min_j \big\| g\,q_i - q_j \big\|^2 \;}
$$

(nearest-neighbour RMS of the transformed cloud against the original, normalized by the
cloud's largest axis extent; $n = 2500$ query points against an 8 000-point base). The
residual for every group generator must be small **and** a *control* transformation that
is **not** a symmetry must score several times worse. For Mobula:
$\rho_\sigma = 0.0136$ versus pure-flip control $0.0455$ ($3.3\times$) — the cloud maps
onto itself under the rotoreflection far better than under the flip, so the symmetry is
genuinely $S_4$ and not an illusion of near-mirror-symmetry.

### 7.7 Novelty (measured and hedged)

Fingerprint distance to all 135 catalog systems (§6.4), reported with nearest and
runners-up, interpreted against the validated noise floor (0.0549) and median spacing
(0.115). Always hedged: the catalog is 135 systems, not the literature.

### 7.8 Naming

Only after the full battery: marine lineage, verified unclaimed by search
("`<name>` attractor").

### 7.9 Extensions introduced for the dihedral specimen

Trinacria (`scripts/verify_trinacria.py`, `results/trinacria_verification.json`) was
certified with a stricter battery. Each item below is a gate unless marked as reported.

- **dt-halving arm** (the §9 lesson made a standing gate): $\lambda_1$ at
  $\Delta t \in \{0.01, 0.005, 0.0025, 0.00125\}$ over the same physical window
  ($T = 2000$ after $T = 100$) from the same 24 initial conditions, compared as ensemble
  means with across-orbit standard errors. Gates: successive relative changes at most 10%,
  the two finest within 5%, and every coarser mean within 3 combined SE of the finest (no
  systematic drift as $\Delta t \to 0$). A $T = 500$ spot check at $\Delta t = 0.000625$ is
  reported. The headline spectrum is the finest-step ensemble mean, not a single orbit.
- **Spectrum versus divergence**: $|\sum_i \lambda_i - \overline{\nabla\!\cdot\!F}| \le 3\%$
  of $|\overline{\nabla\!\cdot\!F}|$, with the divergence averaged along an independent long
  orbit (an internal consistency check of the tangent integrator).
- **Equilibria as group orbits**: the multistart roots are completed by applying all six
  group matrices, then partitioned into $D_3$-orbits. Gates: orbit sizes in $\{1, 2, 3, 6\}$,
  spectra identical within an orbit to $10^{-6}$, every $z = 0$ equilibrium on a half-turn
  axis, and opposite rotation sense ($J_{xy} - J_{yx}$) for the on-axis pair.
- **Symmetry residual with strong controls** (extends §7.6): $\rho_g$ for all five
  non-identity elements against the self floor $\rho_0$ (first half of the orbit against the
  second, both directions; base 60 000 points, 5 000 queries). Gates: $\max_g \rho_g \le
  1.5\,\rho_0$ and $\le 0.015$; every control at least $3\times$ worse (z-flip, y-mirror,
  $C_6$, inversion, a 37° rotation, all 144 elements of the copies of $O$ and $I$ that
  contain this $D_3$, and rotations about $z$ by $120°/k$, $k = 2..6$); $C_6$ and y-mirror at
  least $5\times$.
- **Single attractor versus a $\kappa$-conjugate pair**: a small $\kappa$ residual on one orbit
  is necessary but not sufficient, since a pair of $C_3$ attractors exchanged by $\kappa$ would
  also produce one. Gates: random initial conditions in $[-2,2]^3$ (48, plus a 480 sample)
  land with $\rho_{\text{id}} \le 0.02$ against the reference cloud; orbits launched from
  $\kappa x_0$ and $\rho x_0$ fill the same cloud; on a 64-orbit ensemble, the centroid is
  within 3 SE of the origin and below 1% of the extent, $P(z > 0)$ is within 3 SE of 1/2, the
  pooled KS distance between $z$ and $-z$ is at most 0.02, and the mean angular velocities
  above and below the mid-plane are equal and opposite within 5% (and 3 SE).
- **Basin comparison**: the same 48-IC protocol is applied to every specimen
  (`scripts/basin_comparison.py`, `results/basin_comparison.json`).

Two gates were changed after their first run, and both original verdicts are kept in the
JSON: the $z$-rotation scan control (the first version used the scan minimum at 5°, which is
close to the identity for any cloud) and the strict basin gate (§8). The image-launch
fingerprint distance is reported, not gated, because the fingerprint ignores rotations.

---

## 8. The five certified specimens

All numbers from `results/*_verification.json` (validated novelty ranges from
`results/fingerprint_validation.json` in brackets where they differ from the discovery-run
value):

| | **Naiad** | **Aurelia** | **Cassiopea** | **Mobula** | **Trinacria** |
|---|---|---|---|---|---|
| group | $C_2$ | $C_3$ | $C_4$ | $S_4$ (improper) | $D_3$ (chiral dihedral) |
| spectrum $\lambda_{1,2,3}$ | $+0.296,\ 0.001,\ -2.696$ | $+0.233,\ 0.001,\ -1.512$ | $+0.525,\ 0.000,\ -2.412$ | $+0.329,\ 0.000,\ -1.238$ | $+0.328,\ -0.001,\ -1.492$ (24-orbit mean, $\Delta t = 0.00125$; quote $\lambda_1 = 0.33 \pm 0.01$) |
| $D_{KY}$ | 2.110 | 2.155 | 2.218 | 2.266 | 2.219 |
| mean divergence | $-2.386$ | $-1.290$ | $-1.868$ | $-0.904$ | $-1.164$ ($\sum\lambda_i = -1.165$) |
| nearest catalog | SprottJerk 0.34 [0.288, NN unstable] | Aizawa 0.185 [0.184, stable] | NoseHoover 0.224 [0.173, stable] | NuclearQuadrupole 0.273 | Coullet 0.212 [NN unstable: JerkCircuit 0.241] |
| $\lambda_1$ across 5 ICs | positive | positive | positive | 0.334–0.384 | 0.305–0.358 |
| dt-halving gate | not gated | not gated | not gated | not gated | 0.332, 0.332, 0.332, 0.328 ($\Delta t$ = 0.01 … 0.00125) |
| symmetry residual | — (visual + equilibria) | — | — | $\rho_\sigma{=}0.0136$ vs control 0.0455 | max over $D_3$ 0.0046 vs floor 0.0045; controls $\ge 15\times$ |
| basin, 48 ICs in $[-2,2]^3$ | 48/48 | 48/48 | 43/48 | 41/48 | 24/48 bounded, all bounded orbits land |

Trinacria is the only specimen whose central equilibrium is not a saddle-focus (§4.4), and
the only one certified under a relaxed gate: its pre-registered basin gate (every one of 48
random initial conditions lands on the attractor) fails because half of the box escapes to
infinity (56 ± 2% of a 480-IC sample stay bounded). It was accepted under the gate "every
bounded orbit lands on the symmetric set" (270/270), by an explicit decision recorded with
its reasons in `results/d3_selection.json` (`user_decision`).

Additional structural evidence: in all three $C_n$ specimens the central equilibrium is a
saddle-focus satisfying Shilnikov's eigenvalue (magnitude) condition; the homoclinic hunt
(`scripts/homoclinic_hunt.py`) shoots the 1-D unstable manifold and sweeps the rotation
rate, finding returns to within $7\times10^{-5}$ of the equilibrium for Naiad (with a
ladder of seven near-connection windows), $7\times10^{-4}$ for Cassiopea, and a standoff
of 0.13 for Aurelia — strong evidence (not proof) that homoclinic connections thread the
first two.

---

## 9. The step-size artifact and the dt-robustness gate

**The single most important methodological lesson of the program.** The two-orbit
Benettin estimator at coarse $\Delta t$ does not merely overestimate $\lambda_1$; around
*stable equilibria of strongly radial flows* it **fabricates chaos outright**.

**The discovery.** In the tetrahedral search, batch candidates reported
$\lambda_1 \approx +0.4$ to $+2.3$ at $\Delta t = 0.01$ with bounded, isotropic-looking
orbit samples. Re-certifying the leading candidate
$(\rho, d, a, b) = (7.17,\ 0.491,\ 7.33,\ -2.80)$ with the long fine integrator:

| $\Delta t$ | behaviour |
|---|---|
| 0.01 | apparent attractor, span $\approx [3.79, 3.79, 3.79]$, $\lambda_1 \approx +0.47$ |
| 0.005 | **collapses to the fixed point** $(5.814, 5.814, 5.814)$, span $= 0$ |
| QR spectrum at 0.005 | $\lambda = (-56.9,\ -86.2,\ -86.2)$ — a deeply stable equilibrium |

Every "chaotic" polyhedral candidate behaved the same way (e.g. $-15$ to $-24$ at
$\Delta t = 0.0025$ in the stability-targeted scan). Mechanism: at coarse step, RK4
overshoots around a stiff stable spiral, sustaining a spurious bounded "orbit" whose
shadow-separation grows; refining the step reveals the true sink.

**The gate (now standard for every search).** A candidate counts as chaotic only if:

1. $\lambda_1 > $ threshold at $\Delta t = 0.005$ (time-window-matched transient/window),
   **and**
2. $\lambda_1 > $ threshold at $\Delta t = 0.0025$ with double the steps (same physical
   time), **and**
3. the two estimates agree to within $\sim 45{-}50\%$ relative difference, **and**
4. for finalists: the *continued-orbit* long integrator at $\Delta t \le 0.0025$ (and spot
   checks at $0.00125$, $0.000625$) confirms a converged positive exponent. A $\lambda_1$
   that *changes systematically* as $\Delta t \to 0$ (e.g. the octahedral survivor of
   §11: $0.026 \to 0.051 \to 0.103$) is rejected as non-converged.

**Positive control.** The gate must not kill real chaos: the Halvorsen system
($\dot x = -ax - 4y - 4z - y^2$ cyclic, $a = 1.4$) passes cleanly through the same
machinery with $\lambda_1 = 0.632$ at $\Delta t = 0.005$ and $0.662$ at $0.0025$, and the
$C_n$/$S_4$ specimens' certified exponents are step-stable. The estimator and gate are
sound; the artifact is specific to coarse steps near strongly attracting equilibria.

---

## 10. The polyhedral irreducibility obstruction

The flagship goal was the first chaotic flows with $T$, $O$, $I$ rotation symmetry. The
result is the opposite and sharper: **within every low-degree equivariant ansatz searched,
no dt-robust strange attractor exists.** The evidence has five independent strands.

### 10.1 The structural argument (why there is no engine)

- The 3-D representations of $T, O, I$ are absolutely irreducible $\Rightarrow$ by Schur's
  lemma the equivariant linear part is $\rho\,\mathbb{1}$. No linear rotation, no
  saddle-focus, three equal eigenvalues at the origin.
- The origin is the *only* group-fixed point (no trivial subrepresentation), so a
  $G$-symmetric attractor must organize around the one point where the flow is a pure
  balloon.
- The rotational engine $b\,(X \times \nabla I)$ conserves both $r^2$ and $I$ exactly
  (§4.3), hence is integrable — motion along the closed intersection curves of spheres
  and invariant level sets. Chaos can therefore only come from the *interplay* of the
  gradient terms ($a$, $c$) with the radial terms; numerically that interplay always
  resolves into relaxation.
- When the radial well dominates, the dynamics is slaved to a thin spherical shell, and a
  flow on a 2-sphere cannot be chaotic (Poincaré–Bendixson). Escaping the shell requires
  exactly the radial–angular exchange the irreducible structure starves.

### 10.2 Polynomial campaigns (all with the §9 gate)

| campaign | ansatz / box | samples | dt-robust chaos |
|---|---|---|---|
| degree-3 T, original box | $(\rho,d,a,b)$ | 4 000 | 0 (all $\lambda_1 < -0.3$) |
| degree-3 T, high-radius box | $\rho \le 10$, $d \ge 0.05$, wide $a,b$ | 6 000 | 0 (blow-ups + sinks; survivors were §9 artifacts) |
| 5-dial T (complete degree-3, $+c$) | full box | 4 000 + 16 000 | 1 stage-1 nominee / 16 000; 0 survive |
| two-step-size search loop | 24 rounds × 400 | 9 600 | 0 |
| degree-4/5 widenings | $+e\,r^2\nabla(xyz)$, $+f\,r^2(x^3,y^3,z^3)$ | 8 000 + 8 000 + 12 000 | 0 |
| non-gradient rotational terms | $+g\,X{\times}(x^3,y^3,z^3)$ (deg 4), $+h\,\nabla(xyz){\times}(x^3,y^3,z^3)$ (deg 5), 9-dial | 10k + 10k + 12k | 0 |
| strong-rotation / weak-dissipation regime | $\|b\| \le 35$, $d \le 1.5$ | 30 000 | 0 (3–4 stage-1 nominees per 10k, all §9 artifacts) |
| octahedral family | $(\rho,d,a,b)$, wide | 12 000 | 0 |
| icosahedral family | $(\rho,d,a,b)$, $\nabla I_6$ engine | 10 000 (+ law run) | 0 |

### 10.3 The stability-targeted scan (Poincaré–Hopf logic)

Blind sampling could miss a tiny chaotic window, so the search was inverted: scan for
parameter sets where (i) the origin is a source, (ii) **every** equilibrium (found by
220-seed multi-start `fsolve`) has an unstable direction, (iii) the flow is bounded
(120 000 fine steps). In such a configuration no fixed point can be the attractor, so a
bounded orbit is *forced* onto a cycle or a strange set. Linear analysis of the on-axis
fixed points showed they become saddle-foci precisely for $c < 0$, so the scan was biased
there ($c \in [-3, -0.05]$, 6 000 scans). Result: 6 all-unstable-and-source
configurations, 5 bounded — **all limit cycles** (two-step-size $\lambda_1$ tests
negative). The forced non-equilibrium attractor exists; it is periodic, not strange.

### 10.4 The trigonometric (labyrinth) escape route

The only chaotic-flow mechanism in the literature with comparable symmetry is Thomas's
cyclically symmetric labyrinth ($\dot x = \sin y - bx$, cyclic), which needs no linear
engine. Its $T$-equivariant analogue was built from the trigonometric equivariants

$$
S = (\sin y \sin z,\ \sin z \sin x,\ \sin x \sin y), \quad
G_s = \nabla(\sin x \sin y \sin z), \quad
C = X \times S,
$$

in the family $\dot X = (\rho - \mu)X + aS + gG_s + b_c C$ (plus damped and cell-spanning
variants). Outcome across 24 000+ samples and a dedicated small-damping scan:

- exactly **one** dt-robust chaotic candidate
  $(\rho,\mu,a,g,b_c) = (1.749, 1.776, 3.372, 0.012, 2.911)$: QR spectra
  $(0.039, 0.003, -0.376)$ / $(0.049, -0.000, -0.384)$ / $(0.020, -0.005, -0.364)$ at
  $\Delta t = 0.01 / 0.005 / 0.0025$ — genuinely chaotic, dissipative
  ($\overline{\nabla\!\cdot\!F} = -0.335$), bounded;
- but **tiny and symmetry-broken**: extent span $\approx [0.40, 0.37, 0.17]$, centered
  $\approx 3$ units off-origin — a single cell of the $2\pi$-periodic labyrinth, far from
  the origin. The group's other eleven copies are separate conjugate attractors; no single
  $T$-invariant chaotic set exists.
- The deep trap, stated precisely: a $G$-symmetric attractor must enclose the unique
  $G$-fixed point (the origin), but near the origin $\sin \to$ its Taylor polynomial and
  the flow reduces to the non-chaotic polynomial case of §10.2. *The symmetric place is
  the calm place.*
- Adding dissipation strong enough to merge cells kills the chaos entirely (0 candidates
  in 16 000 damped samples; 0 cell-spanning chaotic orbits in the small-$b$ scan).

### 10.5 Statement of the conclusion

> **Across $\approx 1.5\times10^5$ parameter evaluations spanning the complete degree-3
> equivariant fields, degree-4/5 gradient and non-gradient extensions, stability-targeted
> configurations, and trigonometric labyrinth analogues — all under a step-size-robust
> Lyapunov gate — the tetrahedral, octahedral, and icosahedral equivariant flows exhibit
> no strange attractor. Bounded orbits relax to equilibria or (symmetric) limit cycles.
> The only genuine chaos found is weak ($\lambda_1 \approx 0.04$), spatially localized,
> and symmetry-broken. We advance this, together with the Schur/integrability structure of
> §10.1, as strong numerical evidence — not a proof — that low-degree equivariant 3-D
> flows under an irreducible polyhedral action resist chaos, and as a candidate
> explanation for the literature's gap: chaotic polyhedral maps are classical, chaotic
> polyhedral flows are absent.**

What the polyhedral groups *do* host, verified constructively: symmetric equilibria and
limit cycles — the rungs below chaos.

---

## 11. The chaos-vs-symmetry-group law

### 11.1 Protocol (`scripts/symmetry_law.py`)

Twelve groups (eleven in the original run, plus $D_3$ added with Trinacria), one rule each:

- **Same evaluator settings**: batched search at $\Delta t = 0.005$, transient 16 000,
  Lyapunov window 50 000 steps.
- **Same budget**: 12 000 uniform candidate evaluations per group ($I$: 4 000 — its
  60-direction degree-5 gradient costs ~50× per step; the obstruction is established as
  budget-independent by $T$ and $O$ at full budget).
- **Certification, not estimation**: the top 24 batch nominees per group are re-integrated
  with the continued-orbit fine integrator ($\Delta t = 0.0025$, 40 000-step transient,
  120 000-step window); survivors must exceed $\lambda_1 > 0.02$. The recorded comparable
  is the **maximum certified $\lambda_1$**; at that maximum a full QR spectrum and
  $D_{KY}$ are computed (finite-difference Jacobian tangent RK4).
- The batch **chaotic fraction** is recorded but de-emphasized (§11.3).

### 11.2 Results (`results/symmetry_law.csv`, figure `gallery/symmetry_law.png`)

| group | action | $|G|$ | $n_{\text{eval}}$ | max certified $\lambda_1$ | $D_{KY}$ at max | chaotic fraction | certified of top-24 |
|---|---|---:|---:|---:|---:|---:|---:|
| $C_2$ | reducible | 2 | 12 000 | 0.410 | 2.250 | 0.0176 | 24 |
| $C_3$ | reducible | 3 | 12 000 | 0.589 | 2.348 | 0.0195 | 24 |
| $C_4$ | reducible | 4 | 12 000 | 0.564 | 2.196 | 0.0051 | 23 |
| $S_4$ | reducible, improper | 4 | 12 000 | **0.628** | 2.298 | 0.0021 | 18 |
| $C_5$ | reducible | 5 | 12 000 | 0.509 | 2.266 | 0.0027 | 22 |
| $C_6$ | reducible | 6 | 12 000 | 0.272 | 2.352 | 0.0008 | 9 |
| $D_3$ | reducible, chiral dihedral | 6 | 12 000 | **0.689** | 2.490 | 0.0063 | 23 |
| $C_7$ | reducible | 7 | 12 000 | 0.326 | 2.315 | 0.0006 | 5 |
| $C_8$ | reducible | 8 | 12 000 | 0.176 | 2.415 | 0.0006 | 4 |
| $T$ | **irreducible** | 12 | 12 000 | **0.000** | — | 0.0014 | 0 |
| $O$ | **irreducible** | 24 | 12 000 | 0.026† | (1.94)† | 0.0003 | 2† |
| $I$ | **irreducible** | 60 | 4 000 | **0.000** | — | 0.0013 | 0 |

†The two octahedral "survivors" fail the convergence arm of the §9 gate
($\lambda_1 = 0.026 \to 0.051 \to 0.103$ as $\Delta t$ halves twice: non-converged), have
strongly anisotropic supports (span $\approx [0.15, 2.0, 2.0]$), sit off-center
($|\overline{X}| \approx 1.3$ — symmetry-broken), and carry $a \approx 0.001$ (the
$O$-signature term essentially off). They are reported in the CSV as measured but flagged
as **not** genuine $O$-symmetric attractors. The icosahedral nominee behaved identically
($0.018 \to 0.039 \to 0.050$, off-center, $a \approx -0.01$) and fell below threshold
under the standard certification.

**Reading.** Two regimes, split exactly at the reducible/irreducible boundary: the
$C_n$/$S_4$ groups host strong chaos with a **non-monotonic** dependence on group order
(0.18–0.63; the rotoreflection $S_4$ is the strongest of all, and $C_8 < C_2$), while the
polyhedral groups sit at the floor — the **irreducibility cliff**. The family maxima
exceed the named specimens' $\lambda_1$ (specimens were selected for beauty and novelty,
not maximal stretching).

**The $D_3$ row** (added with Trinacria, run with `PYTHONHASHSEED=0`, 753 s) is the largest
maximum in the table: 0.689, above $S_4$ (0.628) and twice $C_6$ (0.272), the cyclic group of
the same order. The maximizer was spot-checked outside the protocol ($T = 1000$): it is
centred (centroid within 0.6% of the extent), with $\rho$ and $\kappa$ residuals of 0.0135
against an identity floor of 0.0087 and a $z$-flip control of 0.027, so it reads as a
symmetric attractor rather than a broken one. Equal order therefore does not mean equal
chaos: $C_6$ and $D_3$ both have order 6 and differ by a factor of 2.5. One row is not a
trend, and the $D_3$ family has one more parameter than $C_n$ (8 against 7), so the
comparison carries the §11.3 caveats. The named specimen Trinacria ($\lambda_1 = 0.33$) was
selected for its symmetry evidence and look, not for maximal stretching.

### 11.3 Fairness caveats (stated, not buried)

The families have different parameter counts ($C_n$: 7, $S_4$ and $D_3$: 8, $T/O/I$: 4) and
different nonlinearity degrees, so the *chaotic fraction* is only loosely comparable
across families; the robust comparable is **max certified $\lambda_1$ under equal
evaluation budget**. The decline within $C_n$ at high $n$ partly reflects the fixed box
($\bar w^{n-1}$ changes character with $n$); the cliff at T/O/I does not — §10's targeted
searches rule out a mere box artifact there.

### 11.4 Null hypothesis

In the Letellier–Gilmore cover/image construction, an $n$-fold cover of a known attractor
has the *same* Lyapunov exponents as its image (the cover is locally identical to it), so
group-independence is the null. Our families are independently searched, not covers of a
shared image, so the measured group-dependence (both the non-monotonic wandering and the
cliff) is a property of the equivariant design space itself, not a cover artifact.

---

## 12. Numerical environment and known pitfalls

- Python 3.14 (Homebrew), NumPy pinned **2.3.5**, SciPy 1.17, pytest 9. A dysts install
  can silently downgrade NumPy to 2.1.x, which (a) renders 3-D matplotlib scatters pure
  black and (b) has a histogram fast-path bug (integer `bins` + `range=` miscounts
  edge-landing float64 values) — `atlas/fingerprint.py` uses explicit bin edges as the
  permanent fix. If renders come out black, check `numpy.__version__` first.
- Use `np.ptp(arr)`, not `arr.ptp()` (removed in NumPy 2).
- Batched integrations freeze culled members at 0 inside `np.errstate(all="ignore")` to
  keep overflow warnings from masking real signal.
- All stochastic steps (search, fingerprint sampling, multi-start seeds) use seeded
  `np.random.default_rng` for reproducibility.
- Long single-orbit certifications are CPU-bound Python; the law runs them per-group in
  parallel processes with BLAS threads pinned to 1 (except the matmul-heavy icosahedral
  evaluator, which prefers threaded BLAS).

---

## 13. Reproducibility map

| claim | code | artifact |
|---|---|---|
| $C_n$ family + QD atlas | `atlas/family.py`, `atlas/map_elites.py`, `scripts/run_atlas.py` | `results/atlas/archive.json` |
| Aurelia certification | `scripts/verify_chaos.py` | `results/verification.json` |
| Naiad / Cassiopea certification | `scripts/verify_naiad.py`, `scripts/verify_cassiopea.py` | `results/{naiad,cassiopea}_verification.json` |
| $S_4$ equivariance (math gate) | `tests/test_s4_equivariance.py` | pytest |
| Mobula package = family (cross-check) | `tests/test_mobula_matches_family.py` | pytest |
| Mobula certification (incl. symmetry residual) | `scripts/verify_mobula.py` | `results/mobula_verification.json` |
| $D_n$ equivariance and overgroup exclusion (math gate) | `atlas/family_d.py`, `tests/test_d_equivariance.py` | pytest |
| D3 nominee screen and selection | `scripts/search_d3.py` | `results/d3_candidates.json`, `results/d3_selection.json` |
| Trinacria package = family (cross-check) | `tests/test_trinacria_matches_family.py` | pytest |
| Trinacria certification (extended battery, §7.9) | `scripts/verify_trinacria.py` | `results/trinacria_verification.json` |
| shared basin protocol, all specimens | `scripts/basin_comparison.py` | `results/basin_comparison.json` |
| fingerprint metric + validation (noise floor) | `atlas/fingerprint.py`, `scripts/validate_fingerprint.py` | `results/fingerprint_validation.json`, `gallery/fingerprint_map.png` |
| Lyapunov convergence / error bars | `scripts/convergence_study.py` | `results/convergence.json`, `gallery/convergence.png` |
| Shilnikov / homoclinic hunt | `scripts/homoclinic_hunt.py` | `results/homoclinic_hunt.json`, `gallery/homoclinic.png` |
| T family + equivariance gate + obstruction note | `atlas/family_t.py`, `tests/test_t_equivariance.py` | pytest; docstring |
| O family + gate | `atlas/family_o.py`, `tests/test_o_equivariance.py` | pytest |
| I family (60 matrices, $I_6$) + gate | `atlas/family_i.py`, `tests/test_i_equivariance.py` | pytest |
| the law (protocol, certification, figure; D3 row with `PYTHONHASHSEED=0`) | `scripts/symmetry_law.py` | `results/law_*.json`, `results/symmetry_law.csv`, `gallery/symmetry_law.png` |
| renders / rotations | `scripts/render_*.py`, `scripts/render_rotations.py` | `gallery/` |
| literature verdicts | — | `docs/RELATED_WORK.md` |

Suite: `python -m pytest tests/` (34 tests green at the end of the first program; 138 after
the dihedral family and Trinacria).

---

## 14. Limitations

1. **The obstruction is numerical evidence over stated ansätze**, not a theorem: degree
   ≤ 5 polynomial equivariants (complete only through degree 3), one trigonometric
   family, $\sim 1.5\times10^5$ evaluations. A proof, or a counterexample at higher
   degree / different functional form, would each be significant.
2. **Novelty is catalog-relative** (135 systems), with a measured noise floor; it is not
   a literature claim.
3. **The law's chaotic-fraction column is not budget-fair across families** (§11.3);
   only the max-certified-$\lambda_1$ column is the headline comparable. The icosahedral
   row used a reduced budget (4 000) for cost reasons.
4. **Shilnikov chaos is evidenced, not proven** — near-connections at $7\times10^{-5}$
   are not a homoclinic orbit; rigorous verification (e.g. interval arithmetic) is open.
5. Certified exponents carry the error bars of the convergence study; values cited to 3
   decimals are stable to the window lengths used there.
6. **Trinacria's group is not a novelty claim.** Letellier and Gilmore's cover construction
   can lift any chaotic image system to a $D_3$-equivariant flow, and we could not read that
   work in full; the claim is one designed and certified specimen with exactly $D_3$, to our
   knowledge the first, never "the first dihedral attractor" (`docs/RELATED_WORK.md`,
   Part B, Claim 4).
7. **Trinacria has the smallest basin of the five** (about half of $[-2,2]^3$ escapes) and
   mixes slowly between its two halves (dwell up to about 40 time units), so its centring and
   $z$-statistics are only meaningful on the 64-orbit ensemble; two zero-expectation ensemble
   statistics sit at 2.4 and 2.8 SE.

---

## 15. References

The complete annotated bibliography with verdicts is `docs/RELATED_WORK.md`; the compiled
paper draft (`paper/main.tex`, `paper/references.bib`) carries full citations. Key
methodological sources:

- G. Benettin, L. Galgani, A. Giorgilli, J.-M. Strelcyn, *Lyapunov characteristic
  exponents for smooth dynamical systems...*, Meccanica 15 (1980) — the QR spectrum
  method.
- J. Kaplan, J. Yorke, *Chaotic behavior of multidimensional difference equations* (1979)
  — the dimension formula.
- M. Golubitsky, I. Stewart, D. Schaeffer, *Singularities and Groups in Bifurcation
  Theory II* (1988) — equivariant vector fields, the $\mathbb{Z}_n$ resonant coupling.
- M. Field, M. Golubitsky, *Symmetry in Chaos* (1992/2009) — symmetric chaotic *maps*.
- I. Melbourne, M. Dellnitz, M. Golubitsky, *The structure of symmetric attractors*,
  Arch. Rational Mech. Anal. 123 (1993) — admissibility theory.
- C. Letellier, R. Gilmore, *Covers and images of chaotic flows* (2001); R. Gilmore,
  C. Letellier, *The Symmetry of Chaos* (2007) — the cover/image null.
- R. Thomas, *Deterministic chaos seen in terms of feedback circuits* (1999) — the
  cyclically symmetric labyrinth flow.
- R. Osada et al., *Shape distributions*, ACM TOG 21 (2002) — the D2 descriptor.
- J.-B. Mouret, J. Clune, *Illuminating search spaces by mapping elites* (2015);
  J. Lehman, K. Stanley, novelty search — the QD lineage.
- J. C. Sprott, *Some simple chaotic flows* (1994); *Elegant Chaos* (2010) — the genre
  and its standards.
- W. Gilpin, *Chaos as an interpretable benchmark...* (2021) — the dysts catalog.
