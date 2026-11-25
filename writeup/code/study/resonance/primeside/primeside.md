# The one-point law read on primes: exact constants, effect sizes, Pythagorean rays, a draft registered test

Folder: `scratchpad/primeside/`. Nothing in the repository was changed. This round was stopped early because memory was
tight (coordinator's request). **No sum on the fresh range (x > 3·10⁷) was computed.** Labels: **proved**, **computed**,
**heuristic**, **conjecture**, **literature (opened)**, **NOT DONE**.

## 0. Status

| item | status |
|---|---|
| constants c, root numbers, L-values for 15 cases (`ps_constants.py`, output `ps_constants_out.txt`, `constants.json`) | **DONE** |
| independent checks of sieve, enumerations, eigenvalues (`ps_check_lib.py`, `ps_check_lib_out.txt`: 25/25 PASS) | **DONE** |
| predicted constants K (table §1) | **DONE** (from `constants.json`, formula of §2) |
| effect sizes (§3) | **DONE** (formulas, numbers) |
| Pythagorean-ray derivation (§4) | **DONE** (derivation); numerical check **NOT DONE** (sums for x ≤ 2·10⁷ exist in `sums_dry.jsonl`, never evaluated) |
| classical equivalent, what is provable (§5), prior art (§6) | **DONE** with the sources opened |
| prime sums on the seen range (dry run, `ps_sums.py dry`): `sums_dry.jsonl`, 82 s | **DONE** (computed, **not evaluated**) |
| `ps_freeze.py` (predictions + sha256), `ps_test.py` (readings, decisions, hash check) | **written, never run**; no frozen file, no hash |
| `registration.txt` | **DRAFT only** (no hash) |
| Δ on primes, lattice-point (theorem) control of the zones, Sym² beyond the draft | **NOT DONE** |

## 1. Predicted constants (DONE)

Statement for each case, with a(p) = Hecke eigenvalue normalised so that |a(p)| ≤ 2 (no L-functions needed to read it):

    P(x) := Σ_{p ≤ x} a(p)·ln p·e(−d·(p/Q)^{1/d})  ~  K·x^{(d+1)/(2d)}        (conjecture; e(y) = e^{2πiy})

| case | what a(p) is (elementary) | d | Q | ε | c | K | arg K | 1/c ("constant 1" reads) |
|---|---|---|---|---|---|---|---|---|
| A4 | p = a² + b², θ = arctan(b/a): 2cos 4θ; 0 if p ≡ 3 (4) | 2 | 4 | +1 | 2.026446 | 0.95528 − 0.95528i | −45° | 0.493 |
| A8 | 2cos 8θ | 2 | 4 | +1 | 1.266091 | 0.59684 − 0.59684i | −45° | 0.790 |
| A12 | 2cos 12θ | 2 | 4 | +1 | 1.220874 | 0.57553 − 0.57553i | −45° | 0.819 |
| A16 | 2cos 16θ | 2 | 4 | +1 | 1.470547 | 0.69322 − 0.69322i | −45° | 0.680 |
| A20 | 2cos 20θ | 2 | 4 | +1 | 0.821297 | 0.38716 − 0.38716i | −45° | 1.218 |
| B6 | p = a² + 3b², θ = arctan(√3b/a): 2cos 6θ; 0 if p ≡ 2 (3) | 2 | 3 | +1 | 2.370628 | −1.20086 + 1.20086i | +135° | 0.422 |
| B12 | 2cos 12θ | 2 | 3 | +1 | 1.829022 | 0.92650 − 0.92650i | −45° | 0.547 |
| C23 | (roots of x³ − x − 1 mod p) − 1 | 2 | 23 | +1 | 1.853352 | 0.56420 − 0.56420i | −45° | 0.540 |
| C31 | (roots of x³ + x − 1 mod p) − 1 | 2 | 31 | +1 | 1.879242 | 0.53095 − 0.53095i | −45° | 0.532 |
| C47a | 2 if p = x² + xy + 12y²; 2cos(2π/5) if p = 2x² ± xy + 6y²; 2cos(4π/5) if p = 3x² ± xy + 4y²; 0 if (−47/p) = −1 | 2 | 47 | +1 | 1.030813 | 0.26246 − 0.26246i | −45° | 0.970 |
| C47b | same with 2cos(2π/5) ↔ 2cos(4π/5) | 2 | 47 | +1 | 0.314982 | 0.08020 − 0.08020i | −45° | 3.175 |
| D11a | a_p(11a)/√p (point count) | 2 | 11 | +1 | 0.965595 | −0.35347 − 0.35347i | −135° | 1.036 |
| D37a | a_p(37a)/√p | 2 | 37 | −1 | 0.219691 | 0.05938 + 0.05938i | +45° | 4.552 |
| E11a (Sym² 11a) | a_p²/p − 1 (p ≠ 11), 1/11 at p = 11 | 3 | 121 | +1 | 1.944698 (reused) | 0.75727i | +90° | 0.514 |
| DD (Δ) | τ(p)/p^{11/2} | 2 | 1 | +1 | 2.360823 (reused) | 1.57388 + 1.57388i | +45° | 0.424 |

The first two exploratory runs agree with this table: A4 at 2·10⁷ read 0.956 − 0.973i against K = 0.955 − 0.955i; C23 at
3·10⁷ read 0.546 − 0.570i against 0.564 − 0.564i. The dry-run file reproduces the exploratory A4 value at x = 10⁴
(0.8191 − 1.0078i; the old script left out p = 2, whose term is 0.69 in these units: 0.8186 − 1.0082i).

How c was computed (computed; `ps_constants_out.txt`): c = B/L(1, Ad) with B = Π_p E_p absolutely convergent (p ≤ 10⁷ for
A–C, 10⁶ for D) and L(1, Ad) from PARI; checked against the conditionally convergent singular series Π_{p≤P} G_p(1/p) at
P = 10⁴…10⁷ (agreement 10⁻⁴–10⁻⁵ at 10⁷, e.g. A4: 2.02633 at 10⁷ against 2.02645; C47b: 0.31495 against 0.31498), and
C23/C31 against the round-151 constants files (1.8533507, 1.8792428: agree to 2·10⁻⁶). Every L-function passed PARI's
functional-equation check (lfuncheckfeq −63 … −95); every root number is +1 except 37a (−1). L(1, ψ²) for A4…A20:
1.501229, 2.131617, 2.264202, 1.997278, 2.862556; B6, B12: 1.489608, 1.858872; L(1, ρ) for C23, C31: 0.723999/L(1,χ),
etc. (`constants.json`). E11a and DD reuse c from `universal/constants_S11a.json`, `constants_D12.json` (B over p ≤ 3000).

## 2. Derivations (proved unless marked)

**K.** From the prime law of general_law.md §2.1 / degd_proof Corollary G2.2:
K = ε·c·e^{−i(φ₀+π/4)}·(2√d/(d+1))·Q^{−1/(2d)}, φ₀ = (π/2)Σμ_j − dπ/4. Weight k, degree 2: φ₀ + π/4 = (2k−1)π/4. For the
Q(i) cases k = m + 1 ≡ 1 (4), so the phase is e^{−iπ/4} for every m, and (2√2/3)·4^{−1/4} = 2/3: K = (2/3)εc·e^{−iπ/4}.
B6 (k = 7): e^{−13πi/4} = e^{+3πi/4} (+135°). E11a: μ = {1,1,2}, φ₀ = 5π/4, e^{−3πi/2} = i: K = i·c·(√3/2)·121^{−1/6}.

**Local factors of c = Π_p c_p.** Good p (χ(p) = ±1): c_p = 1 − (λ_p² − 1)/p; inert p (λ = 0, Satake ±1): 1 + 1/p.
Ramified prime of a CM or dihedral form (one Satake parameter α, |α| = 1: F_p = 1 − X): **c_p = F_p(1/p)/(1 − 1/p) = 1
exactly** (so "no factor at 2" in the exploratory run was right). λ(2) = (−1)^{m/4} for Q(i), λ(3) = (−1)^{m/6} for
Q(√−3), λ(D) = 1 for the dihedral forms (the prime above D is principal). Multiplicative p (11a, 37a): F_p = 1 − X/p,
c_p = 1 + 1/p. Absolutely convergent form: E_p := c_p·L_p(1, Ad) = (1 − (λ²−1)X)/((1−X)(1 − (λ²−2)X + X²)), X = 1/p, at
every good p (split and inert); E_p = p/(p − 1) at the ramified prime; E_N = (1 + 1/N)(1 − 1/N²)^{−1} for 11a, 37a (PARI's
Sym² factor). Ad f = χ_D ⊕ Ind(ψ²) for f = θ_ψ, so **c = B/(L(1, χ_D)·L(1, ψ²))**, L(1, χ_{−4}) = π/4,
L(1, χ_{−3}) = π/(3√3); ψ² = (z/|z|)^{2m}, the CM form of weight 2m + 1 (PARI mfinit, analytic s = 1 ↔ arithmetic
s = m + 1); for h = 3, Ind ψ² = Ind ψ̄ = ρ; for D = −47, Ad ρ_a = χ ⊕ ρ_b.

**Ray form (proved, algebra).** Π_split(1 − (λ²−1)/p)·Π_inert(1 + 1/p) = Π_{p∤D}(1 − χ_D(p)/p)·Π_split(1 − (λ²−2)/(p−1)), so

    c = (1/L(1, χ_D))·Π_{p split} (1 − λ_{ψ²}(p)/(p − 1)),     Q(i): c = (4/π)·Π_{p≡1(4)} (1 − 2cos(2mθ_p)/(p − 1)).

(Partial products: A4 2.02637 at 10⁷.) The factor 4/π is the density of Gaussian primes per lattice point; the product is
the Pythagorean-ray part (§4).

**First-order kernel factor (proved per term, Lemma 5 of degd_proof).** With the bare kernel the integer resonance has
density (1 − iκ/t), t = 2π(x/Q)^{1/d}, κ = (d² − 1)/(24d) + ½Σ(μ_j − μ_j²) (= 1/16 − (k−1)²/4 in degree 2). Hence the main
term M(x) = K[x^{3/4} − (3iκ√Q/2π)x^{1/4}] (degree 2). For Q(i) this is the Bessel correction (4m² − 1)/(16πr) of
J_m(2πr) (checked by hand: same constant); relative size ≤ 6·10⁻³ at x ≥ 3·10⁷ (m = 20), 1.4 % for Δ at 10⁶. Included in
`ps_lib.main_term`.

**Prime powers (heuristic, but the two models agree).** For N a perfect square (Q(i): N = 4; Δ: N = 1) the squares are
coherent: e(−2√(p²/N)) = 1. Rule used ("λΛ rule"): Σ_n λ(n)Λ(n)·kernel has no x^{1/2} term, so the primes-only sum is
P(x) = c·M(x) − Σ_{p^k ≤ x, k ≥ 2} λ(p^k)·ln p·e(−d(p^k/Q)^{1/d}). For Q(i) λ(p²) has mean 1 (λ(p²) = λ_p² − χ(p)), so the
term is −√x·(1 + o(1)) (real, the same for every m): 0.4 % of the main term at 10⁹, 1.5 % at 2·10⁷. The 2-D (Gaussian)
Hardy–Littlewood model gives the same: the inert q ≡ 3 (4) are Gaussian primes on the axes, the almost-periodic density
C(b) of primes a² + b² extended to b = 0 is C(0) = 1, which equals their weighted density, so Σ over all Gaussian prime
ideals (= Σ_n Λ_f(n) up to o(√x)) has no √x term. For level 1 (Δ) the rule gives P ≈ cM (mean λ(p²) = 0) and a −√T
term in the zero sums, which is what the round-151 review measured. The rule is registered with the alternative "no
prime-power term" as a secondary reading; at the expected precision (~2 % per case) it is not decisive.

## 3. The effect in elementary words (DONE)

(i) **Local correlation.** Over primes near x (weighted by ln p), the mean of a(p)·e(−2√(p/Q)) is
((d+1)/2d)·K·x^{(1−d)/(2d)} = (√2/2)·εc·e^{−iφ}·(Qx)^{−1/4} in degree 2 — a bias of relative size x^{−1/4}, against
x^{−1/2} for the classical prime races (Chebyshev bias ~ x^{−1/2}: 10⁻², 10⁻³, 3·10⁻⁵ at these x).

| case | x = 10⁴ | 10⁶ | 10⁹ |
|---|---|---|---|
| A4 | 0.101 | 0.032 | 0.0057 |
| A8 / A12 / A16 / A20 | 0.063 / 0.061 / 0.074 / 0.041 | 0.020 / 0.019 / 0.023 / 0.013 | 0.0036 / 0.0034 / 0.0041 / 0.0023 |
| B6 / B12 | 0.127 / 0.098 | 0.040 / 0.031 | 0.0072 / 0.0055 |
| C23 / C31 | 0.060 / 0.056 | 0.019 / 0.018 | 0.0034 / 0.0032 |
| C47a / C47b | 0.028 / 0.0085 | 0.0088 / 0.0027 | 0.0016 / 0.00048 |
| D11a / D37a | 0.038 / 0.0063 | 0.012 / 0.0020 | 0.0021 / 0.00035 |
| E11a (degree 3: (2/3)|K|x^{−1/3}) | 0.023 | 0.0050 | 0.00050 |
| Δ | 0.167 | 0.053 | 0.0094 |

(ii) **Gaussian primes.** For a Gaussian prime z at radius r = |z|, let u = 4·arg z − 2π|z| (mod 2π). From the law (A4,
ε = +1): (1/4)Σ_{z} ln N(z)·e^{iu} = P(x), dP = (3/4)K x^{−1/4}dx and Σ_shell ln N = 4dx give mean e^{iu} =
(c/2)e^{−iπ/4}/√r, so the density of u is **(1/2π)[1 + (c/√r)·cos(u + π/4)]**, c = 2.026 (factor checked). For all lattice
points (theorem, constant 1; π dx points per shell): (1/2π)[1 + ((4/π)/√r)·cos(u + π/4)]. Primes carry the lattice's own
resonance amplified by cπ/4 = 1.59 (m = 4). Sizes: amplitude c/√r = 20 % at r = 100 (N z ≈ 10⁴), 6.4 % at r = 1000,
1.1 % at r = 31,623 (N z = 10⁹); lattice points 12.7 %, 4.0 %, 0.7 %. For m = 8, 12, 16, 20 replace c by 1.266, 1.221,
1.471, 0.821 (ratios 0.99, 0.96, 1.15, 0.65). Where the excess sits: §4.

## 4. Pythagorean rays (derivation DONE; numerical check NOT DONE)

**Lattice points (proved by Poisson).** Σ_z f(z) = Σ_{w∈ℤ²} f̂(w), f(z) = e^{im arg z}e^{−2πi|z|}. The phase |z| + ⟨z, w⟩ is
stationary only for |w| = 1, on the ray z ∥ −w; Jacobi–Anger gives f̂(w) = e^{im arg w}·2πi^m∫ρJ_m(2πρ)e^{−2πiρ}dρ, whose
non-oscillating part is e^{−iπ/4}√ρ(1 + i(4m²−1)/(16πρ) + …). The only unit vectors of ℤ² are ±1, ±i: the integer effect
sits in the four **Fresnel zones |b| ≲ √a around the axes** (phase −πb²/a there, since e^{−2πia} = 1).

**Primes (Hardy–Littlewood, heuristic in the sum over moduli; each modulus exact).** Model Λ(a² + b²) ≈ Σ_q μ(q)c_q(n)/φ(q).
For modulus q the lattice is periodic mod q and Poisson runs over the dual (1/q)ℤ²: resonance needs |w| = q, w ∈ ℤ². The
two-dimensional Gauss sum Σ_{z₀ mod q} e((r|z₀|² + ⟨z₀, w⟩)/q) = q·χ₋₄(q) for odd prime q when |w|² ≡ 0 (q), and the mod-2
sum vanishes. So modulus q contributes μ(q)χ₋₄(q)/q × (1/4)Σ_{|w|=q} e^{im arg w}. For a split p = ππ̄ the vectors of
length p are units × {p, π², π̄²}: **besides the axes, modulus p resonates along the directions of π² and π̄², i.e. of the
Pythagorean triple (a² − b², 2ab, p)** (p = 5: 3 + 4i and 4 + 3i). The dual lattice of a class mod q contains unit vectors
exactly for these directions. Summing over moduli: every direction φ of a primitive triple with hypotenuse h carries

    A_h = (4/π)·Π_{p | h} (−1/(p − 1))   if h is squarefree,     A_h = 0 otherwise,

and the axes carry A = 4/π; c = Σ_rays A·e^{imφ} = (4/π)Π_{p≡1(4)}(1 − 2cos(2mθ_p)/(p−1)) (matches §2). Elementary check of
one weight: on the ray through 3 + 4i, |z| ≈ t + s²/2t with t = (3a + 4b)/5, so e^{−2πi|z|} has the factor e(−(3a+4b)/5);
3a + 4b mod 5 is uniform over {0,…,4} for lattice points (mean of e(−·/5) = 0) and uniform over {1,2,3,4} for Gaussian
primes (3a + 4b ≡ 0 iff (2+i) | z), mean −1/4: weight (4/π)(−1/4) = −1/π. Predictions with a sharp sign:
- the ray through 3 + 4i (53.13°) and 4 + 3i: A = −1/π each, i.e. −0.25 of the axis amplitude, phase e^{±im·53.13°};
  for m = 4 the two rays add (2/π)·0.843 = +0.537 to c = 2.026 (26 %);
- **7 + 24i (hypotenuse 25 = 5²): A = 0** — a primitive triple whose ray carries no prime resonance (μ(25) = 0);
- 33 + 56i and 63 + 16i (65 = 5·13): A = +(4/π)/48;
- lattice points: 0 on every Pythagorean ray.
Signal estimate (m = 0, Gaussian window w(u) = e^{−πu²/4} in u = s/√t, β = 2, fixed in `ps_lib`): on primes to 2·10⁷ the
3 + 4i zone sum is ≈ (1/π)|F|(2/3)x^{3/4} ≈ 6·10⁴ against a random-phase sd ≈ 4·10³, i.e. about 15 sd. The sums for the
windows axis, 3+4i, 4+3i, 5+12i, 12+5i, 8+15i, 15+8i, 20+21i, 21+20i, 7+24i, 24+7i, 12+35i, 35+12i, 9+40i, 40+9i and seven
control angles (0.35, 0.55, 0.70, 0.90, 1.05, 1.25 rad, π/4), m = 0 and 4, **are computed for (6·10⁵, 2·10⁷] in
`sums_dry.jsonl`; the ray-model predictions (`ps_lib.zone_prediction`, all primitive directions with h ≤ 5000, overlaps of
neighbouring zones included) are coded in `ps_freeze.py`; the comparison was NOT run.** Verdict on the rays: derived, not
yet verified.

## 5. Classical equivalent and what is provable (literature opened: Huang–Liu–Rudnick arXiv:1903.04005 intro; Stucky arXiv:2008.11325 Thm 1.1)

The axis part of the resonance is: Σ_{a ≤ R}Σ_b Λ(a² + b²)·e^{−2πi(√(a²+b²) − a)}·w(b/√a) ~ (4/π)·F_w·(2/3)R^{3/2},
F_w = ∫w(u)e^{−iπu²}du — primes a² + b² just above the square a², at a distance b² with b ≲ √a, counted with the Fresnel
weight e(−b²/2a). It is the Bateman–Horn density of x² + b² averaged over b (mean of the singular series over b =
Π_{p≡1}(1 − 1/p)Π_{p≡3}(1 + 1/p) = 4/π, computed above) — for fixed b this is Landau-type and open, but averaged over the
zone it is a smooth count of Gaussian primes in sectors of angular width N^{−1/4}.
- Ricci: asymptotic for every sector of angular width > X^{−3/10+ε} (unconditional); Stucky (2021): polar boxes of area
  δx^θ ≥ x^{7/10+ε}, θ > 7/10, with smooth weights; under GRH width X^{−1/2+ε}; almost all sectors of width X^{−3/5}
  (Huang–Liu–Rudnick). The axis zone has width X^{−1/4} > X^{−3/10}: **the axis part is within reach of Ricci/Stucky**
  (decompose in dyadic norm shells and polar boxes of width X^{−1/4−0.02}; the smoothing step was not written: a sketch).
- Each single Pythagorean ray needs the same count in a residue class mod the Gaussian prime(s) dividing h (Hecke
  characters of fixed modulus times angular characters): the method extends, but no statement with congruence
  conditions was found in the two sources opened.
- **What blocks the full law:** the sum over all rays — uniformity in the hypotenuse h up to ≈ x^{1/4} (beyond that the
  zones of width √t overlap) and the region between zones (the minor arcs). Per-ray GRH errors x^{1/2+ε} summed over the
  ≍ H rays exceed x^{3/4} for H ≈ x^{1/4}; this is the minor-arc bound o(x^{3/4}) of STUDY.md (known x^{5/6}).
- Baier–Zhao (n² + k on average over k ≤ K, K ≥ N^{1+ε}) does not apply: here k = b² is a sparse set of length N^{1/2}.
  Fouvry–Iwaniec (b prime), Friedlander–Iwaniec (a² + b⁴), Lam–Schindler–Xiao, Merikoski, Pratt (thin sets of b):
  **not opened** (no time); they concern existence/asymptotics for thin sets of b, not a Fresnel-weighted zone.

## 6. Prior art (searches with mathematical terms only)

| statement | verdict | basis |
|---|---|---|
| prime-side law for one form, with c = Π G_p(1/p) | NOT FOUND | Zubrilina 2310.07681 (opened via ar5iv): family over levels, trace-formula density, "the visible text does not discuss individual forms"; Cowan 2408.12723 (opened): ratios-conjecture murmurations, families only; no singular series or B/L(1, Sym²) |
| fact 2054 (fixed level, harmonic weights) = murmuration density integrated against the kernel | KNOWN IN SUBSTANCE | the Petersson c = N term gives the density ∝ J_{k−1}(4π√(p/N)), whose product with e(−2√(p/N)) has a non-oscillating part; Kuan–Lesesvre 2507.11418 (abstract) prove murmurations with the Petersson formula in the weight aspect; Bober–Booker–Lee–Lowry-Duda (level 1, weight aspect; search snippet only) |
| integer resonance Σλ(n)e(−α n^{1/d}) | KNOWN | Kaczorowski–Perelli, Twists and resonance I/II (1304.4734, 1507.07177; abstracts); Ren–Ye |
| von Mangoldt version of the resonance | NOT FOUND | the K–P abstracts do not mention Λ |
| Iwaniec–Luo–Sarnak 2000, Zhao 2006 | not re-opened here | the family note of the repo opened ILS; Zhao's x^{5/6} as quoted in STUDY.md |
| Gaussian primes with a radial phase e(|z|), Pythagorean rays | NOT FOUND | one search; sources on angles of Gaussian primes (Rudnick–Waxman, HLR) treat sectors, not radial phases |

## 7. Draft of the registered test (scripts written; freeze NOT DONE)

- Statistic (primary): smooth sum S_W = Σ_p a(p) ln p·p^{−1/2+1/(2d)}·e(−d(p/Q)^{1/d})·W(p^{1/d}), W a Hann window in
  v = p^{1/d} over the range; main term (3/2)K∫W(1 − iκ√Q/(2πv))dv (degree 2), minus the prime-power term. Reading
  ρ = (S_W + PP_W)/MT_W, expected 1 + 0i; "constant 1" gives 1/c (table §1), "no effect" 0. Sharp increments reported too.
  Which has the smaller fluctuation: in the random-prime model the smooth relative error is 1.15× the sharp one; not
  measured (dry run not evaluated).
- Errors: rms of 8 control frequencies (factors 0.6713 … 1.3711, α² far from small rationals) on the same primes, per
  series; 8 half-overlapping sub-windows as a check. Expected relative error per component on the fresh range (random
  model, smooth): ≈ 440/(19,600|K|): 1.7 % (A4), 2.8 % (C23), 20 % (C47b); D/E on (3·10⁷, 3·10⁸]; Δ on (3·10⁵, 2·10⁶] ≈ 5 %.
- Ranges: fresh (3·10⁷, 10⁹] for A, B, C; (3·10⁷, 3·10⁸] for 11a, 37a, Sym² 11a (ellap ≈ 3–5·10⁻⁵ s/prime, ≈ 500 s per
  curve); (3·10⁵, 2·10⁶] for Δ (mfcoefs; ramanujantau is 0.12 s/prime at 10⁶). Dry: (6·10⁵, 2·10⁷] (seen), Δ (10⁴, 3·10⁵].
- Timing: the whole dry run took 82 s. Fresh estimate: A–C ≈ 50–80 chunks of 2·10⁷ at a few s each (5–10 min), D ≈ 2 × 500 s,
  Δ < 1 min: three or four bounded runs of `ps_sums.py fresh 540`.
- Decisions (draft, `decisions.json`): pooled slope within 1 ± 0.05 and 3 errors; each family within 3 errors; at most one
  series beyond 3 errors; law beats constant 1 and no effect by Δχ² ≥ 25; control rms in [0.6, 1.6] errors; Pythagorean
  windows pooled within 4 errors, zero-weight/control windows within 4 sd. Falsifier: pooled slope outside 1 ± 0.10, or
  χ²(law) > χ²(constant 1).

## 8. What remains (exact commands; run from the folder, each under the memory cap)

    P=.venv/bin/python; cd <scratchpad>/primeside
    1. review decisions.json, then freeze:   (ulimit -v 1500000; timeout 600 nice $P -I -B ps_freeze.py)
       -> predictions_primeside.json, freeze_sha256.txt  (ps_freeze.py was never run: expect small bugs)
    2. put the hash into ps_test.py (FROZEN = '...'), then the dry readings (not blind):
       (ulimit -v 1500000; timeout 600 nice $P -I -B ps_test.py dry)
    3. finish registration.txt with the hash and the dry-run numbers; commit before step 4
    4. fresh sums, repeat until ALL DONE:    (ulimit -v 1500000; timeout 600 nice $P -I -B ps_sums.py fresh 540)
    5. (ulimit -v 1500000; timeout 600 nice $P -I -B ps_test.py fresh)
    Not written: the lattice-point (theorem) control of the zone sums; a referee pass.

## 9. Doubts

1. The λΛ rule for prime powers (§2) is heuristic; it flips the sign of the √x term relative to "primes only = cM" for
   Q(i). The difference is 0.4–1.5 % of the main term.
2. The ray model ignores e^{im(arg z − φ_h)} inside a zone (relative m²/r: 0.3 % for m = 4) and the zones' radial cutoff.
3. ps_freeze.py, ps_test.py and the zone-prediction code are untested; ps_sums.py ran on the dry range only.
4. E11a and Δ constants are reused (B over p ≤ 3000: ~10⁻⁴); D/E fresh range is shorter (to 3·10⁸), so 37a (|K| = 0.084)
   will read with ≈ 40 % error.
5. Δ to 2·10⁶ via mfcoefs needs ~0.3–0.6 GB of PARI stack (not tested above 10⁶).
6. The dry-run file includes chunks from 0; the readings use only (6·10⁵, 2·10⁷]. The fresh range has never been touched.
