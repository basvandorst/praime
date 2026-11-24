# The echo law: the error of the prime number theorem against its own tail

> **In short:** For E(x) = ψ(x) − x, the logarithmic mean of E(ax)E(x)/x exists under RH and equals Σ m(λ)² a^(½+iλ)/(¼ + λ²) over distinct ordinates. With simple zeros it equals −a ∫_a^∞ E₀(t) dt/t², E₀ the oscillating part: the error's autocorrelation at ratio a is its own tail integral. Proof: (1) a diagonal argument for almost periodic means; (2) the Riesz-mean explicit formula. Conversely, RH with simple zeros holds exactly when the upper mean of E(x)²/x is at most 2 + γ − ln 4π. The law holds for ψ(x, χ), with a cross law; it fails for μ.

Proof note, companion to `writeup/code/study/variance/D_delta_proof.md` (cited below as **[PN]**). Every number in this note is recomputed and asserted by one script, [check_echo.py](check_echo.py), with its output in [check_echo_output.txt](check_echo_output.txt) (187 checks, 0 failed, about 1 minute). The script imports nothing from the repo; its repo input is the zero table `writeup/data/riemann/zeros_40000_flint.txt`. Zeros of L(s, χ₋₄) and L(s, χ₋₃) come from PARI (`lfunzeros`), the primes to 10⁸ from the script's own sieve. Each step carries a label: **proved** (a complete argument from the cited inputs), **computed** (numerical, with the error stated), **heuristic**, or **conjecture**.

---

## 1. Statements

**Notation.** γ is Euler's constant, γ₁ = −0.0728158… the first Stieltjes constant. A nontrivial zero of ζ is ρ = ½ + iγ_ρ (counted with multiplicity); λ runs over the *distinct* ordinates, with multiplicity m(λ). A height is written H. As in [PN] the trivial-zero function is always written with its argument:

```
E(x)  = ψ(x) − x,                 E₀(x) = E(x) + ln 2π + ½ ln(1 − x⁻²)   (x > 1),
T(a)  = Σ_{k≥1} a^(−2k)/(2k(2k+1)) = 1 − ½ ln(1 − a⁻²) − a·artanh(1/a),     T(1) = 1 − ln 2,
J(a)  = Σ_ρ a^ρ/(ρ(1 − ρ))        (absolutely convergent, real, no hypothesis),
C     = J(1) = 2 + γ − ln 4π = 0.0461914179322421,        B = 1 − γ − ln 2π = −1.41509273131,
K₂    = Σ_ρ 1/(ρ(1 − ρ))² = 6 + 2γ + 2γ² + 4γ₁ − π²/4 − 2 ln 4π = 7.42012728749297·10⁻⁵.
```

By von Mangoldt's formula, E₀(x) = −lim_H Σ_{|γ_ρ|≤H} x^ρ/ρ for x > 1 not a prime power: E₀ is E without its constant −ln 2π and its trivial-zero part −½ ln(1 − x⁻²). For a > 0 put F_a(y) = E(a e^y) E(e^y) e^(−y). Then (1/L) ∫_1^(e^L) E(ax)E(x) dx/x² = (1/L) ∫_0^L F_a(y) dy: the *logarithmic mean* of E(ax)E(x)/x. A function F has *uniform mean* M if sup_{c≥c₀} |(1/L) ∫_c^(c+L) F − M| → 0 as L → ∞. We write ⟨G(x)⟩ for the uniform logarithmic mean of G, i.e. the uniform mean of G(e^y).

**Theorem 1 (the echo law; RH).** Assume RH and fix a > 0.

- (a) *Existence.* sup_{c≥c_a, L≥1} (1/L) ∫_c^(c+L) |F_a| < ∞ (c_a = max(0, −ln a)), and the uniform mean A(a) = ⟨E(ax)E(x)/x⟩ exists. In particular (1/L) ∫_1^(e^L) E(ax)E(x) dx/x² → A(a) for a ≥ 1. For every bounded Riemann-integrable φ on [0, 1] (for example the Hann weight), (1/L) ∫_c^(c+L) φ((y − c)/L) F_a(y) dy → A(a) ∫_0^1 φ, uniformly in c.
- (b) *Value.* A(a) = Σ_λ m(λ)² a^(½+iλ)/(¼ + λ²), absolutely convergent, and A(1/a) = A(a)/a.
- (c) *Constants.* Replacing E by E + u, with ∫_1^∞ u(x)² dx/x² < ∞, in either factor does not change A(a). So ⟨E₀(ax)E₀(x)/x⟩ = A(a): the constant ln 2π and the trivial zeros contribute nothing to the left side.
- (d) *Closed form (RH and simple zeros).* If all zeros are simple, then for a ≥ 1

  ```
  A(a) = J(a) = −a ∫_a^∞ E₀(t) dt/t²  =  a [1 + γ + ∫_1^a E(t) dt/t²] − ln 2π + T(a).          (1.1)
  ```

  In general A(a) = J(a) + Σ_λ m(λ)(m(λ) − 1) a^(½+iλ)/(¼ + λ²).
- (e) *No hypothesis.* The three right-hand expressions in (1.1) are equal for every a ≥ 1 without RH; ∫_1^∞ E(t) dt/t² converges and equals −1 − γ.

So, under RH with simple zeros, the oscillating part of the error obeys the closed equation

```
⟨E₀(ax) E₀(x)/x⟩ = −a ∫_a^∞ E₀(t) dt/t²          (a ≥ 1),                                    (SC)
```

quadratic on the left, linear on the right ("echo law"). In terms of E itself the right side is −a ∫_a^∞ E(t) dt/t² − ln 2π + T(a): the constant and the trivial zeros, invisible on the left, enter the right side because it is linear in E.

**Theorem 2 (converse; no hypothesis).** Let I(X) = ∫_1^X E(x)² dx/x².

- (i) If ζ has a zero with real part β > ½, then I(X) ≠ O(X^θ) for every θ < 2β − 1; in particular limsup I(X)/ln X = ∞.
- (ii) RH ⟺ limsup_{X→∞} I(X)/ln X < ∞ ⟺ the logarithmic mean A(1) of E(x)²/x exists as a finite limit.
- (iii) **RH and all zeros simple ⟺ limsup_{X→∞} I(X)/ln X ≤ 2 + γ − ln 4π** ⟺ A(1) exists and equals 2 + γ − ln 4π.
- (iv) Assume RH. Then A(a) − J(a) = √a·D(ln a), D(b) = Σ_λ m(m − 1) e^(iλb)/(¼ + λ²), a Bohr almost periodic function with coefficients ≥ 0. These are equivalent: all zeros simple; D(0) = 0; D = 0 on a half-line [b₁, ∞). If only finitely many zeros are multiple, D = 0 on any nonempty open interval already forces all zeros to be simple.

**Theorem 3 (what the equation determines).**

- (a) *Bootstrap (proved).* Let Λ ⊂ ℝ be countable and symmetric with #(Λ ∩ [t, t + 1]) ≪ ln(|t| + 2), and c: Λ → ℂ∖{0} with c(−λ) = c̄(λ) and |c(λ)| ≤ K/(1 + |λ|). Let Q_c be the L²_loc limit of Σ_{|λ|≤H} c(λ) e^(iλy) and R(x) = −√x Q_c(ln x). Then for every a ≥ 1, ⟨R(ax)R(x)/x⟩ = √a Σ |c(λ)|² a^(iλ) and −a ∫_a^∞ R(t) dt/t² = √a Σ c(λ) a^(iλ)/(½ − iλ), and these are equivalent: (SC) holds for R for all a in some half-line [a₁, ∞); c(λ) = 1/(½ + iλ) for every λ ∈ Λ.
- (b) *Consequence.* (SC) does not determine the zeros: R = −Σ_{λ∈Λ} x^(½+iλ)/(½ + iλ) solves it for *every* admissible Λ, so deleting zeros from the explicit formula keeps a solution. What (SC) fixes is the shape of each term: every frequency carries the amplitude 1/ρ exactly, in modulus *and* phase, i.e. it is a simple zero on the critical line. The left side alone sees only |c|²; the right side is linear in c and pins the phase. E₀ is in the class (Λ the ordinates, c = m/ρ), so (SC) holds for E₀ iff m ≡ 1, in line with Theorem 2.
- (c) *Primes ⟷ echo (proved; RH, simple zeros).* For every a₀ > 1, the restriction of A to [1, a₀] and the restriction of ψ to [1, a₀] determine each other.

**Theorem 4 (characters; GRH).** Let χ be primitive mod q > 1, χ(−1) = (−1)^κ, κ ∈ {0, 1}, with E_χ(x) = ψ(x, χ) (no main term), and assume GRH for L(s, χ). Then ⟨E_χ(ax) Ē_χ(x)/x⟩ = Σ_λ m_χ(λ)² a^(½+iλ)/(¼ + λ²), in the senses of Theorem 1(a); with simple zeros it is

```
J_χ(a) = Σ_ρ a^ρ/(ρ(1 − ρ)) = −a ∫_a^∞ E_χ,₀(t) dt/t²
       = a [L′/L(1, χ) + ∫_1^a ψ(t, χ) dt/t²] − L′/L(0, χ) + T_o(a)          (χ odd),
       = a [L′/L(1, χ) + ∫_1^a ψ(t, χ) dt/t²] − ln a − 1 − B₀(χ) + T(a)      (χ even),          (1.2)
```

with E_χ,₀ = ψ(x, χ) − A_χ(x) the oscillating part (A_χ in §5), T_o(a) = Σ_{k≥1} a^(1−2k)/((2k − 1)2k) = artanh(1/a) + (a/2) ln(1 − a⁻²), T_o(1) = ln 2, and B₀ = lim_{s→0}(L′/L(s, χ) − 1/s). The slope of J_χ jumps by χ(n)Λ(n)/n at each prime power n (by Re χ(n)Λ(n)/n for the real part). *Cross law:* for primitive χ ≠ χ′ (χ′ may be trivial, with E for ζ), under GRH for both, ⟨E_χ(ax) Ē_χ′(x)/x⟩ = Σ_{λ common} m_χ(λ) m_χ′(λ) a^(½+iλ)/(¼ + λ²); it vanishes for all a iff the two L-functions have no common zero, and at a = 1 it is ≥ 0 with equality iff they have no common zero. The criterion of Theorem 2(iii) holds with C_χ = J_χ(1) in place of C.

**Where the hypotheses are used.** RH (GRH) gives |x^ρ| = √x, real frequencies and |ρ|² = ρ(1 − ρ): Theorems 1(a)–(d), 2(iv), 3(b) for E₀, 3(c), 4. Simple zeros are used only to identify the mean with J (Theorem 1(d), 3(c), the J_χ form in Theorem 4, the race formula); in Theorem 2(iii) they are the conclusion. Not used anywhere: linear independence of ordinates (LI), except in §6.2, where it is the stated hypothesis. Theorem 1(e), Theorem 2(i)–(iii) and the identities (1.1), (1.2) between the right-hand sides use no hypothesis.

**Cited inputs.** [PN] §1 (I1)–(I3) (explicit formula, zero counting, growth of ζ′/ζ, Davenport ch. 15–19); [PN] Lemma 1 (pair sum), Lemma 2 (mean square; its proof uses only |β(λ)| ≤ K w(λ), w(λ) = m(λ)/(1 + |λ|), so it holds for any such β), the Riesz-mean formula (4.2), Step 3 of §4.1, §5 (characters) and C7(a). Bohr's theory of almost periodic functions (uniqueness of Fourier–Bohr coefficients, relative density of ε-almost periods).

---

## 2. Proof of Theorem 1

Assume RH. Let β(λ) = m(λ)/(½ + iλ). Since (1 + |λ|)² ≤ 8(¼ + λ²), |β(λ)| ≤ 2√2 w(λ). Let Q be the L²_loc limit of Q_H(y) = Σ_{|λ|≤H} β(λ) e^(iλy) ([PN] Lemma 2(b)); Q is real since β(−λ) = β̄(λ).

**Lemma 5 (shifted cross means; proved).** Let Λ*, Λ′* be two sets of real ordinates, each with ≪ ln(|t| + 2) points (with multiplicity) per unit interval, and β, β′ coefficients with |β| ≤ K w, |β′| ≤ K w′. Let Q, Q′ be their L²_loc limits. Then, uniformly in b ∈ ℝ and c,

```
(1/L) ∫_c^(c+L) Q(y + b) Q̄′(y) dy  →  Σ_{λ ∈ Λ*∩Λ′*} β(λ) β̄′(λ) e^(iλb)      (L → ∞).
```

*Proof.* Q(· + b) is the L²_loc limit of the series with coefficients β(λ)e^(iλb), which have the same moduli, so [PN] Lemma 2(b) applies to it on every interval. For I = [c, c + L] and k_I(θ) = (1/L) ∫_I e^(iθy) dy, (1/L) ∫_I Q_H(y + b) Q̄′_H(y) dy = Σ β(λ) β̄′(λ′) e^(iλb) k_I(λ − λ′). As H → ∞ the left side tends to (1/L) ∫_I Q(· + b) Q̄′ (L² convergence and Cauchy–Schwarz on I). On the right, |β β̄′ k_I| ≤ K² w_∪(λ) w_∪(λ′) κ(λ − λ′), with w_∪ = (m + m′)/(1 + |λ|) the weight of the union multiset, which also satisfies the counting hypothesis; so [PN] Lemma 1 for the union (its proof uses only the counting hypothesis) makes the double series absolutely convergent. The terms λ = λ′ give the claimed sum (k_I(0) = 1). The others are bounded by Σ_{λ≠λ′} |β||β′| min(1, 2/(L|λ − λ′|)), independent of b and c, which tends to 0 by dominated convergence. ∎

With Q′ the single frequency 0 (coefficient 1), Lemma 5 gives: the uniform mean of Q is β(0), which is 0 when 0 is not an ordinate.

**Lemma 6 (identification; proved).** For almost every y > 0, E(e^y) e^(−y/2) = −Q(y) + e(y) and E₀(e^y) e^(−y/2) = −Q(y), with e(y) = −e^(−y/2)[ln 2π + ½ ln(1 − e^(−2y))] ∈ L²(0, ∞).

*Proof.* [PN] Lemma 3 with δ removed: by (I1), E(x) = −lim_H Σ_{|γ|≤H} x^ρ/ρ − ln 2π − ½ ln(1 − x⁻²) for x > 1 not a prime power, and under RH Σ_{|γ|≤H} x^ρ/ρ = √x Q_H(ln x). A subsequence of Q_H converges almost everywhere to Q. Near y = 0, e(y) ≪ 1 + |ln y|; for large y, e(y) ≪ e^(−y/2). Both are square-integrable. ∎

*Proof of (a)–(c) (proved).* Let b = ln a. For y > c_a,

```
F_a(y) = √a [−Q(y + b) + e(y + b)] [−Q(y) + e(y)].
```

By Lemma 5 with Q′ = Q, (1/L) ∫_c^(c+L) √a Q(y + b) Q(y) dy → √a Σ_λ |β(λ)|² e^(iλb) = Σ_λ m(λ)² a^(½+iλ)/(¼ + λ²), uniformly in c. The cross terms are at most √a L^(−1)(∫_{I+b} Q²)^(½) ‖e‖₂ ≤ √a B L^(−½) ‖e‖₂, with B² = sup_{|I|=L≥1} (1/L) ∫_I Q² < ∞ ([PN] Lemma 2(a)), and the e·e term is at most √a ‖e‖₂²/L. This proves (b) and the existence in (a). The bound on (1/L) ∫ |F_a| follows from |F_a| ≤ (√a/2)(|Q(· + b) − e(· + b)|² + |Q − e|²). The same estimates apply to any u with u(e^y)e^(−y/2) ∈ L²(0, ∞), i.e. ∫_1^∞ u² dx/x² < ∞, which gives (c); E₀ = E + u with u = ln 2π + ½ ln(1 − x⁻²). A(1/a) = A(a)/a: substitute x ↦ ax and use uniformity in c.

*Weights.* Let φ be Riemann-integrable on [0, 1] and ε > 0. Take step functions D ≤ φ ≤ U on M equal subintervals with ∫(U − D) < ε. On each subinterval (length L/M → ∞) the uniform mean applies, so the means with weights D and U tend to A(a)∫D and A(a)∫U. The difference of the φ-mean and the D-mean is at most Σ_j (U_j − D_j)(1/L) ∫_{I_j} |F_a| ≤ B′ ∫(U − D) < B′ε, with B′ the bound of (a). Let ε → 0. ∎

*Proof of (d) and (e) (proved).* Under RH, |ρ|² = ρ(1 − ρ), so with simple zeros (b) reads A(a) = Σ_ρ a^ρ/(ρ(1 − ρ)) = J(a); with multiplicities, Σ_λ m² = Σ_λ m + Σ_λ m(m − 1). For (e): the Riesz-mean formula [PN] (4.2), proved without hypothesis, is

```
P(y) := Σ_{n≤y} Λ(n)(y/n − 1) = y ln y − (1 + γ) y + ln 2π + J(y) − T(y)          (y ≥ 1).
```

Since P(y) = y ∫_1^y ψ(t) dt/t² and y ∫_1^y t dt/t² = y ln y, this is y[1 + γ + ∫_1^y E(t) dt/t²] = ln 2π + J(y) − T(y), the head form of (1.1). Divide by y and let y → ∞: J(y)/y → 0 by [PN] C7(a) (no hypothesis) and T(y) → 0, so ∫_1^∞ E dt/t² converges to −1 − γ. Then −a ∫_a^∞ E dt/t² = a[1 + γ + ∫_1^a E dt/t²]. Finally E₀ − E = ln 2π + ½ ln(1 − t⁻²), and with u = 1/t, (a/2) ∫_a^∞ ln(1 − t⁻²) dt/t² = (a/2) ∫_0^(1/a) ln(1 − u²) du = −Σ_k a^(−2k)/(2k(2k + 1)) = −T(a). So −a ∫_a^∞ E₀ dt/t² = −a ∫_a^∞ E dt/t² − ln 2π + T(a) = J(a). The script checks the last integral by quadrature at a = 1.5, 3, 10 (to 10⁻²⁵) and the closed form of T against its series. ∎

*Remarks.* (1) The case a = 1 is Cramér's mean value theorem (Math. Z. 12 (1922)), with the constant written in the classical form Σ 1/|ρ|² = 2 + γ − ln 4π. (2) **Relation to [PN].** Polarisation gives ⟨(E((1 + δ)x) − E(x))²/x⟩ = (2 + δ)A(1) − 2A(1 + δ), using ⟨E(yx)²/x⟩ = yA(1). By [PN] Lemma 4 the left side is δ·M(δ). So Theorem 1(b) and [PN] Theorem A imply each other (A(1) is recovered from Theorem A as δ → ∞, since |c_δ(λ)|²/δ → 1/(¼ + λ²) boundedly), and (1.1) at one point is [PN] Theorem B written for one ratio instead of a difference: δ𝔇(δ) = (2 + δ)C − 2J(1 + δ). The script checks 𝔇(1) = 3C − 2J(2) = ln(27/(4π)) − γ = 0.187596954133505. What the echo form adds: the reading as an equation in E (SC), the criterion of Theorem 2, the bootstrap of Theorem 3, the imaginary part and the cross law for characters (polarisation of |Δ|² gives only Re), and §6. (3) *The tube echo (proved, RH).* For the tube function J of fact 2048, J(e^y)e^(−y/2) = Σ_λ m e^(iλy)/(¼ + λ²) converges absolutely, so it is Bohr almost periodic and ⟨J(ax)J(x)/x⟩ = Σ_λ m² a^(½+iλ)/(¼ + λ²)², which with simple zeros is J₂(a) = Σ_ρ a^ρ/(ρ(1 − ρ))² of [PN] C8 (kernel (s(1 − s))⁻²). The script checks the closed form of J₂ against the 40,000 zeros at y = 1, 2, 5, 10.5, 40 (to 2.6·10⁻¹⁴).

**Without RH.** The right side of (1.1) is defined and equal to the tail integral for every a ≥ 1 (Theorem 1(e)). The left side fails at a = 1: by Theorem 2(i), (1/ln X) ∫_1^X E(x)² dx/x² has limsup ∞. For a ≠ 1 we prove nothing (§9).

---

## 3. Proof of Theorem 2

*(i) (proved).* For Re s > 1, −ζ′/ζ(s) − s/(s − 1) = s ∫_1^∞ E(x) x^(−s−1) dx. Suppose I(X) ≤ c X^θ for X ≥ 1, with 0 < θ. By Cauchy–Schwarz on [X, 2X], ∫_X^(2X) |E| x^(−σ−1) dx ≤ I(2X)^(½) (∫_X^(2X) x^(−2σ) dx)^(½) ≪ X^((θ + 1 − 2σ)/2), which is summable over X = 2^k when σ > (1 + θ)/2. So the Mellin integral converges absolutely and locally uniformly there and is holomorphic, and so is −ζ′/ζ(s) − s/(s − 1). Hence ζ has no zero with real part > (1 + θ)/2. If a zero has real part β > ½, this is false for every θ < 2β − 1. Taking θ ∈ (0, 2β − 1), I(X) = O(ln X) would give I(X) = O(X^θ); so limsup I(X)/ln X = ∞. ∎

*(ii), (iii) (proved).* Under RH, Theorem 1 gives I(X)/ln X → A(1) = Σ_λ m²/(¼ + λ²) < ∞. Conversely limsup < ∞ excludes zeros off the line by (i), and existence of a finite limit implies limsup < ∞. For (iii): if limsup ≤ C, then RH holds by (ii), the limit exists and equals Σ m²/(¼ + λ²), while C = J(1) = Σ m/(¼ + λ²) under RH. So Σ m(m − 1)/(¼ + λ²) ≤ 0, and every m = 1. The converse is Theorem 1(d) at a = 1. ∎

*(iv) (proved).* The difference is Theorem 1(b) minus (d). D(0) = Σ m(m − 1)/(¼ + λ²) is a sum of terms ≥ 0, zero iff all m = 1. If D = 0 on [b₁, ∞): the series converges absolutely (m ≪ ln(|λ| + 2)), so (1/B) ∫_{b₁}^(b₁+B) D(b) e^(−iμb) db → the coefficient at μ, by dominated convergence term by term; all coefficients vanish. If only finitely many m > 1, D is a trigonometric polynomial, which is real-analytic, so vanishing on an open interval makes it vanish identically. ∎

*What is not covered.* An almost periodic function with coefficients ≥ 0 can vanish on a bounded interval (a periodised triangle function has Fejér coefficients ≥ 0), so (iv) gives nothing for "A = J on a bounded interval not containing a = 1" with infinitely many multiple zeros. Whether the existence of A(a) on a half-line a ≥ a₁ > 1 forces RH is also left open (§9).

**The same for characters (proved).** −L′/L(s, χ) = s ∫_1^∞ ψ(x, χ) x^(−s−1) dx (no pole), and C_χ = J_χ(1) = Σ_ρ 1/(ρ(1 − ρ)) is real without hypothesis ([PN] §5). The proof is word for word: GRH(χ) with simple zeros ⟺ limsup (1/ln X) ∫_1^X |ψ(x, χ)|² dx/x² ≤ C_χ.

---

## 4. Theorem 3: the law as an equation in E

### 4.1 The special cases

- *a = 1 (proved).* (SC) is Cramér's mean value with the classical constant: ⟨E²/x⟩ = −∫_1^∞ E dt/t² − ln 2π + T(1) = 1 + γ − ln 2π + 1 − ln 2 = C.
- *Near a = 1 (proved; the expansion from [PN] C2).* With r(a) = A(a)/(√a A(1)) the correlation coefficient of E(ax)/√(ax) and E(x)/√x, and simple zeros, r(1 + δ) = 1 − (δ/(2C))(ln(1/δ) + B) + O(δ² ln(1/δ)). The slope of J at 1⁺ is −∞, since T′(y) = 1/y − ½ ln((y + 1)/(y − 1)) (sympy). Computed: r(1.001) = 0.9406, r(1.01) = 0.6540 (the two-term expansion gives 0.9405 and 0.6547), r = ½ at a = 1.01743, and r first vanishes at a = 1.06585. The error forgets its sign within 7% of x.
- *Increments (proved).* For 1 ≤ a₁ < a₂, Theorem 1 at a₁ and a₂ gives ⟨[ψ(a₂x) − ψ(a₁x) − (a₂ − a₁)x] E(x)/x⟩ = J(a₂) − J(a₁): the excess of prime powers in (a₁x, a₂x] correlates with E(x) by an explicit prime sum. *The derivative (heuristic).* Formally, dA/da = ⟨x(ψ′(ax) − 1) E(x)/x⟩ = J′(a) = Σ_{n≤a} Λ(n)/n − ln a + γ + T′(a) = −∫_a^∞ E₀ dt/t² + E₀(a)/a: the Λ-fluctuation at ax against E at x. The limit a₂ → a₁ inside the mean is not justified here. What is proved is the corner: J′ jumps by Λ(n)/n at every prime power n (from P′(a) = Σ_{n≤a} Λ(n)/n) and is smooth elsewhere on (1, ∞). Between prime powers J″(a) = −1/a + T″(a) = −1/a + 1/(a²(a² − 1)) (sympy). This is negative for a > a* = 1.32472 (the real root of a³ = a + 1). So beyond a* the function J is concave on each piece: its local minima lie at prime powers and its local maxima inside the pieces. On (1, a*) J is convex, and r has one interior minimum there, r(1.16482) = −0.20888.
- *a → ∞ (proved, RH and simple zeros).* |A(a)| ≤ C√a. b ↦ r(e^b) = Σ_λ e^(iλb)/(C(¼ + λ²)) is Bohr almost periodic with coefficients ≥ 0 summing to 1 and no constant term (ζ(½) ≠ 0). Hence: its mean over ln a is 0; its mean square is K₂/C², RMS √K₂/C = 0.1865; limsup_{a→∞} r(a) = 1 (Dirichlet's simultaneous approximation for the first N ordinates); and liminf_{a→∞} r(a) ≤ min_{[1,1000]} r = r(3) = −0.45400, because every value of an almost periodic function recurs to within ε on a relatively dense set of shifts. Under LI, Kronecker's theorem gives liminf = −1 (conjecture-level, since LI is). So the error at x and at ax never decorrelates for good: for arbitrarily large a it is again almost perfectly correlated.

### 4.2 Proof of Theorem 3

*(a) (proved).* Lemma 5 applies to c (Λ satisfies the counting hypothesis, |c| ≤ K/(1 + |λ|)), so the left mean is √a Σ|c|² a^(iλ), as in §2. For the right side, ∫_a^∞ R dt/t² = −∫_{ln a}^∞ Q_c(y) e^(−y/2) dy converges absolutely, because the unit-interval L² norms of Q_c are bounded ([PN] Lemma 2(a)). Termwise integration is allowed: |∫_{ln a}^∞ (Q_c − Q_{c,H}) e^(−y/2) dy| ≤ Σ_k e^(−k/2) sup_{|I|=1} ‖Q_c − Q_{c,H}‖_{L²(I)} → 0 by [PN] Lemma 2(b). With ∫_{ln a}^∞ e^((iλ−½)y) dy = a^(iλ−½)/(½ − iλ), −a ∫_a^∞ R dt/t² = √a Σ c(λ) a^(iλ)/(½ − iλ), absolutely convergent. If c = 1/(½ + iλ), both sides are √a Σ a^(iλ)/(¼ + λ²). Conversely, if they agree for a ≥ a₁, the absolutely convergent series f(b) = Σ (|c(λ)|² − c(λ)/(½ − iλ)) e^(iλb) vanishes on [ln a₁, ∞). One-sided Bohr means give |c(μ)|² = c(μ)/(½ − iμ) for every μ ∈ Λ. So c(μ) = |c(μ)|²(½ − iμ), whose modulus gives |c(μ)| = 1/|½ − iμ|, hence c(μ) = 1/(½ + iμ). ∎ (At a = 1 alone this fails in the class: Σ(|c|² − c/(½ − iλ)) = 0 is one real equation. For E itself a = 1 suffices because the coefficients m/ρ have integer m ≥ 1.)

*(c) (proved).* By (1.1), (A(a) + ln 2π − T(a))/a = 1 + γ + ∫_1^a E dt/t² on [1, a₀]. The right side is absolutely continuous with derivative E(a)/a² almost everywhere, so A on [1, a₀] determines E almost everywhere on [1, a₀], hence ψ (right-continuous) everywhere there. Conversely A on [1, a₀] is a function of ψ on [1, a₀]. ∎

*Remark (heuristic).* In the zero picture, A(e^b)e^(−b/2) on 0 ≤ b ≤ ln a₀ is the Fourier transform of the measure Σ δ_λ/(¼ + λ²) seen through a window of length ln a₀. By the uncertainty principle this resolves ordinates only to about 2π/ln a₀: the classical duality "primes up to a₀ ↔ zeros at resolution 1/ln a₀", here in the two-point function. Not made precise.

---

## 5. Characters (Theorem 4)

*Setting.* χ is primitive mod q > 1, χ(−1) = (−1)^κ. By (I1) for χ, ψ₀(x, χ) = −lim_H Σ x^ρ/ρ + A_χ(x) with A_χ(x) = −(1 − κ) ln x − b(χ) + Σ_{k≥1} x^(κ−2k)/(2k − κ); b(χ) = L′/L(0, χ) for odd χ and B₀ for even χ. For odd χ the series is artanh(1/x); for even χ it is −½ ln(1 − x⁻²). Put E_χ,₀ = ψ(x, χ) − A_χ(x).

*Proof of Theorem 4 (proved).* As §2: e_χ(y) = e^(−y/2) A_χ(e^y) ∈ L²(0, ∞) (O(y e^(−y/2)) at infinity, logarithmic at 0). Lemma 5 needs no symmetry of the ordinates, so for complex χ the mean of E_χ(ax)Ē_χ(x)/x is √a Σ_λ |β_χ(λ)|² e^(iλb). Under GRH, 1 − ρ = ρ̄, so with simple zeros it is J_χ(a). The closed forms (1.2) are [PN] §5 (J_χ(y) = P_χ(y) + y L′/L(1, χ) − L′/L(0, χ) + T_o(y) for odd χ; −ln y − 1 − B₀ + T(y) in place of the last two terms for even χ), with P_χ(y) = y ∫_1^y ψ(t, χ) dt/t². The tail form follows as in §2. ∫_1^∞ ψ(t, χ) dt/t² = −L′/L(1, χ) without hypothesis, because J_χ(y)/y → 0 (zeros have Re ρ < 1, and Σ 1/|ρ(1 − ρ)| < ∞). The tail of A_χ: a ∫_a^∞ artanh(1/t) dt/t² = T_o(a) and a ∫_a^∞ ln t dt/t² = ln a + 1 (checked by quadrature). The cross law is Lemma 5 with Λ* ≠ Λ′*; at a = 1 the sum Σ m m′/(¼ + λ²) has terms > 0, and for all a its vanishing forces every coefficient to vanish by Bohr uniqueness. ∎

*Remarks.* (1) "No common zero" (under GRH, no common ordinate) is the two-function case of LI. It is much weaker than LI, which also gives simple zeros and no additive relations. (2) With χ′ = χ̄, ⟨E_χ(ax)E_χ(x)/x⟩ (no conjugate) is the sum over ordinates λ with both λ and −λ zeros of L(s, χ); for complex χ it vanishes iff L(s, χ) and L(s, χ̄) share no zero. (3) The constants: L′/L(1, χ₋₄) = 0.2456095848, L′/L(0, χ₋₄) = 0.7831887854, L′/L(1, χ₋₃) = 0.3682816160, L′/L(0, χ₋₃) = 0.9481988267 (mpmath by Hurwitz–Stieltjes and Lerch's formula, equal to PARI's to 10⁻¹⁴, and satisfying L′/L(0) + L′/L(1) = −ln(q/π) + γ + ln 2 to 10⁻³⁰). C_χ₋₄ = 0.155568, C_χ₋₃ = 0.113230. For the even character (·/5): B₀ = −0.0220247, C = 0.156557 (fact 2049).

**The race (proved; GRH for L(s, χ₋₄), simple zeros for the closed form).** Let χ = χ₋₄ and N(x) = (ln x/√x)(π(x; 4, 3) − π(x; 4, 1)).

1. θ(x; 4, 3) − θ(x; 4, 1) = −θ(x, χ) = −ψ(x, χ) + Σ_{p^k≤x, k≥2} χ(p)^k ln p = −ψ(x, χ) + θ(√x) − ln 2 + O(x^(1/3)), since χ(p)² = 1 for odd p. θ(√x)/√x → 1 by the prime number theorem, with no RH for ζ needed. So R_θ(x) := (θ(x; 4, 3) − θ(x; 4, 1))/√x = 1 − ψ(x, χ)/√x + o(1). The prime squares give the 1.
2. By partial summation, N(x) = R_θ(x) − (ln x/√x) ∫_2^x θ(t, χ) dt/(t ln² t). Under GRH, ∫_2^t ψ(u, χ) du = O(t^(3/2)) (integrated explicit formula), and one more integration by parts makes the integral O(√x/ln² x). So N(x) = 1 − ψ(x, χ)/√x + O(1/ln x) + o(1).
3. The o(1) terms do not change log-means of products (Cauchy–Schwarz against the bounded mean square). L(½, χ₋₄) = 0.667691 ≠ 0, so 0 is not an ordinate and ⟨ψ(x, χ)/√x⟩ = 0 (Lemma 5 with the constant). Hence ⟨N⟩ = 1, ⟨N(ax)N(x)⟩ = 1 + J_χ(a)/√a, the covariance is J_χ(a)/√a, the variance C_χ, and the correlation coefficient is r_χ(a) = J_χ(a)/(√a C_χ).

The prime squares shift the mean to 1 and do not enter the covariance. Values: covariance J_χ(a)/√a = −0.0214575, 0.071325, −0.083571 at a = 2, 3, 5; correlation coefficients:

| a | 2 | 3 | 5 |
|---|---|---|---|
| r(a) for ζ: E(ax) vs E(x) | −0.3752 | −0.4540 | −0.4533 |
| r_χ(a) for χ₋₄: the race at ax vs x | −0.1379 | +0.4585 | −0.5372 |

At a = 3 the two differ in sign because χ₋₄(3) = −1 turns the corner at 3 upside down (§7.3).

---

## 6. Higher order, and what does not close

### 6.1 Three points (heuristic for E; proved for the tube)

Formally, under RH, with β = m/(½ + iλ),

```
⟨E₀(ax) E₀(bx) E₀(x)/x^(3/2)⟩ = −√(ab) Σ_{λ₁+λ₂+λ₃=0} β(λ₁)β(λ₂)β(λ₃) a^(iλ₁) b^(iλ₂).
```

Nonzero terms need an additive relation ±γ₁ ± γ₂ ± γ₃ = 0 among the ordinates (2γ = γ′ included). Under LI there is none, and the three-point mean is 0. Without LI, the sum over relations is not a residue sum of any generating function we know, so it does not close. For E the interchange of the mean with the triple sum needs a bound on the third moment of Q that is not proved here (the square-integrability of Lemma 5 is not enough). For the tube function J (absolutely convergent series) the same formula with +√(ab) in place of −√(ab) and β = m/(¼ + λ²) is a theorem under RH, by Bohr's product rule for almost periodic functions.

### 6.2 Four points under LI (conjecture for E; theorem for the tube)

Assume RH and LI: the positive ordinates, counted with multiplicity, are linearly independent over ℚ (this includes simple zeros). Put P(c) = e^(−c/2) J(e^c) and K(c) = e^(−c/2) J₂(e^c), both even in c, with J₂ the closed form of [PN] C8. Under LI the only solutions of λ₁ + λ₂ + λ₃ + λ₄ = 0 are the three pairings; tuples that satisfy two pairings at once have the form (μ, −μ, −μ, μ) and its permutations. Inclusion–exclusion gives, for Q = −E₀(e^y)e^(−y/2) and shifts α_j,

```
⟨Π_j Q(y + α_j)⟩ = P(α₁−α₂)P(α₃−α₄) + P(α₁−α₃)P(α₂−α₄) + P(α₁−α₄)P(α₂−α₃)
                    − K(α₁−α₂−α₃+α₄) − K(α₁−α₂+α₃−α₄) − K(α₁+α₂−α₃−α₄).
```

The coincidence term Σ_μ |β_μ|⁴ e^(iμc) = Σ_ρ e^(ρc)e^(−c/2)/(ρ(1 − ρ))² is again a residue sum, with kernel (s(1 − s))⁻², so the four-point function also closes. At α = 0: ⟨E₀⁴/x²⟩ = 3C² − 3K₂ = 0.0061783 and the kurtosis is 3 − 3K₂/C² = 2.8957. The fourth cumulant is −3K₂, so under LI the limiting distribution of E(x)/√x has excess kurtosis −3K₂/C² = −0.1043 against a Gaussian. In general, under LI every joint cumulant of order 2r of (E(a_j x)/√(a_j x))_j is a linear combination of J_r(Π a_j^(±1)) = Σ_ρ (·)^ρ/(ρ(1 − ρ))^r, with kernels (s(1 − s))^(−r); odd cumulants vanish. For E this is a conjecture: it needs LI and the moment interchange of §6.1. For the tube function it is a theorem under RH and LI. §7.4 checks the four-point formula at α = 0 on primes.

### 6.3 Why it closes for Λ, and for what else

Under RH the frequency λ of E₀(x)/√x carries the amplitude −m/ρ, and m/ρ = Res_{s=ρ}[(ζ′/ζ)(s)/s]: an integer (the residue of a logarithmic derivative) times the elementary kernel 1/s. The two-point mean keeps only the diagonal, |m/ρ|² = m²·(1/ρ)(1/ρ̄). On the critical line ρ̄ = 1 − ρ, so the conjugate kernel is the reflected kernel 1/(1 − s). When m = 1, m² = m, so the squared amplitude is again a residue of the same logarithmic derivative: 1/(ρ(1 − ρ)) = Res_{s=ρ}[(ζ′/ζ)(s)/(s(1 − s))]. The kernel 1/(s(1 − s)) is elementary, so the contour can be moved to Re s = 2, where −ζ′/ζ is the Dirichlet series of Λ. The sum of residues becomes the Riesz sum Σ_{n≤a} Λ(n)(a/n − 1) plus the residues at s = 1, s = 0 and the trivial zeros: a linear functional of ψ. Nothing in this uses the values of ζ at its zeros, only the integers m. The mechanism needs three things: (M1) the frequencies are poles of a Dirichlet series whose residues are all equal (a logarithmic derivative −L′/L with simple zeros); (M2) the poles lie on a line of symmetry, so that conjugating the kernel is reflecting it (GRH and the functional equation); (M3) finite order and a functional equation, so that the contour can be moved and leaves only elementary residues. It therefore works for Λ and for Λχ with χ primitive (Theorems 1, 4). The same proof, given the explicit formula with its gamma-factor residues (not written out here), covers Λ_K when ζ_K has simple zeros and the von Mangoldt coefficients of any L-function of the Selberg class with simple zeros on the critical line. It fails for μ and for Liouville's λ: their residues at ρ are 1/ζ′(ρ) and ζ(2ρ)/ζ′(ρ), values rather than integers, so squaring creates new transcendental weights. It fails for d(n) and r₂(n), whose oscillation does not come from poles in the strip at all. And it fails for Λ itself when a zero is multiple, since m² ≠ m.

*Mertens (heuristic).* Under RH, simple zeros and Ng's hypothesis Σ_{0<γ≤H} |ζ′(ρ)|⁻² ≪ H (N. Ng, 2004, quoted from memory), M(x)/√x has a limiting distribution. The same diagonal argument would give ⟨M(ax)M(x)/x⟩ = Σ_ρ a^ρ/|ρζ′(ρ)|². The analogue of Lemma 1 for these weights is not checked here. On the line, |ζ′(ρ)|² = ζ′(ρ)ζ′(1 − ρ), so each term is Res_{s=ρ}[a^s/(s(1 − s)ζ(s)ζ′(1 − s))]. But this kernel has poles at 1 − ρ′ for every zero ρ′ of ζ′, and it is not a Dirichlet series times an elementary factor in any right half-plane: by the functional equation, ζ′(1 − s) = −g(s)[ζ′(s) + h(s)ζ(s)] with g(s) = 2(2π)^(−s) cos(πs/2) Γ(s) and h = g′/g, which mixes arithmetic and archimedean terms. Moving the contour produces sums over the zeros of ζ′ instead of a finite prime sum. The two-point function of M is not a linear functional of M.

*Divisor problem (heuristic).* For Δ(x) = Σ_{n≤x} d(n) − x ln x − (2γ − 1)x, ζ(s)² has no poles in the strip. The oscillation comes from the functional equation (Voronoi's series, frequencies 4π√n in √x, amplitudes d(n)n^(−3/4)). The diagonal of Δ(ax)Δ(x) pairs n with na, so the two-point mean vanishes unless a is rational. At a = 1 it is a Rankin–Selberg constant built on Σ d(n)² n^(−3/2) = ζ(3/2)⁴/ζ(3) (Cramér, Tong; quoted from memory), a different arithmetic function, d², and not a linear functional of Δ.

*Number fields with forced multiple zeros (follows by the proof of Theorem 4 under the stated hypotheses, given the explicit formula for Artin L-functions; not written out, not checked numerically).* If K is Galois over ℚ with a non-abelian group, then ζ_K = Π_π L(s, π)^(dim π) (Artin), and every zero of L(s, π) with dim π ≥ 2 is a multiple zero of ζ_K. So the simple-zero form of the law fails for ψ_K(x) − x. Assume GRH, Artin's conjecture for these π, simple zeros of each L(s, π) and no common zeros. Since Λ_K = Σ_π (dim π)Λ_π, Theorem 4's cross law gives ⟨E_K(ax)E_K(x)/x⟩ = Σ_π (dim π)² J_π(a), not J_K(a) = Σ_π (dim π) J_π(a). Example: the splitting field of x³ − 2, where ζ_K = ζ·L(s, χ₋₃)·L(s, ρ₂)² and L(s, ρ₂) is a Hecke L-function of ℚ(√−3), hence entire.

---

## 7. Numerics (computed)

### 7.1 J(a) from 40,000 zeros of ζ, 1 ≤ a ≤ 50

Zero side: 2√a[Σ_{n≤N} cos(γ_n u)/(¼ + γ_n²) + ∫_H^∞ cos(tu) θ′(t)/(π(¼ + t²)) dt − cos(uH)S(H⁺)/(¼ + H²)], with u = ln a, H = γ_N = 33190.012 and S(H⁺) = N − 1 − θ(H)/π = 0.4648. The tail integral is done by a 10-term integration by parts (u ≥ ln 1.01) or by quadrature (u = 0). Grid: a = 1 and a = 1.01, 1.02, …, 50 (4,901 values).

- At a = 1: zeros − closed = −3.9·10⁻¹⁴.
- Away from prime powers (|ln a − ln n| > 0.004 for every n): max |zeros − closed| = 2.5·10⁻⁷.
- *Resonance at prime powers* (heuristic, as [PN] C9): with the explicit series of S(t), the truncated sum exceeds the closed form by (Λ(n)/π)√(a/n)·I(H), I(H) = cos(wH)/H − |w|(π/2 − Si(|w|H)), w = ln a − ln n. At a = n this is Λ(n)/(πH) (Landau's formula seen in J). At each prime power n ≤ 50 the observed difference equals Λ(n)/(πH) to within 0.6% (Λ(2)/(πH) = 6.648·10⁻⁶, observed 6.648·10⁻⁶ at a = 2). After subtracting the resonance terms of all prime powers ≤ 200, max |zeros − closed − resonance| = 5.0·10⁻⁹ over the whole grid. Example: at a = 6, zeros − closed = −3.6·10⁻⁹, residual 5.7·10⁻¹⁰.
- Closed values: J(1.5) = 0.0065637, J(2) = −0.0245114, J(2.5) = 0.0157764, J(6) = 0.0328967, J(10) = 0.0245927, J(30.5) = −0.0614892, J(50) = 0.0417043.

### 7.2 Correlation coefficients r(a) = A(a)/(√a A(1)) (closed form)

| a | 1.01 | 1.1 | 1.5 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 100 | 1000 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ζ | +0.6540 | −0.1369 | +0.1160 | −0.3752 | −0.4540 | −0.0457 | −0.4533 | +0.2907 | −0.4418 | −0.0445 | +0.0285 | +0.1684 | −0.2515 | −0.2101 |
| χ₋₄ | +0.8534 | +0.2822 | −0.2671 | −0.1379 | +0.4585 | −0.1305 | −0.5372 | −0.0066 | +0.4536 | +0.2296 | +0.0371 | +0.1161 | −0.4695 | −0.2642 |
| χ₋₃ | +0.8112 | +0.1345 | −0.2300 | +0.3122 | −0.1021 | −0.1850 | +0.5034 | −0.0591 | −0.5185 | −0.0414 | +0.1301 | +0.2879 | +0.1209 | +0.0258 |

Extremes on [1, 1000] (grid of 2,000,001 points uniform in ln a, refined at the prime powers and by bounded minimisation):

| | first zero | minimum | largest value after the first zero | RMS over ln a on [1, 1000] | asymptotic RMS √K₂/C |
|---|---|---|---|---|---|
| ζ | 1.06585 | −0.45400 at a = 3 | 0.39341 at a = 208.115 | 0.1854 | 0.1865 |
| χ₋₄ | 1.19747 | −0.53720 at a = 5 | 0.60762 at a = 491 | 0.2730 | 0.2754 |
| χ₋₃ | 1.13650 | −0.51852 at a = 7 | 0.51542 at a = 257 | 0.2304 | 0.2314 |

The maximum is r(1) = 1. For ζ the minimum sits on a corner (prime powers are upward kinks), the large maxima inside pieces, and the first dip, r(1.16482) = −0.20888, inside the convex piece (1, 1.32472). For χ the corners at n with χ(n) = −1 are downward kinks, so maxima can sit on them: 491 ≡ 3 (mod 4), 257 ≡ 2 (mod 3). The asymptotic RMS for χ uses K₂,χ = Σ 1/|ρ|⁴ from the PARI zeros to height 1000 (1.835740·10⁻³ for χ₋₄, 6.866039·10⁻⁴ for χ₋₃).

### 7.3 Primes to 10⁸: the law, its corners, the cross law, the race

**Setup.** Sieve to 10⁸ (π(10⁸) = 5,761,455, ψ(10⁸) − 10⁸ = −1757.2034). The test ∫_1^X E dt/t² + 1 + γ at X = 10⁸ gives −1.1·10⁻⁶, inside the bound (C√X + ln 2π + 1)/X = 4.6·10⁻⁶. Mean over y = ln x ∈ [y₀, y₀ + L] = [ln 10², ln 10⁷] (L = 11.5129) with the Hann weight w(y) = (1 − cos(2π(y − y₀)/L))/L, of the oscillating parts E₀(ax)E₀(x)/x and E_χ,₀(ax)E_χ,₀(x)/x. The step functions are integrated exactly between consecutive jumps of ψ(x) and ψ(ax), with 3-point Gauss–Legendre for the smooth factors; ax ≤ 10⁸ for a ≤ 10. *Finite window (identity proved under RH and simple zeros as in Lemma 5; sum computed).* By Lemma 6 the window mean equals √a Σ_{λ,λ′} β(λ)β̄(λ′) a^(iλ) K(λ − λ′) exactly, with K(θ) = ∫w(y)e^(iθy)dy = e^(iθy_c)(sin z/z)·π²/(π² − z²), z = θL/2, y_c = y₀ + L/2. Its diagonal is the law; the off-diagonal ("window term") is computed from the zeros. For ζ we use all 40,000 zeros (pairs with |λ − λ′| < 60; widening the band to 120 changes the term by less than 10⁻⁶ at a = 1, 2, 10); for the characters, the PARI zeros to height 1000 (868 for χ₋₄, 823 for χ₋₃). Their completeness was checked by the mean of S(γ_n⁺) = n − θ_χ(γ_n)/π in each block of height 100 (0.50 to within 0.005 in every block) and by C_χ from the zeros, which matches the closed form to 2·10⁻⁸.

| a | ζ: primes | ζ: law J(a) | χ₋₄: primes | χ₋₄: law | χ₋₃: primes | χ₋₃: law |
|---|---|---|---|---|---|---|
| 1 | +0.046168 | +0.046191 | +0.155546 | +0.155568 | +0.113120 | +0.113230 |
| 2 | −0.024498 | −0.024511 | −0.030498 | −0.030346 | +0.049976 | +0.049988 |
| 3 | −0.036279 | −0.036323 | +0.123648 | +0.123539 | −0.019783 | −0.020029 |
| 10 | +0.024159 | +0.024593 | +0.057319 | +0.057107 | +0.102915 | +0.103070 |

Over 11 ratios (1, 1.5, 2, 2.9, 3, 3.1, 4.9, 5, 5.1, 7, 10) and the three functions:

- max |primes − law| = 7.6·10⁻⁴ (the window term);
- ζ: max |primes − law − window| = 3.8·10⁻⁶ with 40,000 zeros. The window term itself changes by up to 2.8·10⁻⁵ between zeros to height 1000 and all 40,000, which measures the truncation;
- χ₋₄ and χ₋₃ (zeros to height 1000): max |primes − law − window| = 4.8·10⁻⁵ and 4.7·10⁻⁵, the size of the truncation effect seen for ζ.

**Corners.** Second differences A(n + 0.1) − 2A(n) + A(n − 0.1):

| | primes | law | law without its corner (minus 0.1·χ(n)Λ(n)/n) |
|---|---|---|---|
| ζ at 3 | +0.03340 | +0.03343 | −0.00319 |
| ζ at 5 | +0.02956 | +0.03021 | −0.00198 |
| χ₋₄ at 3 (χ(3) = −1) | −0.03635 | −0.03620 | +0.00042 |
| χ₋₄ at 5 (χ(5) = +1) | +0.03054 | +0.03227 | +0.00008 |
| χ₋₃ at 3 (χ(3) = 0) | −0.00002 | +0.00042 | +0.00042 |
| χ₋₃ at 5 (χ(5) = −1) | −0.03131 | −0.03211 | +0.00008 |

Primes and law + window term agree to within 10⁻⁴ in all six cases. The corners are there with the sign of χ(n): at a = 3, ζ has an upward kink, χ₋₄ a downward one and χ₋₃ none.

**Cross law.** ⟨E₀(ax)E_χ₋₄,₀(x)/x⟩ and ⟨E_χ₋₄,₀(ax)E_χ₋₃,₀(x)/x⟩ in the same window. The law says 0 in the limit: no common ordinate below 1000, with closest pairs 5.2·10⁻⁴ (ζ, χ₋₄) and 1.6·10⁻⁴ (χ₋₄, χ₋₃) apart. The window value is the cross double sum over zeros to height 1000:

| a | ζ × χ₋₄: primes | window prediction | χ₋₄ × χ₋₃: primes | window prediction |
|---|---|---|---|---|
| 1 | −0.001063 | −0.001029 | +0.005724 | +0.005705 |
| 2 | −0.007410 | −0.007409 | +0.005138 | +0.005096 |
| 3 | +0.003865 | +0.003892 | +0.005583 | +0.005543 |

The cross means are of the order of the window term, at most 0.0074, and are explained by it to within 4.3·10⁻⁵ (the truncation error). In this window the cross correlation is fully accounted for by pairs of nearby zeros of the two L-functions; nothing remains that would indicate a common zero.

**Race** (same window). Window covariance Cov_w(N(ax), N(x)) of the normalised race N, and of the θ-race R_θ:

| a | 1 | 2 | 3 | 5 |
|---|---|---|---|---|
| N (π-race) | 0.15696 | −0.02010 | 0.07114 | −0.08349 |
| R_θ (θ-race) | 0.15997 | −0.01761 | 0.07317 | −0.08168 |
| law J_χ(a)/√a | 0.15557 | −0.02146 | 0.07133 | −0.08357 |

Window means: N 1.2233, R_θ 0.9214 (limit 1 for both). The race covariances agree with the law to within 0.0044. This is much worse than the 5·10⁻⁵ of ψ(x, χ), as §5 predicts. The prime squares give θ(√x)/√x, which tends to 1 only with fluctuations of order x^(−1/4) (θ(10)/10 − 1 = −0.4653 at x = 10², where the window starts), and π differs from θ/ln x by O(1/ln x). Both vanish in the limit and are not small in a window that starts at 10².

### 7.4 Third and fourth moments (a = 1)

Window moments of E₀/√x on the primes, against the same moments of the zero series Q on a grid (60,001 points, 120,001 for zeros to 4000):

- ζ, third moment: primes −1.689·10⁻⁴. The window value from zeros to height 1000, 2000, 4000 is −1.751, −1.703, −1.668 (·10⁻⁴), converging toward the primes. The LI limit is 0; in this window the near-relations γ₁ + γ₂ ≈ γ₃ of the low zeros dominate.
- ζ, fourth moment: primes 0.0059981. LI formula 3C² − 3K₂ plus the window part from zeros to 4000 (3,474 zeros; window part −1.808·10⁻⁴): 0.0059975 (0.0059949 with zeros to 2000). Without the −3K₂ coincidence term: 0.0062201, off by 2.2·10⁻⁴. Window kurtosis 2.8140; the LI limit is 2.8957.
- χ₋₄: third moment −9.452·10⁻⁴ (zeros to 1000: −8.870·10⁻⁴); fourth moment 0.070489 against 3C_χ² − 3K₂,χ + window part = 0.070487; window kurtosis 2.9134, LI limit 2.7724.

These check the combinatorics of §6.2 (pairings minus coincidences) on real data. They do not test LI: in a window of length 11.5 the near-relations among the low zeros contribute as much as the coincidence term.

---

## 8. Prior art, and what is new here

Classical and used as such: Cramér's mean value theorem (the case a = 1, Math. Z. 12 (1922); exact form not rechecked); the almost periodic (B²) structure of E(e^y)e^(−y/2) and the limiting distributions of Wintner and Rubinstein–Sarnak, whose covariance at shift ln a is Σ|β|²a^(iλ), i.e. Theorem 1(b); the Riesz-mean explicit formula and Σ 1/|ρ|² = 2 + γ − ln 4π; the Mellin–Cauchy–Schwarz argument of Theorem 2(i) (Landau-type); the LI cumulant structure of §6.2 (Rubinstein–Sarnak). The combinations stated here are: the reading of Theorem 1(b) through (1.1) as the closed equation (SC); the one-sided criterion "RH ∧ simple ⟺ limsup ≤ 2 + γ − ln 4π"; the bootstrap Theorem 3; the complex and cross laws with closed forms; the mechanism of §6.3. A search of the literature (24 November 2026; about 25 papers opened) gives: *known in substance*. The closed form through the primes is in Gun, Murty and Rath, "Transcendental sums related to the zeros of zeta functions" (Mathematika 64 (2018) 875–897, §3: Σ_ρ x^ρ/(ρ(1 − ρ)) in terms of ψ₀(x) and Σ Λ(n)/n, by the same derivation; they use it for transcendence questions). The tail form is in Ramaré (Acta Arith. 159 (2013), Lemma 2.2, signs corrected in the corrigendum of 2023), with the two-sided remainder of Kaczorowski–Ramaré (2003). The covariance function of the normalised remainder as a sum over zeros is in Iksanov, Kabluchko and Marynych (arXiv:2502.04969, 2025), whose plotted correlation function shows the corners at ln n without mentioning them. Not found in those sources: the two results stated as one law for the autocorrelation of ψ(x) − x, the reading of its corners as the prime powers, the criterion, the bootstrap, the character and cross laws, and a check on prime data. Suzuki's screw function (arXiv:2206.03682), with weights 1/γ² and kinks at ln n, is the nearest published analogue. Not opened: Ingham 1932, Cramér 1920, Goldston–Montgomery 1987. The repo's own [PN] already contains the equivalent variance form (Remark 2 of §2).

---

## 9. Open points

1. *Bounded intervals.* Does A(a) = J(a) for all a in a bounded interval [a₁, a₂] ⊂ (1, ∞) force simple zeros? Proved only when finitely many zeros are multiple (Theorem 2(iv)).
2. *Without RH, a ≠ 1.* If RH fails, does the mean A(a) fail to exist for some (or every) a ≠ 1? Only a = 1 is proved (Theorem 2(i)); the liminf of I(X)/ln X is not determined either.
3. *Higher moments of E.* The three- and four-point formulas need the interchange of the mean with the triple and quadruple zero sums; a third or fourth moment bound for Q under RH is not proved here. The LI limits (third moment 0, kurtosis 2.8957) are conjectural.
4. *Mertens.* The analogue of Lemma 1 for the weights 1/|ρζ′(ρ)| under Ng's hypothesis was not checked; §6.3 is heuristic there.
5. *Derivative form.* dA/da as the mean of the Λ-fluctuation against E (§4.1) is formal; only the increment form is proved.
6. *Numerics.* The resonance term of §7.1 is heuristic (the series for S(t) does not converge). The character residuals 4.8·10⁻⁵ are attributed to truncating the zeros at height 1000, by analogy with ζ (where 40,000 zeros bring the residual to 3.8·10⁻⁶); this is not proved. The completeness of the PARI zero lists rests on the Riemann–von Mangoldt count, not on a Turing-method proof.
7. *Citations quoted from memory and not rechecked:* Cramér 1922 (exact form), Ng 2004 (hypothesis), Cramér/Tong on the divisor problem, Rubinstein–Sarnak.
