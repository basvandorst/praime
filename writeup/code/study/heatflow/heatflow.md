# The zeros under the heat flow: a gas of primes with a pair interaction

> **In short:** The de Bruijn–Newman heat flow H_t multiplies the n-th Dirichlet term of ζ by exp((t/4) ln² n): the primes become a gas with pair interaction (t/2)·ln p·ln q, free at t = 0. (1) Its cluster weights κ_t(n) are positive at every integer for t > 0, and the zeros of H_t show lines at composite numbers with exactly these weights (139 composites, agreement to 10⁻⁵). (2) For Haar unitary matrices under the same flow the second moment of every power sum up to the matrix size is an exact sum over partitions (proved), with limit law K_τ(u). (3) At heights 10¹²–10²² the flowed zeros of ζ and flowed matrices have the same spacing variance and form factor at every τ tested.

**Lab log:** [LOG.md](LOG.md) (loops 1–8, each with its prediction written before the run, and the failed ones). **Proof:** [theory/proof.md](theory/proof.md), checked by [theory/check.py](theory/check.py). **Date:** 27 November 2026.

**Status:** computed, with one part proved. Section 4 (matrices: the exact formula and its limit) is proved in [theory/proof.md](theory/proof.md) and refereed once with a second code ([theory/referee/](theory/referee/): no failing step in the main theorems; one side remark of the note, that for k > N the full Hall norm is strictly larger than the true moment, is wrong for t < 0, with counterexamples at N = 1, k = 2, t = −0.346574 and N = 2, k = 3; for t ≥ 0 it holds). Sections 3, 5 and 6 are numerical: the line formula for H_t is a heuristic consequence of an effective formula in print, confirmed to 10⁻⁵; the common family of ζ and matrices is a conjecture of ours, supported by loops 6 and 8. Nothing here proves or approaches a proof of the Riemann hypothesis (section 7).

## 1. The question

H_t(z) = ∫₀^∞ e^(tu²) Φ(u) cos(zu) du, with H_0(z) = ξ(½ + iz/2)/8. The zeros of H_t are real exactly for t ≥ Λ (de Bruijn, Newman); Λ ≥ 0 (Rodgers–Tao 2018) and the Riemann hypothesis is Λ ≤ 0. The zeros move by dz_k/dt = 2 Σ_{j≠k} 1/(z_k − z_j) (Csordas–Smith–Varga). In ζ's heights γ = z/2 this is

  dγ_k/dt = ½ Σ_{j≠k} 1/(γ_k − γ_j).                                  (1)

The repo's zero-gas study ([../zerogas/LOG.md](../zerogas/LOG.md)) reads the zeros of ζ as a gas whose slow modes are the waves of the prime powers. Question: what does the flow (1) do to these modes, what are the statistics of the zeros of H_t on the scale where the gas turns into a crystal, and does ζ go through the same states as a random unitary matrix under the same flow?

## 2. In print and ours

After a literature search by one agent (Rodgers–Tao, Polymath 15, Dobner, Kabluchko, Hall–Ho, Hall–Ho–Jalowy–Kabluchko, Ki–Kim–Lee, Csordas–Smith–Varga, Newman 1976 read in full; Zenodo and GitHub searched; **not opened**: the Polymath wiki and blog threads and the original papers of Landau, Gonek and Diaconis–Shahshahani):

| item | in print | ours |
|---|---|---|
| H_t ≈ Dirichlet series with coefficients b_n = exp((t/4) ln² n), real part of the exponent shifted by (t/4) ln(T/2π) | Polymath 15 (arXiv 1904.12438, eqs. 14–15); Dobner (arXiv 2005.05142) | — |
| damping factor of the n-th term, exp(−(t/4) ln n (L − ln n)) | Rodgers–Tao (arXiv 1801.05914, p. 17, factor K_{t,n}) | the same rate from the flow (1) around a lattice (section 3.1) |
| b_n is not multiplicative for t ≠ 0 | Polymath 15, section 8.1 (noted, then ignored for an Euler-product mollifier) | the cluster weights κ_t(n) as a law, and the lines of the zeros at composite numbers, measured (section 3) |
| the flow on polynomials as a Curie–Weiss model, Lee–Yang | Newman 1976; Kabluchko (arXiv 2203.05533) | the reading for primes: pair interaction (t/2) ln p ln q |
| relaxation time ~ 1/log² T; "liquid" below exp(C/t), "solid" above | Rodgers–Tao p. 27 (heuristic); Polymath 15, Thm 1.5 and Remark 1.6; Ki–Kim–Lee (gaps tend to 1 for fixed t > 0) | statistics on that scale: the exact second moment for matrices, the limit law K_τ(u), the measurements on ζ (sections 4–5): not found |
| share of zeros off the line for t < 0 | not searched separately | the leading law (2/(9π))(2\|τ\|)^(3/2), measured on matrices (section 6) |

So the status of "ours" is: not found in the sources read, with the unopened ones named. The first term of K_τ(u) is Rodgers–Tao's factor squared; each ingredient is elementary.

## 3. The primes as a gas with interaction

**3.1 The rate of a mode (derived; agrees with the factor in print).** Put the zeros near a lattice of density ρ = L/2π, L = ln(T/2π), with a small displacement δ_k = ε cos(qγ_k). Linearising (1): dδ_k/dt = −½ Σ_{m≠0} (δ_k − δ_{k+m})/(m/ρ)² and Σ_{m≠0} (1 − cos(qm/ρ)) ρ²/m² = πρq − q²/2 for 0 ≤ q ≤ 2πρ. So the mode decays at the rate q(L − q)/4. With q = ln n = uL:

  amplitude of the wave of n at time t = e^(−τ u(1−u)),  τ = t L²/4.            (2)

This is exactly |b_n n^(−s*)| / n^(−1/2) = exp((t/4) ln² n − (t/4) L ln n) of the effective formula. The rate u(1 − u) is that of a log-gas around a lattice; for a matrix of size N it is k(N − k)/N².

In unfolded units x (mean spacing 1) and in τ, ζ and matrices obey the same equation (derived from (1) and from dθ_k/dt = Σ cot((θ_k − θ_j)/2), t = τ/N²):

  dx_k/dτ = (1/(2π²)) Σ_{j≠k} 1/(x_k − x_j).                                  (3)

**3.2 The interaction (elementary).** With ln n = Σ_i l_i over the prime factors of n (with multiplicity), b_n = ∏_i e^((t/4) l_i²) · ∏_{i<j} e^((t/2) l_i l_j). So Σ b_n n^(−s) is the grand partition function of the "primon gas" with an added pair energy −(t/2) ln p ln q (in units where s is the inverse temperature); at t = 0 it is the Euler product. The logarithm of such a sum is a sum over connected graphs (Mayer):

  log Σ_n b_n n^(−s) = Σ_n κ_t(n) n^(−s),  κ_t(pq) = e^((t/4)(l_p + l_q)²) − e^((t/4)(l_p² + l_q²)),  κ_t(p²) = e^(t l_p²) − ½ e^(t l_p²/2),

and in general κ_t(n) is computed by the recursion κ_n ln n = b_n ln n − Σ_{d | n, 1 < d < n} κ_d ln d · b_{n/d}. For t > 0 every bond e^((t/2) l_i l_j) − 1 is positive, so κ_t(n) > 0 at every integer n ≥ 2; for t < 0, κ_t(pq) < 0. At first order in t the weight at pq is (t/2) ln p ln q, which is the convolution Λ∗Λ in ζ″/ζ: the classical case.

**3.3 The lines (heuristic; measured).** Landau's formula says Σ_{γ ≤ T} x^(iγ) = −(T/2π) Λ(x)/√x + O(ln T): the zeros have lines at prime powers only. Applying the same residue argument to the effective formula for H_t gives, for the zeros γ_k(t) of H_t in a window:

  Σ_k e^(iγ_k(t) ln n) = −(1/2π) · ln(n) κ_t(n) · n^(−1/2) · ∫ (γ/2π)^(−(t/4) ln n) dγ · e^(−i s ln n),  s = t(−π/8 + 7/(8γ)).   (4)

So for t > 0 **the zeros have a line at every integer**, with the cluster weights. The phase is a rigid shift of all zeros by −πt/8 (the imaginary part of the exponent shift), with a correction 7/(8γ) (reported by the compute agent, from the smooth density near the origin).

Measured (loop 4 and the compute agent; the first 12,000 zeros moved by (1), all zeros and their mirror images included, the zeros above as a smooth density):

- The moved zeros are zeros of H_t: checked without the flow equation, from H_t(2γ) = (1/8)·E_G ξ(½ − √(t/2) G + iγ) (G standard normal; Gauss–Hermite, mpmath), six zeros at t = 0.3 and 0.6: within 7.4·10⁻⁶ (the zeros had moved by up to 0.36).
- Sharp window (zeros 2,000–10,000; noise taken at half-integers): at t = 0 nothing at n = 6, 15, 21 (−1.3, 1.5, −2.5); at t = 0.3: −48.1, −82.5, −91.7 against −47.1, −82.9, −90.4 predicted; pooled slope over the 139 composites n ≤ 200 that are not prime powers: 0.984 ± 0.016 (t = 0.1), 0.985 ± 0.009 (t = 0.3), 0.980 ± 0.009 (t = 0.6); primes 0.996–0.999, prime powers 0.992–1.001.
- Smooth window (agent's code [agent_zeta/](agent_zeta/), Kaiser taper over the middle 60 %, 7,200 zeros): pooled slope for the composites 0.999992 ± 0.0000074 at t = 0.2 and 0.999982 ± 0.0000048 at t = 0.4; primes 0.999998 ± 0.0000018 and 0.999993 ± 0.0000020. The agent reports the same with 30,000 zeros moved (deviation at most 3.3·10⁻⁵, halving with the height, so of order 1/T; not rerun by the validation script). So the 1.5 % deficit of the sharp window was the window.

## 4. Matrices: an exact formula and its limit (proved)

Let U be Haar unitary N × N with eigenvalues z_i, E(x) = ∏(1 + x z_i) = Σ e_j x^j, and E_t(x) = Σ_j e^(−t j(N−j)) e_j x^j (the same flow in the angle, up to a constant). Let p_k(t) be the k-th power sum of the negated reciprocal roots of E_t. Proofs in [theory/proof.md](theory/proof.md).

**Theorem 1 (cluster expansion; any N, k, complex t).** p_k(t) = e^(−tNk) Σ_{λ ⊢ k} c_λ(t) p_λ with c_λ = (−1)^(r−1) (k/z_λ) ∏_s e^(t i_s²) Φ_λ(t), Φ_λ = Σ over connected graphs on the r parts of ∏_edges (e^(2t i_s i_s′) − 1). Partitions play the role of factorisations: a part i is a "prime" of size i.

**Circle.** For t ≥ 0 all roots of E_t lie on a circle (Lee–Yang; also from the flow of the angles, which repel), so p_k(t) are power sums of points on the unit circle after the radius is divided out.

**Theorem 2 (exact second moment; k ≤ N).** E|p_k(t)|² = e^(−2tNk) Σ_{λ ⊢ k} c_λ(t)² z_λ. (From E p_λ p̄_μ = δ_λμ z_λ for |λ|, |μ| ≤ N, Diaconis–Evans 2001, Thm 2.1(a).)

**Theorem 3 (limit).** With t = τ/N², k = uN, N → ∞ and 0 < u ≤ 1:

  E|p_k(t)|²/N → K_τ(u) = u² e^(−2τu) Σ_{r≥1} (1/r!) ∫_{a_i > 0, Σa_i = u} ∏_i (e^(2τa_i²)/a_i) · Φ_r(a)² da₁…da_{r−1},
  Φ_r(a) = Σ over connected graphs on r points of ∏_edges (e^(2τ a_i a_j) − 1).

The first term is u e^(−2τu(1−u)); to second order K = u − 2τu²(1−u) + τ²[2u³(1−u)² + u⁵/3] (from the proof; checked against the series at τ = 0.1: 0.47571 against 0.47573 at u = 0.5). The same integral comes out of the integers: Σ_p f(ln p/L)/p → ∫ f(a) da/a, and 1/z_λ → ∏ da_i/a_i / r!. Beyond u = 1 no formula is known to us.

Computed (loops 1–3):

| check | result |
|---|---|
| Theorem 2 against simulation, N = 30, 4,000 matrices, k = 1…12, τ = 1, 3, 6 | 36 values, largest deviation 2.7 standard errors |
| the limit against finite N = 100 (exact moments, 20,000 Gaussian samples), τ = 1, 4, u = 0.1…0.9 | 18 values, largest 2.5 standard errors; e.g. τ = 4, u = 0.9: 1.1480 ± 0.0084 against 1.1541 (clusters carry 62 %) |
| the peak at u = 1 (N = 100, 1,500 matrices) | τ = 1: 1.246 ± 0.033 (series 1.240); τ = 4: 3.310 ± 0.081 (series 3.259) |
| spacing variance of the flowed eigenvalues (N = 100) | 0.1795, 0.1140, 0.0395, 0.0025 at τ = 0, 1, 4, 16; every matrix keeps N zeros on the circle |

## 5. One family for ζ and matrices (conjecture; measured)

**Conjecture.** In the scaling t = 4τ/ln²(T/2π) the zeros of H_t near height T have the local statistics of Haar unitary matrices flowed to t = τ/N²: one family of point processes from the sine process (τ = 0) to the lattice (τ = ∞), with form factor K_τ(u) for 0 < u ≤ 1. At τ = 0 this is Montgomery's conjecture.

Test (loops 6 and 8): 9 blocks of 3,333 zeros from Odlyzko's tables near zero 10¹², 10²¹, 10²² (L = 24.5, 44.6, 46.8), moved by (3) with the zeros outside a block as a uniform continuum; statistics on the middle 2,000 zeros of each block; matrices N = 200.

| τ | spacing variance, ζ | matrices | line of p = 2, 3, 5 relative to τ = 0 (prediction (2)) |
|---|---|---|---|
| 0 | 0.1779 ± 0.0019 | 0.1772 ± 0.0016 | |
| 0.5 | 0.1397 ± 0.0014 | 0.1393 ± 0.0012 | 0.9918 (0.9905), 0.9854 (0.9852), 0.9786 (0.9787) |
| 1 | 0.1130 ± 0.0011 | 0.1129 ± 0.0010 | 0.9836 (0.9811), 0.9711 (0.9705), 0.9577 (0.9579) |
| 2 | 0.0770 ± 0.0007 | 0.0772 ± 0.0007 | 0.9678 (0.9625), 0.9432 (0.9420), 0.9174 (0.9176) |
| 4 | 0.0388 ± 0.0004 | 0.0393 ± 0.0004 | 0.9373 (0.9266), 0.8901 (0.8877), 0.8420 (0.8426) |

Form factor of the same blocks in bins of u (loop 8), against the series: the seven bins from 0.15 to 0.85 agree within 1.2 errors at τ = 0, 1 and 4 (e.g. bin 0.45–0.55: 0.3082 ± 0.0080 against 0.3097 at τ = 1, 0.0840 ± 0.0024 against 0.0848 at τ = 4). The bin 0.85–0.95 is low at every τ, also at τ = 0 against Montgomery's u (0.8766 ± 0.0106 against 0.9; at τ = 4: 1.1419 ± 0.0176 against 1.2056, the registered reading failed there); the lowest bin likewise (0.0936 ± 0.0025 against 0.1 at τ = 0). The peak at u = 1: 1.2237 ± 0.0192 (τ = 1), 2.9828 ± 0.0334 (τ = 4) in the bin 0.95–1.05.

At low height the comparison fails (loop 5): for zeros 6,000–10,000 (L = 7) the spacing variance is 0.1565 at t = 0 against 0.1796, and the ratio to matrices falls from 0.871 to 0.720 at τ = 7.7. There the primes are few and discrete; the conjecture is a statement about the limit.

## 6. The other side, t < 0 (matrices measured; ζ a prediction)

Backwards two neighbours attract; alone they obey g²(τ) = g₀² + 2τ/π² (g in mean spacings; from (3)). A pair with g₀ < √(2|τ|)/π has collided and left the line. With the small-gap law of the sine process, density (π²/3)s²:

  share of zeros off the line at τ < 0 ≈ (2/(9π)) (2|τ|)^(3/2).                     (5)

Measured on Haar matrices (N = 100, 4,000 matrices, loop 7): ratio to (5): 0.92 ± 0.07, 0.96 ± 0.04, 0.935 ± 0.024 at τ = −0.03, −0.06, −0.12; then 0.906, 0.837, 0.713, 0.523 at τ = −0.25, −0.5, −1, −2 (the neighbours slow the collisions; at τ = −2 the share is 0.296). For ζ, if the Riemann hypothesis and the small-gap law hold, (5) with τ = t ln²(T/2π)/4 would count the zeros of H_t off the line for t < 0: a quantitative form of Λ ≥ 0. Not computed for ζ here.

## 7. What this says about the Riemann hypothesis

Stated plainly: it gives no proof and no new route to one.

- A reading: the zeros are real for t ≥ Λ, and the cluster weights κ_t(n) are all ≥ 0 exactly for t ≥ 0. With Λ ≥ 0 known, the Riemann hypothesis says the two thresholds coincide: the zeros stay real as long as the interaction between the primes is not repulsive, and the free gas (the Euler product) sits on the edge.
- The tempting principle "for deformations of ζ by a universal factor, real zeros ⟺ non-negative log-coefficients" is **refuted** (LOG.md, between loops 2 and 3): the factor 1 + cu² preserves real zeros but gives the term n the weight 1 + c(L/2 − ln n)²/4, whose log-coefficient at pq is negative for small p, q and large c; and Dirichlet characters have signed weights χ(n)κ_t(n). The coincidence at t = 0 is special to the Gaussian semigroup.
- What the flow does give: an exactly solvable family around ζ in which the arithmetic (integers, factorisations) and the matrix side (partitions, cycle types) produce the same series term by term.

## 8. What failed (all recorded in LOG.md)

- Loop 1: "clusters carry more than 10 % at τ = 3, k = 12": 6 %.
- Loop 2, first run: roots of degree-100 polynomials by numpy lost all digits; replaced by exact Gaussian moments (loop 2b) and by FFT with sign changes (loop 3).
- Loop 5: ζ at L = 7 against matrices: off by 13–28 %.
- Loop 7: the direction of the correction to (5) (guessed higher, is lower).
- Loop 8: one cell of 14 beyond 3 errors (bin 0.85–0.95 at τ = 4), a bin that is already 2.2 errors low at τ = 0.
- The compute agent's registered decision rule fails formally at t = 0.6 (errors of 3·10⁻⁶ against a 1/T term of 10⁻⁵).

## 9. Open

- A proof of (4) with an error term, for fixed t > 0 (the effective approximation of Polymath 15 should give it; unconditional where all zeros are real, t ≥ 0.2).
- The form factor beyond u = 1 and the full pair correlation of the family.
- The correction to (5), and a count on H_t itself for t < 0.
- Correct the side remark on k > N in [theory/proof.md](theory/proof.md) (the referee's counterexamples for t < 0) and cite Diaconis–Evans eq. (2.1) for the truncated formula.
- A validation script that reruns every loop and checks each quoted number (not written yet).

## Validation

No single validation script yet. Each number above was printed by the script named with it, with fixed seeds: [loop1.py](loop1.py) (`30 4000 12 1 3 6`), [loop2b.py](loop2b.py) (`100 20000 0 1 4`), [loop3.py](loop3.py) (`100 1500 0 1 4 16`), [loop4_lines.py](loop4_lines.py) and [loop4_check.py](loop4_check.py) on the run of [agent_zeta/run.py](agent_zeta/run.py) (`12000 0.3 0.005 runA.npz 0.6`, about 10 s), [agent_zeta/main.py](agent_zeta/main.py) (`runA.npz 0.2 0.8 0.2,0.4`, rerun on 27 November with the numbers of section 3.3), [loop5.py](loop5.py), [loop6.py](loop6.py), [loop7.py](loop7.py) (`4000 0 -0.03 -0.06 -0.12 -0.25 -0.5 -1 -2`), [loop8.py](loop8.py). Proof checks: [theory/check.py](theory/check.py) (41 checks, 0 failed) and the referee's own code [theory/referee/ref.py](theory/referee/ref.py) (85 checks, 0 failed). Every run under `ulimit -v 1500000` with `timeout`. Zeros: [odlyzko_zeros1.txt](../../../data/riemann/odlyzko_zeros1.txt), [odlyzko_zeros3.txt](../../../data/riemann/odlyzko_zeros3.txt), [odlyzko_zeros4.txt](../../../data/riemann/odlyzko_zeros4.txt), [odlyzko_zeros5.txt](../../../data/riemann/odlyzko_zeros5.txt).
