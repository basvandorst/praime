# Do consecutive primes remember their non-abelian labels?

**In short.** The label of a prime is its Frobenius class in a non-abelian extension (number of roots of x³ − x − 1 or x³ + x − 1 mod p; the ideal class of p in Q(√−47) and Q(√−39); the angle sector of p = a² + b²; the sign, the parity and a_p mod 3 of the curve 37a). Each label has an abelian part (a congruence class: (p/23), p mod 4, p mod 3, ...) and a non-abelian part.
1. The Hardy–Littlewood expansion that gives Lemke Oliver and Soundararajan's (LOS) biases cannot couple the non-abelian part. Its main terms are products of major-arc pieces. Λ(n)a_ψ(n) has none when the theta series (or Hecke or Sym^k form) behind a_ψ is cuspidal, which is the case exactly when the label is not a congruence condition. So the null is **H0 = N2**: LOS on the abelian part, and independent Chebotarev draws inside each fibre. It holds at every lag, every gap and every order in 1/log x. Heuristically, the corrections are of relative size x^(−1/2).
2. On the discovery range (all 5,761,454 odd primes below 10⁸, eight label systems) nothing departs from H0:
   - 75 fibre-interaction cells: max |z| = 2.37, Σz² = 61.9;
   - the full tables against N2 with their own abelian part: all p ≥ 0.39;
   - the gap-conditioned tables, the run lengths and the triples agree;
   - the within part of the block transition tables is 1.0002 ± 0.0030 (non-abelian), against 0.91–0.95 for the abelian sign (the LOS repulsion).
   Any non-abelian memory is below about 4 % of the abelian one.
3. The only large departures from "Poisson" are in the column part, that is the class counts per block. There the non-abelian contrasts are Poissonian up to a crossover H* = e^B (q x^(d−1))^(1/d) (B = 1 − γ − ln 2π, d the degree, q the conductor) and rigid beyond it. This is the degree-d variance law of Bui–Keating–Smith (2016), made local. It fits every scheme, block length and height bin with no fitted parameter. It is known in substance and is not a memory of consecutive primes.
4. The holdout is frozen: predictions_consec.json, sha256 f3e3599d43a6d5d9e73eefa1a46a1a96686ca88c58b3ce07623188096b4f3cb1. The ranges are (10⁸, 10⁹] for the five lattice and Hecke schemes and (10⁸, 4·10⁸] for 37a. A non-abelian correlation of |r| ≈ 1.5·10⁻³ would show at 5σ.

Status tags: **proved** (algebra), **conjecture (HL)** (from the Hardy–Littlewood k-tuple conjecture or LOS), **heuristic**, **computed** (a number from a script here), **known** (in print).

## 1. The labels (computed; checked)

| scheme | field / object | label codes | abelian part (fibre) | non-abelian functions (mean 0, mean square 1 in the fibre) |
|---|---|---|---|---|
| C23 | splitting field of x³ − x − 1 = Hilbert class field of Q(√−23), S₃ | 0 inert (1 root), 1 = x²+xy+6y² (3 roots), 2 = 2x²±xy+3y² (0 roots) | (p/23) | f = a_ψ(p)/√2 ∈ {√2, −1/√2}, ψ of order 3 |
| C31 | x³ + x − 1, Q(√−31), S₃ | 0 inert, 1 = x²+xy+8y², 2 = 2x²±xy+4y² | (p/31) | same |
| C47 | Hilbert class field of Q(√−47), D₅ | 0 inert, 1 principal (1,1,12), 2 = (2,±1,6), 3 = (3,±1,4) | (p/47) | a_ψ/√2, a_ψ²/√2 (ψ of order 5) |
| C39 | Hilbert class field of Q(√−39), D₄ | 0 inert, 1 principal (1,1,10), 2 = (3,3,4), 3 = (2,±1,5) | genus: (p/3), (p/13) | ±1 on the principal genus (ψ of order 4) |
| G | Q(i): p = a² + b², a > b > 0, θ = atan(b/a) ∈ (0, π/4) | 0 for p ≡ 3 (4), 1..4 = sector of width π/16 | p mod 4 | √2 cos 4θ, √2 cos 8θ (Hecke characters ξ⁴, ξ⁸) |
| E37s | 37a1: y² + y = x³ − x, non-CM | sign of a_p (a_p = 0 excluded from the fibre) | none | sign, U₁ = a_p/√p, U₂ = U₁² − 1 |
| E37m2 | a_p mod 2 (Gal Q(E[2]) = S₃, quadratic subfield Q(√37)) | 0 if (p/37) = −1 (a_p even), 1 even, 2 odd | (p/37) | (1[even] − 1/3)/√(2/9) |
| E37m3 | a_p mod 3 (GL₂(F₃)) | 3·[p ≡ 2 (3)] + (a_p mod 3) | p mod 3 | centred 1[a_p ≡ 0], 1[a_p ≡ 1] − 1[a_p ≡ 2] |

Chebotarev inside the fibres: P(3 roots | split) = 1/3; D₅: 1/5, 2/5, 2/5; D₄: 1/2, 1/2; GL₂(F₃): P(a_p ≡ 0, 1, 2) = (1/4, 3/8, 3/8) for p ≡ 1 (3) and (1/2, 1/4, 1/4) for p ≡ 2 (3). These are counted from the 48 matrices, by determinant and trace.

**Checks.** The script is [check_lib.py](check_lib.py) → [check_lib_out.txt](check_lib_out.txt): 0 errors. It covers all odd p < 3·10⁴ and every third prime in (10⁸ − 10⁵, 10⁸):
- the primes, against sympy;
- the cubic root counts, by PARI polrootsmod;
- C47 and C39, by the roots of the Hilbert class polynomial polclass(D) and by brute-force representation with plain integers;
- the Gaussian angles, by a² + b² = p;
- a_p of 37a, by a point count for p < 6000.

## 2. The null, derived

### 2.1 Decomposition (proved)
Every label indicator is a class function on G = Gal(L/Q). So it is a combination of irreducible characters: the trivial one, the abelian ones (Dirichlet characters, which give congruence conditions) and the non-abelian ones of dimension ≥ 2. For the dihedral cases L is the Hilbert class field of K = Q(√D). The 2-dimensional characters are Ind ψ with ψ ∈ Cl(K)^ and ψ² ≠ 1, and tr Ind ψ(Frob_p) = a_ψ(p) = ψ(𝔭) + ψ(𝔭̄) for split p, 0 for inert p. For ψ² = 1, ψ ≠ 1 (a genus character), Ind ψ splits into two Dirichlet characters: the abelian part. The joint law of (label(p_n), label(p_(n+k))) is therefore a sum of correlations E[χ_π(Frob p_n) χ_π′(Frob p_(n+k))] over pairs of irreducibles. The pairs in which both are abelian are LOS's correlations. The others are the non-abelian memory.

### 2.2 N1 and N2 (conjecture (LOS) for the abelian part)
- **N1:** P(a, b) = P(a)P(b) with the Chebotarev marginals.
- **N2:** P(a, b) = P_LOS(χ(a), χ(b)) · P(a | χ(a)) · P(b | χ(b)).

For a real character χ of prime conductor q (and q = 4), LOS (arXiv 1603.03720v4, Main Conjecture p. 2–3) give π(x; q, a) = li(x)/φ(q)^r · (1 + c₁ lnln x/ln x + c₂/ln x + O((ln x)^(−7/4))), with c₁(q; a) = (φ(q)/2)((r − 1)/φ(q) − #{i : a_i ≡ a_(i+1)}). Projected on χ ⊗ χ, this is their Conjecture 1.3:

  Σ_(p_n ≤ x) χ(p_n)χ(p_(n+1)) = −x/(2 (log x)²) · log(2π log x/q) + O(x/(log x)^(11/4)),

and Conjecture 1.4 gives −li(x)/(2(k − 1) log x) for k ≥ 2. (Check, proved: the c₁ part projects to −½·lnln x/ln x. The c₂ part, from (1.1)–(1.2), projects to −½·log(2π/q)/ln x. Both were redone for q = 4 by hand.) Over a range the cumulative forms are differenced. With P(χ_n = s, χ_(n+k) = t) = (1 + st·b_k)/4 this gives ([consec_null.py](consec_null.py), computed):

| q | b₁ on [3, 10⁸] | b₁ on (10⁸, 10⁹] | b₂, b₃ on (10⁸, 10⁹] |
|---|---|---|---|
| 23 | −0.04386 | −0.04157 | −0.02374, −0.01187 |
| 31 | −0.03576 | −0.03448 | same |
| 47 | −0.02446 | −0.02460 | same |
| 4 (Gaussian) | −0.09134 | −0.08310 | same |
| 37 | −0.03096 | −0.03050 on (10⁸, 4·10⁸] | −0.02454, −0.01227 |
| 3 (E37m3) | −0.09915 | −0.09215 on (10⁸, 4·10⁸] | same |

For C39 (genus characters mod 3 and 13; composite modulus 39) the cross terms need LOS's (2.23). They were not evaluated, so C39 uses only N2-emp: the measured abelian table times Chebotarev inside the fibres.

### 2.3 Why there is no non-abelian term (heuristic, from conjecture (HL))
- **The expansion.** Write the "consecutive" condition as LOS do: by inclusion–exclusion over the sites between p_n and p_(n+k). Every term is then a correlation Σ_n F₀(n) F_h(n + h) Π_(t∈T) Λ(n + t), with F₀, F_h ∈ {Λ, Λχ, Λa_ψ}. The Hardy–Littlewood main term of such a correlation is the average of the product of the major-arc approximants of its factors. For Λ that approximant is Σ_q μ(q)c_q(n)/φ(q), the expansion behind the project's resonance constants.
- **Λa_ψ has none.** Σ_(n ≤ x) Λ(n)a_ψ(n)e(an/q) is a combination of Σ Λ(n)a_ψ(n)χ(n) over Dirichlet characters χ mod q. Now θ_ψ ⊗ χ = θ_(ψ·(χ∘N)). This is a cusp form unless ψ·(χ∘N) factors through the norm, which needs ψ² = 1. So for ψ of order ≥ 3, L(s, θ_ψ ⊗ χ) has no pole and does not vanish on Re s = 1 (proved), and every such sum is o(x).
- **Consequence.** Every term carrying a factor Λa_ψ has main term 0: in c₁, in c₂ and at every further order of the LOS expansion. The same holds for the Gaussian Hecke characters ξ^(4m), m ≥ 1, which are cusp forms of weight 4m + 1, and for U_k(cos θ_p) of 37a: Sym^k is cuspidal automorphic (Newton–Thorne), and the sign of a_p is an L² series in the odd U_k. For a_p mod ℓ it holds for the non-abelian irreducibles of GL₂(F_ℓ): Artin L-functions have neither pole nor zero on Re s = 1.
- **The genus characters.** For ψ² = 1, θ_ψ is an Eisenstein series, a_ψ(p) = χ_(D₁)(p) + χ_(D₂)(p), and the term is LOS's abelian term.
- **Lattice cross-check.** #{(u, v): Q_C(u) = p, Q_C′(v) = p + h, both prime} is a count on the quaternary quadric Q_C′(v) − Q_C(u) = h. Its singular series with the primality conditions, and its singular integral, depend only on the ℓ-adic and real classes of Q_C, Q_C′. All forms in one genus are locally equivalent, so the main term is the same for every pair of classes, and it cancels in Σ_C ψ(C). (The total over classes is what LOS use.) For several genera (−39) the main term depends on the genus: that is the abelian part.
- **What is left.** Only minor-arc terms remain, of heuristic relative size x^(−1/2+ε). One explicit piece is known in the project: the radial resonance of fact 2055, Λa_ψ(n) ∋ (3/4)K n^(−1/4) e(2√(n/q)) + c.c., |K| ≈ 0.4–0.6. Its autocorrelation at a shift h ≪ √x is about |K|² x^(−1/2)/log x ≲ 10⁻⁵ at 10⁸. This is far below the errors here (8·10⁻⁴).

**H0** (conjecture): N2 holds exactly at every lag k, for every gap and every pattern of the intermediate primes, to all orders in 1/log x.

**H1**: none. No non-factorising singular series exists at fixed shifts.

### 2.4 What H0 does not cover: long range, the column part (known in substance)
H0 is a statement at fixed shifts. Summed over shifts up to H, the zeros of L(s, θ_ψ) make the block counts of a non-abelian contrast rigid once H exceeds a crossover. Take a primitive L-function of degree d and conductor q. The variance of Σ_(x<n≤x+H) Λ_F(n), divided by its diagonal H ln x, is
- for d = 1 (Montgomery–Soundararajan, as in round 135): R = [ln(x/H) + ln q + B]/ln x at every H;
- for d ≥ 2 (Bui–Keating–Smith 2016, Theorem C1, converted to a local form by d/dX): R = [d ln(x/H) + ln q + dB]/ln x above the crossover, and R = 1 below it (H0: only the diagonal). In the integrated form, their constant −(γ₀ + log 2π)d becomes dB locally.

So R = min(1, ·), with crossover H* = e^B (q x^(d−1))^(1/d). At x = 5·10⁷: H* ≈ 3,400 for G (q = 4), 8,200 for C23, 10,500 for the U₁ of 37a; and ≈ 2·10⁶ for Sym² (d = 3, q = 37²). The sign of a_p mixes U₁ (weight c₁² = (8/3π)² = 0.72) with higher odd U_k. Script: [consec_bks.py](consec_bks.py).

### 2.5 The abelian within part (conjecture (LOS); round 137)
For the 2 × 2 table of abelian signs, the within part of the per-block χ² is 1 + 2ρ₁, with ρ₁ = (b₂ − b₁²)/(1 − b₁²) (project derivation within/derivation.md §2.5). With the LOS values this is 0.942 (q = 23) on the discovery range and 0.949 (q = 23), 0.950 (31), 0.951 (47), 0.938 (4) on the holdout. For the non-abelian fibre tables H0 gives exactly 1, since every τ_j vanishes when the labels are independent given the abelian data.

## 3. Measurements on the discovery range (computed)
Scripts: [consec_tables.py](consec_tables.py) (segments of 5·10⁶ → disc/seg_*.npz), [consec_analyze.py](consec_analyze.py) → [readings_disc.json](readings_disc.json), [consec_report.py](consec_report.py) → [report_disc.txt](report_disc.txt), [consec_fulltable.py](consec_fulltable.py) → [fulltable_disc.txt](fulltable_disc.txt), [consec_bks.py](consec_bks.py) → [bks_disc.txt](bks_disc.txt).

**(a) Full tables (lag 1–3): N1 against N2.** χ² of the label tables:

| scheme | N1 (k = 1) | N2-LOS (k = 1) | N2-emp k = 1, 2, 3 (p) |
|---|---|---|---|
| C23 | 11,833 / 8 | 13.4 / 8 | 1.0/5 (0.96), 3.1/5 (0.69), 5.2/5 (0.39) |
| C31 | 5,461 / 8 | 147 / 8 | 0.55, 0.88, 0.86 |
| C47 | 4,143 / 15 | 38.7 / 15 | 0.84, 0.67, 0.79 |
| G | 73,628 / 24 | 2,746 / 24 | 0.98, 0.99, 0.89 |
| E37m2 | 5,647 / 8 | 4.4 / 8 | 0.60, 0.56, 0.59 |
| E37m3 | 86,152 / 35 | 3,134 / 35 | 0.83, 0.61, 0.69 |

N1 is rejected by the abelian memory. N2 with LOS's asymptotic numbers removes 90–99.9 % of the χ², but not all of it: measured b₁/LOS is 1.03 (23), 0.86 (31), 1.10 (47), 1.24 (q = 4), 1.01 (37), 1.23 (q = 3), within LOS's O((log x)^(−7/4)) remainder. With the data's own abelian table and Chebotarev inside the fibres (N2-emp) every table fits. The summed N2-emp χ² is 176 on 240 dof, which is low. Its source is the marginals: all 28 marginal z-values have |z| < 0.7. This is the column rigidity at H ≈ x (§2.4), the same "too good" fit of single Frobenius frequencies that Ash–Bate–Gross observed (heuristic explanation).

**(b) Fibre interactions (H0: 0).** r is the correlation of f(p_n) and g(p_(n+k)) over the pairs in which both lie in the fibre, z = r√n:

| scheme | lag-1 pairs | r (lag 1, f₀f₀) | all cells k = 1..3 (|z| max) |
|---|---|---|---|
| C23 | 1,374,755 | −0.00053 ± 0.00085 | 1.66 |
| C31 | 1,395,846 | +0.00046 ± 0.00085 | 0.82 |
| C47 | 1,401,647 | −0.00001 ± 0.00084 | 2.02 |
| C39 | 294,771 | +0.00178 ± 0.00184 | 0.97 |
| G | 1,277,323 | −0.00081 ± 0.00088 | 1.34 |
| E37s | 5,760,276 | −0.00022 ± 0.00042 | 1.85 |
| E37m2 | 1,395,121 | +0.00052 ± 0.00085 | 0.62 |
| E37m3 | 5,761,450 | −0.00036 ± 0.00042 | 2.37 |

Over all 75 cells Σz² = 61.9 and max |z| = 2.37. For comparison, the abelian b₁ is −0.045 (q = 23) and −0.11 (q = 4). In plain words: given that p_n and p_(n+1) both split in Q(√−23), whether each has 3 roots or none is decided afresh. Any memory is below 0.0017 in correlation (2σ), against 0.045 for the abelian part.

The same holds in four further forms:
- for the k-th next fibre prime (corrM, all |z| < 1.9);
- as contingency χ² inside every abelian fibre pair (102 dof pooled, χ² 99.7, p = 0.55);
- against the neighbour's abelian datum (66 cells, max |z| = 1.97);
- conditioned on the gap: Σ over the gaps 2–198 with ≥ 400 pairs, p = 0.08–0.93; the gap ≡ 0 mod q contrasts |z| ≤ 2.05. For G every fibre pair has gap ≡ 0 mod 4, so that contrast is empty.

**(c) Second order.**
- Triples (first–third interaction given the middle label): p = 0.031 (C23), 0.23–0.93 for the others.
- Run lengths in the fibre subsequence against the geometric law: p = 0.011 (C31, class "0 roots": χ² 40.0/23), 0.10–0.99 for the others. The G sector-2 runs give 20.3/8 on their own. Mean run lengths agree with the geometric ones to 10⁻³.

**(d) Block transition tables, H = 10⁴, 10,000 blocks (r135b method).** The within part of the non-abelian fibre tables is 1 under H0:

| scheme | within (fibre) | within (abelian sign) | LOS value 1 + 2ρ₁ |
|---|---|---|---|
| C23 | 0.991 ± 0.014 | 0.932 ± 0.013 | 0.942 |
| C31 | 0.995 ± 0.014 | 0.932 ± 0.013 | 0.943 |
| C47 | 1.010 ± 0.007 | 0.913 ± 0.013 | 0.945 |
| C39 | 0.970 ± 0.013 | 0.889 ± 0.013 | – |
| G | 1.004 ± 0.005 | 0.951 ± 0.013 | 0.928 |
| E37m2 | 1.003 ± 0.014 | 0.938 ± 0.013 | 0.944 |
| E37s / E37m3 | 0.976 ± 0.014 / 1.000 ± 0.007 | – | – |

Pooled over the non-abelian fibre tables: **1.0002 ± 0.0030**. Pooled over the abelian sign tables (without C39): about 0.933 ± 0.006. That is the LOS repulsion between p_n and p_(n+2) acting on the abelian half only.

**(e) Column part (block dispersion R = ΣT²/ΣV).**
- H = 100: pooled R = 0.99980 ± 0.00038 over all non-abelian functions.
- H = 1000: 0.9977 ± 0.0012. The BKS law gives 0.996–0.999 here, because the low heights are already past the crossover.
- H = 10⁴, 10⁵, 10⁶: R falls to 0.6–0.8 at 10⁵ and to 0.3–0.6 at 10⁶. This follows the law of §2.4 in every scheme and height bin (bks_disc.txt): for bins with ≥ 50 blocks, |z| ≤ 3 except at the crossover itself (U₁ of 37a at H = 10⁴ on [3·10⁷, 10⁸): −3.6).
- The degree shows. Sym² of 37a (U₂, d = 3) stays at 0.93 ± 0.04 at H = 10⁵, where U₁ (d = 2) is at 0.68 ± 0.03. The sign sits in between (0.78 ± 0.03; law 0.79).
- The abelian contrasts follow the d = 1 law from H = 1000 on. At H = 100 they lie 0.01–0.02 above it, a known lower-order effect of very short blocks, not pursued. One abelian cell misses by −4.7 errors (E37m2 χ, H = 10⁵, [10⁷, 3·10⁷), 400 blocks): not pursued. The delta-method errors may be too small for heavy-tailed blocks.

## 4. Candidates
- No reading of the consecutive-prime memory departs from H0 by 4 errors or more; the largest is 2.37.
- The only readings beyond 4 errors are column-part readings at large H. They are the BKS rigidity (§2.4), predicted without a parameter, and they are not memory.
- Four marginal items are registered as watch items, each with a real-effect and a fluctuation prediction:
  - W1: C31 runs, χ² 40.0/23;
  - W2: C23 triple, 8.9/3;
  - W3: E37m3 lag-2 f₀f₀, z = +2.37;
  - W4: G sector-2 runs, 20.3/8.
  A real effect scales its χ² excess (or its z) with the sample: 7.8 times on the main holdout, √2.7 for 37a. A fluctuation returns to χ²(dof) or N(0, 1).

## 5. Frozen predictions for the holdout
Holdout ranges: (10⁸, 10⁹] for C23, C31, C47, C39, G (45,086,079 primes, 7.83 times the discovery sample) and (10⁸, 4·10⁸] for E37s, E37m2, E37m3 (15,574,871 primes, 2.70 times). [predictions_consec.json](predictions_consec.json), sha256 `f3e3599d43a6d5d9e73eefa1a46a1a96686ca88c58b3ce07623188096b4f3cb1`, written by [consec_freeze.py](consec_freeze.py) from readings_disc.json and the two law scripts only. [consec_test.py](consec_test.py) carries the hash and refuses any other file (checked: a one-character change is refused).

Decisions (fixed in predictions_consec.json):
- **D1a.** Pooled contingency χ² inside the fibres: pass if p ≥ 10⁻³. *Falsifier:* p < 10⁻⁵.
- **D1b.** Every interaction cell has |z| < 4.5. *Falsifier:* any |z| ≥ 5.
- **D1c.** The pooled within part is within 4 errors of 1, each scheme within 4.5. *Falsifier:* the pooled value is more than 5 errors from 1.
- **D2.** Gap-conditioned χ²: p ≥ 10⁻³; gap ≡ 0 mod q contrasts |z| < 4.
- **D3.** Neighbour-datum cells: |z| < 4.5.
- **D4.** Runs: p ≥ 10⁻³.
- **D5.** Triples: p ≥ 10⁻³.
- **D6.** Short-shift dispersion: H = 100 and 1000 give R within 4 errors of 1, per function and pooled.

H0 stands if all pass. Failures of D2–D6 without a falsifier are reported as open anomalies.

Secondary readings, which are known laws and not part of the verdict:
- b_k against LOS: the ratio must lie in [0.6, 1.5].
- The abelian within part against 1 + 2ρ₁: within 4 errors + 0.015.
- The column law at H = 10⁵ and 10⁶ per height bin: within 4√(err² + 0.03²). For example, C23 at H = 10⁵ on (10⁸, 3·10⁸] gives 0.808; G 0.716; 37a U₁ 0.833; Sym² 1.000. At H = 10⁶: C23 0.649 on (3·10⁸, 10⁹], Sym² 0.976, sign 0.705.

**Power.** The expected lag-1 fibre pairs are 10.8M (C23), 10.9M (C31, C47), 2.3M (C39), 15.6M (E37s, E37m3), 3.8M (E37m2), 10.0M (G). A non-abelian correlation of |r| = 1.5·10⁻³ (C23 type; 3.3·10⁻³ for C39, 1.3·10⁻³ for E37s and E37m3, 2.6·10⁻³ for E37m2) would give 5σ. That is 1/30 of the abelian b₁ for q = 23. The within part per scheme gets errors of 0.0025–0.0087. A non-abelian lag-2 repulsion of the LOS size (0.05) would be visible at > 10σ.

Dry run on the discovery range (not blind; [test_dry_out.txt](test_dry_out.txt)): D1a p = 0.55; D1b max |z| = 2.37; D1c 1.0002 ± 0.0030; D2–D5 pass; D6 H = 100: 0.99980 ± 0.00038, H = 1000: 0.9977 ± 0.0012 (z −1.9). H0 stands.

## 6. Prior art (mathematical queries only)
- **Ash, Bate, Gross, "Frequencies of successive tuples of Frobenius classes", Exp. Math. 18 (2009) 55–63 (read).** They tested k-serial frequencies of consecutive Frobenius classes against independent Chebotarev draws (that is, N1) for C₂, A₃, S₃ (among them x³ + x − 1, our C31), D₅ and A₅ fields, on the first 5·10⁷ primes and 10⁷ primes above 10¹⁰⁰.
  - Their finding is "nonrandom when K has an abelian subfield of small discriminant", A₅ fields random. They leave open Question 3: "Is the existence of an abelian subfield the only source of nonrandom behavior?". They also state that L-functions "do not appear to give information concerning frequencies of pairs of consecutive Frobenius classes".
  - They also report single Frobenius frequencies fitting Chebotarev "better than true randomness" (§4(a) here explains that heuristically). Their D₅ field 42 gives a 3-serial p > 0.999 on the 5-cycle subsequence: one test among many, no follow-up.
  - **Verdict:** the observation that the abelian part carries the bias is KNOWN IN SUBSTANCE (empirical, against N1). The conditional null N2 with LOS's terms, its derivation (H0), and a test of conditional independence given the abelian data were NOT FOUND.
- **Ash, Beltis, Gross, Sinnott, Exp. Math. 20 (2011)** (read): abelian only. **Lemke Oliver–Soundararajan, PNAS 2016 / arXiv 1603.03720v4** (read): abelian only. They mention sums of two squares and other sifted sets, nothing on Galois labels.
- **Variance of prime sums in short intervals for degree d: Bui–Keating–Smith, arXiv 1506.03741** (read in part; the C2 derivation was not shown in what was read). The two-regime law (transition at X^(1−1/d)), with numerics for τ and the conductor-37 curve at X = 10⁶. **KNOWN** for §2.4 and the column part. Their short-regime constant, read here as a local off-diagonal term, would give R ≈ 1 − 0.42/ln x = 0.979 at H = 100. The data give 0.9998 ± 0.0004, so either that reading of their normalisation is wrong or the constant does not apply to local unweighted counts (doubt 3).
- **Consecutive Gaussian primes in sectors:** NOT FOUND. Related one-point work: Devin, arXiv 2105.02492 (biases of Gaussian-prime angles; abstract only), Rudnick–Waxman (angles, variance; not opened).
- **Signs of a_p or λ_f(p) at consecutive primes:** NOT FOUND. Related work:
  - Matomäki–Radziwiłł (sign changes at consecutive integers);
  - Gillman–Kural–Pascadi–Peng–Sah, arXiv 1907.08285 (bounded gaps between Sato–Tate primes in an interval; abstract only; existence, no frequencies);
  - Kimura, arXiv 2101.05193 (spacing of Sato–Tate angles; not opened).
- **"LOS for Chebotarev / number fields":** NOT FOUND (several searches). Related work: Thorner, bounded gaps between primes in Chebotarev sets (existence; not opened); David–Devin–Nam–Schlitt, arXiv 2105.05048 (LOS bias for consecutive sums of two squares; abelian; not opened).
- **Function-field analogues:** NOT FOUND for consecutive irreducibles with Frobenius conditions. Related work: Bary-Soroker–Gorodetsky–Karidi–Sawin, arXiv 1810.06201 (Chebotarev in short intervals over F_q(T); one-point; abstract only); Bary-Soroker–Stix, arXiv 1711.05564 (cubic twin-prime polynomials counted with a modular-form term, anomalies for q ≡ 1 mod 3; abstract only, and not a statement about Frobenius labels).
- **Named sources not read in full** (they void a NOT FOUND only if they turn out relevant): Thorner (2014); David–Devin–Nam–Schlitt (2021); Devin (2021); Rudnick–Waxman (2019); Kimura (2021); Gillman et al. (2019, full text); the full BKS derivation of C2; Kaneko-type or Pólya (1959) heuristics cited by ABG. None of the abstracts read suggests a treatment of non-abelian labels of consecutive primes.

## 7. Doubts
1. H0 is a heuristic (HL plus the major-arc principle), not a theorem. The size x^(−1/2) of the remainder is a heuristic too. Only the resonance piece was estimated.
2. N2-LOS is not a sharp null at these heights: LOS's remainder is 15–25 % of b₁ (measured/LOS from 0.86 to 1.24). The sharp tests use the data's own abelian part (N2-emp, the within-fibre statistics). LOS's (2.23) was not evaluated for the composite modulus 39.
3. My local reading of BKS C2: the constant (3 + 8 log 2)/6 may belong to the integrated, Λ-weighted form with a convention I did not see. D6 tests "R = 1" (H0), not BKS.
4. The column law uses a kink min(1, ·); the true crossover is smooth, so H = 10⁴ is only described. Conductors used: 23, 31, 47, 39, 4 (Hecke, Q(i)), 37, 37² (Sym²), 37³ (lumped higher odd Sym^k for the sign), 148 (the Artin conductor of the S₃ field of 37a[2], = |disc| of z³ − 4z + 2). The archimedean shifts of the Hecke and Sym^k gamma factors are ignored. Errors at H = 10⁶ in bins with fewer than 50 blocks are not reliable (e.g. two discovery cells at −6σ and −21σ with 9 blocks); such bins are excluded from the registered band.
5. The pooled χ² and max-|z| rules treat cells as independent. They are not: E37s, E37m2 and E37m3 are functions of the same a_p, and sign and U₁ are correlated at 0.85. The thresholds (p ≥ 10⁻³, |z| < 4.5) leave room for that.
6. The labels at holdout heights rely on the same code as the checked ranges (int64 is safe to 5·10¹⁰, isqrt exact). There is no independent check above 10⁸ by design.
7. Pairs whose partner lies beyond X1 are dropped (nothing above X1 is computed), and so are runs that touch a segment edge: biases of order 10⁻⁵.
8. 37a runs only to 4·10⁸ in the holdout (ellap costs about 60 µs per prime; the other five schemes about 7 µs per prime).
9. The number π(4·10⁸) = 21,336,326 in predictions_consec.json is from memory, not computed. A mismatch only prints a warning.

## 8. Files and commands
All files are in this folder (the segment files of the tables are scratch files and are not stored). Run every job as `(ulimit -v 1200000; timeout 600 nice /home/bas/praime/.venv/bin/python -I -B <script> ...)`.
- consec_lib.py (labels), check_lib.py → check_lib_out.txt
- consec_tables.py OUTDIR X0 X1 [--cases ...] (segments of 5·10⁶, resumable, stops after 540 s): discovery `consec_tables.py disc 0 1e8` (20 segments, about 15–25 s each)
- consec_analyze.py DIR OUT.json → readings_disc.json; consec_report.py → report_disc.txt; consec_fulltable.py DIR X0 X1 → fulltable_disc.txt/json; consec_bks.py → bks_disc.txt
- consec_null.py X0 X1 (LOS numbers); consec_freeze.py → predictions_consec.json; consec_test.py (dry → test_dry_out.txt/json)
- holdout order: see registration.txt.
