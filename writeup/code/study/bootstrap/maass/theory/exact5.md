# The even real character problem: q*_∞ ≥ 5 − 2.8·10⁻¹⁴ from Ξ_{χ₅}²·B, unconditionally (27 November 2026)

**In short.** This is the problem of [../maass.md](../maass.md), section 3, for even characters: Gamma factor Γ_R(s), no pole, conductor q, one probability measure on [−1, 1] per prime. The construction of [exact3.md](exact3.md) is applied with χ = χ₅ = (·/5) and ĝ = Ξ_{χ₅}²·B.

1. **Proved.** (log 5) g(0) + A₀(g) = 2 Σ_p P_p(χ₅(p)) for every even B. This is the classical explicit formula applied to the kernel Ξ²B, with the zero side identically 0; it needs no hypothesis on the zeros. A family that is positive at conductor q for all f ∈ C_c^∞(ℝ) satisfies log q ≥ log 5 − 2 Defect(g)/g(0).
2. **Computer-assisted, one code.** Best B found (B = 1, then degree 16):
   - q*_∞ ≥ 4.9057 (B = 1);
   - q*_∞ ≥ 5 − 7.35·10⁻⁴ (degree 4);
   - q*_∞ ≥ 5 − 2.21·10⁻⁶ (degree 8);
   - q*_∞ ≥ 5 − 1.15·10⁻⁹ (degree 12);
   - q*_∞ ≥ 5 − 2.77·10⁻¹⁴ (degree 16).
   Only this ε is shown, not every ε > 0.
3. **Caveats.**
   - Without GRH for L(s, χ₅) no feasible family is known, so the statement may be vacuous.
   - Exactness (q*_∞ = 5) needs inf Defect = 0, which is unproved, and GRH for the upper bound.
   - The repo's LP value 4.9999999994 at support 24 is a finite-support number. It is not compared with, or bounded by, anything here.

## 1. What changes against exact3.md (each change written out; everything not listed is word for word the same)

| item | odd, χ₋₃ (exact3.md) | even, χ₅ (here) |
|---|---|---|
| χ | (n mod 3) ↦ 0, 1, −1 | (n mod 5) ↦ 0, 1, −1, −1, 1 (Legendre symbol mod 5) |
| Λ(s) | 3^{s/2} Γ_R(s+1) L(s, χ₋₃) | 5^{s/2} Γ_R(s) L(s, χ₅) |
| A_κ in the Weil form | A₁: kernel e^{−3u/2} | A₀: kernel e^{−u/2}, i.e. A₀(g) = −(γ + log π) g(0) + ∫₀^∞ [2e^{−2u} g(0) − 2e^{−u/2} g(u)]/(1 − e^{−2u}) du |
| theta function | θ(x) = Σ n χ(n) e^{−πn²x/3} | θ₀(x) = Σ_{n≥1} χ(n) e^{−πn²x/5} (no factor n) |
| (C1) transformation | θ(1/x) = x^{3/2} θ(x) | θ₀(1/x) = x^{1/2} θ₀(x): root number τ(χ₅)/√5 = 1 (classical, Davenport ch. 9, recalled; checked at x = 0.6, 1.7) |
| Mellin | Λ(s) = 3^{−1/2} ∫ θ(x) x^{(s+1)/2} dx/x | Λ(s) = ∫₀^∞ θ₀(x) x^{s/2} dx/x: termwise ∫ e^{−πn²x/5} x^{s/2−1} dx = Γ(s/2)(πn²/5)^{−s/2}; θ₀(x) = O(x^{−1/2} e^{−π/(5x)}) at 0 |
| φ, Ξ = φ̂ | (2/√3) e^{3u/2} θ(e^{2u}) | φ(u) = 2 e^{u/2} θ₀(e^{2u}), even by (C1) |
| φ*φ | (4/3) e^{3w/2} Σ χ(N) d(N) N K₀(2πNe^w/3) | **4 e^{w/2} Σ_N χ(N) d(N) K₀(2πN e^w/5)**: from φ(v+w)φ(v) = 4e^{w/2} Σ χ(nm) exp(−(π/5)(n²e^{2w}e^{2v} + m²e^{−2v})) and ∫ exp(−αe^{2v} − βe^{−2v}) dv = K₀(2√(αβ)); no factor nm, so the weight is d(N), not N d(N) |
| derivative recurrence (D = z d/dz) | D^j[z^{3/2}K₀]: p′ = (3/2)p + zp′ − zq, q′ = (1/2)q + zq′ − zp | D^j[z^{1/2}K₀] = z^{1/2}(p_j K₀ + q_j K₁): **p_{j+1} = (1/2)p_j + z p_j′ − z q_j, q_{j+1} = −(1/2)q_j + z q_j′ − z p_j** |
| h_j(n) = (d/dw)^j (φ*φ)(log n) | (4/3) n^{3/2} Σ χ d N [p_j K₀ + q_j K₁] | 4 n^{1/2} Σ_N χ(N) d(N) [p_j(z) K₀(z) + q_j(z) K₁(z)], z = 2πNn/5 |
| Lemma 1, prime side | log 3 | log 5; Λ′/Λ(s) = (1/2) log 5 + ψ_R(s) + L′/L(s, χ₅) on Re s = 2 |
| Lemma 1, Gamma side | ψ_R(s+1), poles at s = −1, −3, … | ψ_R(s), poles at s = 0, −2, −4, …, all outside 1/2 ≤ Re s ≤ 2. The pole of Γ_R(s) at s = 0 is cancelled by the trivial zero L(0, χ₅) = 0, so Λ is entire (also from the Mellin form) |
| Gauss's integral | z = 3/4 + it/2 | z = 1/4 + it/2: 2 Re ψ_R(1/2+it) = −log π − γ + ∫₀^∞ 2[e^{−2u} − e^{−u/2} cos tu]/(1 − e^{−2u}) du. The bound O(1 + min(1/u, t²u)) near u = 0 is unchanged |
| ramified prime | 3: P_3(0) = 0 | 5: P_5(0) = 0, Defect_5 = −min P_5 |
| decay of g_B | exp(−(2π/3)e^u) | exp(−(2π/5)e^u) (slower), so more primes enter and the tails are larger (section 2) |

Lemma 3 (cut F = b(i d/du)φ, not g) and Theorem 4 hold word for word with log 5 in place of log 3. The decay |F| + |F′| ≤ Ce^{−4|u|} holds because φ decays like exp(−(π/5)e^{2|u|}). Constants: 4 = 2² (from φ), c = 2π/5 (conductor 5 in θ₀), exponent 1/2 (from Γ_R(s)), d(N) (from χ*χ = χ·d).

## 2. Certificate and table (computer-assisted, [validate_exact5.py](validate_exact5.py), arb at 600 bits)

The method is that of exact3.md, section 5: Sturm count for B > 0, balls for g at prime powers n ≤ 400, branch and bound for min P_p (with the repaired branch `m > 0`), and one tail bound for everything above 400 (here below 10⁻¹⁶⁵). Search: [c5opt.py](c5opt.py) for degrees 0-10, [c5refine2.py](c5refine2.py) for degrees 12-18; certificates in cert5_d*.json.

| deg B | Defect/g(0), certified interval | q*_∞ ≥ 5 − … | largest share | b_k |
|---|---|---|---|---|
| 0 | [9.512, 9.513]·10⁻³ | 9.423·10⁻² | p = 5 (then 11) | + |
| 2 | [1.593, 1.594]·10⁻³ | 1.591·10⁻² | 7 | all + |
| 4 | [7.348, 7.349]·10⁻⁵ | 7.349·10⁻⁴ | 11 | all + |
| 6 | [8.372, 8.373]·10⁻⁶ | 8.373·10⁻⁵ | 13 | all + |
| 8 | [2.206, 2.207]·10⁻⁷ | 2.207·10⁻⁶ | 19 | all + |
| 10 | [8.227, 8.228]·10⁻⁹ | 8.228·10⁻⁸ | 23 | all + |
| 12 | [1.147, 1.148]·10⁻¹⁰ | 1.148·10⁻⁹ | 29 | all + |
| 14 | [1.918, 1.919]·10⁻¹³ | 1.919·10⁻¹² | 37 | all + |
| 16 | [2.763, 2.764]·10⁻¹⁵ | 2.764·10⁻¹⁴ | 41 | all + |
| 18 | [2.763, 2.764]·10⁻¹⁵ | 2.764·10⁻¹⁴ | 41 | all +, b_9/b_0 < 10⁻³⁰ |

(B = 1: q*_∞ ≥ 4.905776.) log₁₀(Defect/g(0)) = −2.02, −2.80, −4.13, −5.08, −6.66, −8.08, −9.94, −12.72, −14.56, −14.56.

**Structure (computed; the reading is heuristic).**
- p = 2 and 3 never carry cost: χ₅(2) = χ₅(3) = −1, g(log 2), g(log 3) > 0, and Defect_2, Defect_3 < 10⁻¹⁰⁰ g(0).
- The largest share outside p = 5 falls on the class χ₅(p) = (−1)^D at every degree except 18. This is the same parity pattern as for χ₋₃, from the eventual sign (−1)^D of g_B with b_D > 0.
- **Degrees 18 and 20 stall at the degree-16 value**, exactly as for χ₋₃. The degree-20 search returned b_10 = −6.6·10⁻³⁹ < 0 (B < 0 for large t; not admissible, not certified). The same stall at the same degrees for both characters makes an artefact of the search (start from the padded optimum, b_D ≈ 0) at least as likely as a property of the problem. This is not decided.

## 3. Checks (computed; [validate_exact5.py](validate_exact5.py), section A)

- θ₀ transformation at x = 0.6, 1.7, to 10⁻²⁵.
- Ξ(t) = ∫ φ cos at t = 0, 6.3 against mpmath's L(s, χ₅), to 10⁻²⁰.
- Closed form of φ*φ against quadrature at log 1, 2, 3, to 20 digits.
- g_k(0) for k ≤ 2 against the L-function side.
- K₀, K₁ against mpmath at 130 digits: 3.6·10⁻¹²⁵.
- Lemma 1 with g(0) and A₀ from the L-function side: residual −2.5·10⁻³¹ (B = 1, g(0) = 0.619) and 4.9·10⁻³² (B = 1 + t²/7, g(0) = 0.875).

Result: 137 checks, 0 failed ([validate_exact5_out.txt](validate_exact5_out.txt)). Not checked by a second code: the branch and bound and the tails (the χ₋₃ referee's check covered exact3 only).

## 4. Files

[validate_exact5.py](validate_exact5.py) (`quick` skips section A), [c5core.py](c5core.py), [c5basis.py](c5basis.py) (→ c5basis_D16.json), [c5opt.py](c5opt.py), [c5refine2.py](c5refine2.py), [c5mkcert.py](c5mkcert.py), cert5_d*.json, the logs c5refine2_d12_16.log and c5refine2_d18_20.log.
