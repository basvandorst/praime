# Bad primes in the Bogomolny–Keating two-point function, and the number variance of Dirichlet L-zeros

Files (all in this folder): `bklib.py` (formulas), `predict.py` → `predict_out.txt`, `predict.json`,
`predictions_per_lfunction.json`; `diag.py` → `diag_out.txt` (decomposition, noise); `validate.py` → `validate_out.txt`
(independent recomputation); `paper/Dirichlet_final.tex` (arXiv e-print of the paper, read in full).

Labels: **proved** = follows from the paper's formula by algebra/calculus; **heuristic** = an approximation or the
paper's own heuristic status; **computed** = numerical; **retrospective** = compared with values already known.

## 0. The source, checked

E. Bogomolny, J. P. Keating, "Two-point correlation function for Dirichlet L-functions", J. Phys. A 46 (2013) 095202,
arXiv:1307.6010. I read the TeX source (gzipped e-print). The brief's summary is correct. The definitions:

- χ a **primitive** character of modulus k; zeros 1/2 + iE_n; d(E) = Σ_n δ(E − E_n).
- Mean density d̄(E) = (1/2π) ln(kE/2π) (eq. `d_dirichlet`), from N(E) = (E/2π) ln(kE/2πe) + …
- R2(ε) = ⟨d(E − ε/2) d(E + ε/2)⟩ (written for ζ; the same for L). **ε is in raw height units**, and the average is
  over E in a window. The self term d̄ δ(ε) is not part of the formula (it describes ε ≠ 0).
- R2 = d̄² + R2diag + R2off (Conclusion section), with
  - R2diag(ε) = −(1/4π²) ∂²/∂ε² ln[ |ζ(1+iε)|² Φdiag(ε) Ψ(ε,k) ],
  - Φdiag(ε) = Π_p exp( 2 Σ_m (1−m) cos(mε ln p)/(m² p^m) ),
  - Ψ(ε,k) = exp( −Σ_{p|k} Σ_m e^{imε ln p}/(m² p^m) + c.c. ) (eq. `psi_diag_dirichlet`),
  - R2off(ε) = (1/4π²) e^{2πiε d̄} |ζ(1+iε)|² Φoff(ε) Ψoff(ε,k) + c.c. (eq. `R_off_dirichlet`;
    e^{2πiε d̄} = e^{iε ln(Ek/2π)}),
  - Φoff(ε) = Π_p ( 1 − (p^{iε} − 1)²/(p − 1)² ),
  - Ψoff(ε,k) = Π_{p|k} ( 1 + (p^{iε/2} − p^{−iε/2})²/(p − p^{−iε}) )^{−1} (eq. `psi_off_dirichlet`).
- Status in the paper: heuristic. It uses the Hardy–Littlewood conjecture for prime pairs in progressions (their
  α(h,k) = α(h) β(h,k)), swaps sums and averages without justification, and claims precision O(E^{−1/2}) ("optimistic").
- The character enters only through |χ(n)|² (diagonal) and the Gauss-sum identity Σ_r g(r) e^{2πirs/k} = |χ(s)|² |τ(χ)|²,
  with |τ(χ)|² = k (off-diagonal). Both hold for **every primitive χ, real or complex**. So the predicted R2 depends on
  χ only through its conductor k and the height E. **Proved** from the paper's derivation (section "Off-diagonal
  contribution").
- Sign convention, traced through the paper's derivation: the sum over R gives ζ(1+iε); the sum over q gives
  ζ(1−iε) Φoff(ε) × (factors for p | k). So e^{+iεL} goes with Φoff(ε) built from p^{+iε}. **Proved** (internal
  consistency), and checked on ζ's zeros in §7.1.

## 1. Unfolding: the changes of variable

- L := 2π d̄ = ln(kE/2π). This is the repo's L = ln(q t/2π), with q = k and t = E.
- r := d̄ ε, the distance in mean spacings, so ε = 2πr/L and dε = dr/d̄.
- Pair density per unit unfolded length: R(r) = R2(ε)/d̄². (R2 is a density in (E, E′); each dE = dr/d̄.)
- The phase: εL = (r/d̄)(2π d̄) = 2πr, so **e^{iεL} = e^{2πir}**, exactly period 1 in r.
- All other factors (ζ, Φ, Ψ, F_p) take the argument ε = 2πr/L.

**Proved** (a change of variables). The phase is exact in the repo's unfolded units, because
u = (t/2π) ln(qt/2πe) has du/dt = d̄.

## 2. Number variance from R

Let n(x) = #{j : x ≤ u_j < x + ℓ} in unfolded units (density 1). Then ⟨ρ(a)ρ(b)⟩ = δ(a − b) + R(a − b), so

  Var n = ∫_x^{x+ℓ}∫_x^{x+ℓ} [δ(a − b) + R(a − b) − 1] da db = ℓ + 2 ∫_0^ℓ (ℓ − r)(R(r) − 1) dr.

**Proved** for a stationary window. The measured Σ² (r123.py) averages over windows on a uniform grid in u, so over
heights weighted by the zero density ∝ L(t). Its prediction is therefore Σ²(ℓ) averaged over t ∈ [t1, t2] with weight
ln(qt/2π) (Gauss–Legendre, 8 nodes). **Heuristic** (local stationarity). In practice the average differs from the value
at the mean L by at most 1e-4 in the per-prime off-diagonal part (§7.4, diag_out.txt).

## 3. The GUE part and the regular remainder

As ε → 0, ζ(1+iε) = 1/(iε) + γ + O(ε). Hence (ζ'/ζ)'(1+iε) = −1/ε² + O(1) and |ζ(1+iε)|² = 1/ε² + O(1).
Also d²/dε² ln ζ(1+iε) = −(ζ'/ζ)'(1+iε). So:

- −(1/4π²) ∂² ln|ζ|² = (1/2π²) Re (ζ'/ζ)'(1+iε) ≈ −1/(2π² ε²);
- R2off ≈ (1/2π²) cos(εL)/ε² (since Φoff(0) = Ψoff(0) = 1);
- their sum divided by d̄² is −(4π²/L²)(1 − cos εL)/(2π² ε²) = −(sin πr/πr)².

Therefore

  R(r) − 1 = −(sin πr/πr)² + C(r),
  C(r) = (2/L²)[ Re(ζ'/ζ)'(1+iε) + 1/ε² + Re( e^{2πir} (|ζ(1+iε)|² Φoff Ψoff − 1/ε²) ) ]
         − (1/L²)[ ∂²_ε ln Φdiag + ∂²_ε ln Ψ ],

  Σ²(ℓ) = Σ²_GUE(ℓ) + 2 ∫_0^ℓ (ℓ − r) C(r) dr.

**Proved.** C(0) = 0 for every k, so R(0) = 0. Analytically: regd(0) + Z2m1(0) = 0 (both ±(γ² + 2γ₁) = ±0.18755),
(Φoff − 1)/ε² → Σ_p (ln p)²/(p−1)² = ½ ∂² ln Φdiag(0), and (Ψoff − 1)/ε² → Σ_{p|k}(ln p)²/(p−1) = ½ ∂² ln Ψ(0).
Numerically, |C(10⁻⁴)| < 1.4e-8 for k with bad sets ∅, {2}, {3,5}, {2,7,13}.

## 4. One prime p | k: the diagonal part (exact)

At fixed L, adding p to the prime divisors of k multiplies Ψ by Ψ_p = exp(−2 Σ_m cos(mε ln p)/(m² p^m)). Then:

1. ∂²_ε ln Ψ_p = 2 Σ_m (ln p)² cos(mε ln p)/p^m.
2. ΔR2diag(ε) = −(1/4π²) × (step 1) = −(1/2π²) Σ_m (ln p)² p^{−m} cos(mε ln p). These are exactly the terms
   n = p^m removed from (1/4π²) Σ_{(n,k)=1} Λ²(n)/n e^{iε ln n} + c.c.
3. ΔR_diag(r) = ΔR2diag/d̄² = −(2/L²) Σ_m (ln p)² p^{−m} cos(2π m r ln p/L).
4. ∫_0^ℓ (ℓ − r) cos(ar) dr = (1 − cos aℓ)/a² = 2 sin²(aℓ/2)/a², with a = 2πm ln p/L.
5. ΔΣ²_diag(ℓ) = 2 × (−2/L²) Σ_m (ln p)² p^{−m} × 2 sin²(πmℓ ln p/L) × L²/(4π² m² (ln p)²),
   so **ΔΣ²_diag(ℓ; p, L) = −(2/π²) Σ_m sin²(π m ℓ ln p/L)/(m² p^m)**.

This is the repo's (Berry's) formula exactly. **Proved.** Checks:
- the closed form against a quad integral of step 3: difference ≤ 2e-18 (validate_out.txt);
- our independent code reproduces r123.json's `pred` to 2.8e-17.

Remark (**proved**): for p | k, the p-factor of |ζ|², which is |1 − p^{−1−iε}|^{−2}, times Φdiag's p-factor times Ψ_p
equals 1. A bad prime drops out of the diagonal term completely.

## 5. One prime p | k: the off-diagonal part

**5.1 Algebra (proved).** Write F_p(ε) for the p-factor of Ψoff. Since (p^{iε/2} − p^{−iε/2})² = p^{iε} − 2 + p^{−iε},

  F_p(ε) = (p − p^{−iε})/(p − 2 + p^{iε}),  F_p(0) = 1,
  F_p(ε) − 1 = (2 − p^{iε} − p^{−iε})/(p − 2 + p^{iε}) = 4 sin²(ε ln p/2)/(p − 2 + p^{iε}).

To leading order in 1/p this is (4/p) sin²(ε ln p/2).

*Local-factor reading (our rewriting of the paper's product).* The p-factor of |ζ|², times Φoff's p-factor, times F_p
is (1 − 1/p)^{−2}, a constant. Put differently, in the off-diagonal amplitude e^{2πir}|ζ|²ΦoffΨoff:
- a good prime contributes A_p(ε) = (1 − 1/p)^{−2}/F_p(ε);
- a bad prime contributes the constant A_p(0) = (1 − 1/p)^{−2}.

The conductor freezes the prime's local factor at its ε = 0 value. In the diagonal term it removes the prime instead
(§4). One special case: A_3(ε) vanishes at ε = π/ln 3, where 3^{iε} = −1. This is the only zero of Φoff. When 3 | k it
is lifted, which shows in the p = 3 column near ℓ ≈ L/(2 ln 3) ≈ 4.3.

**5.2 The change in R and Σ² (proved).** Let S be the other primes dividing k. At fixed L,

  ΔR2off(ε) = (1/4π²) e^{iεL} |ζ(1+iε)|² Φoff(ε) Π_{q∈S}F_q(ε) (F_p(ε) − 1) + c.c.

Define H_p(ε; S) := |ζ(1+iε)|² Φoff(ε) Π_{q∈S}F_q(ε) (F_p(ε) − 1). It is regular at 0, with
H_p(0) = (ln p)²/(p−1), real and independent of S, because |ζ|² ≈ 1/ε² and F_p − 1 ≈ ε²(ln p)²/(p−1). Since
H(−ε) = conj H(ε), Re H is even and Im H is odd. Then

  ΔR_off(r) = (2/L²) Re[ e^{2πir} H_p(2πr/L; S) ],
  **ΔΣ²_off(ℓ; p, L, S) = (4/L²) ∫_0^ℓ (ℓ − r) Re[ e^{2πir} H_p(2πr/L; S) ] dr.**

**5.3 R(0) is unchanged (proved).** ΔR_diag(0) + ΔR_off(0) = −(2/L²)(ln p)²/(p−1) + (2/L²)(ln p)²/(p−1) = 0. A bad
prime does not change the repulsion at r = 0. The diagonal term alone would violate this.

**5.4 Shape (heuristic).** Integrate by parts twice, treating H as slowly varying on the scale L/(2π ln p); for
p = 2 that scale is only ≈ 2. With H(0) real and dH/dr(0) imaginary,

  ΔΣ²_off(ℓ) ≈ (1/(π² L²)) [ H_p(0) − Re( H_p(2πℓ/L) e^{2πiℓ} ) ] + (smaller terms from Im H and H″).

So:
- **period 1 in ℓ**;
- at half ℓ, ≈ (H(0) + Re H(ε_ℓ))/(π²L²) > 0 (too positive at half ℓ);
- at whole ℓ, ≈ (H(0) − Re H(ε_ℓ))/(π²L²). This is ≈ 0 while H(ε_ℓ) ≈ H(0), and turns negative once Re H(ε_ℓ) > H(0)
  (too negative at whole ℓ);
- the envelope |H_p(2πℓ/L)|/(π²L²) grows with ℓ for p = 2, since |F_2 − 1| = 4 sin²(ε ln 2/2) rises to 4 at
  ε = π/ln 2, that is ℓ = L/(2 ln 2) ≈ 6.8 at L = 9.38.

For small ℓ (ε ln p ≪ 1, H ≈ H(0)) the total is ΔR ≈ −(4H(0)/L²) sin²(πr), and

  ΔΣ²(ℓ) ≈ −(2 (ln p)²/((p−1) L²)) (ℓ² − sin²(πℓ)/π²), against the diagonal alone −(2 (ln p)²/((p−1) L²)) ℓ².

The ratio total/diagonal is 1 − 4/π² = 0.595 at ℓ = 1/2 and exactly 1 at ℓ = 1. The exact values at L = 9.38 are
0.53–0.56 at ℓ = 0.5 and 1.002–1.017 at ℓ = 1 (predict_out.txt). **Heuristic** (leading order), checked against the
exact numbers. This is why the per-prime law of rounds 117–122 held at ℓ = 1, where sin²(πℓ) = 0, and why it fails at
ℓ = 0.5.

**5.5 First order against exact.** Three levels:
- (a) leading 1/p: F_p − 1 → (4/p) sin²(ε ln p/2), with S = ∅;
- (b) single prime: exact F_p, S = ∅;
- (c) exact: S = the other prime divisors of k.

The off-diagonal part is **not additive** over bad primes; the diagonal part is (§4). Example, p = 2 at L = 9.0,
ℓ = 7: −0.0070 alone, −0.0083 with 29 | k, −0.0192 with 15 | k (validate_out.txt). In round 123's sample this
non-additivity, not the L-averaging, is what separates level (b) from level (c) (§7.4).

## 6. Numerics (computed)

- **Grid.** bklib.py tabulates the regular functions on ε ∈ [0, 7.6] with step 0.0025: ε²|ζ(1+iε)|²,
  (ε²|ζ|² − 1)/ε², Re(ζ'/ζ)'(1+iε) + 1/ε² (mpmath, 30 digits), ln Φoff (direct sum of log(1 − w_p) over p < 10⁶, with
  an accurate complex log1p) and ∂² ln Φdiag (analytic sum over p^m < 10¹⁸). Cubic splines in between.
- **Quadrature.** Gauss–Legendre, 8 nodes on each piece of width 0.025 (40 pieces per period of e^{2πir}), over
  r ∈ [0, 8]. This gives Σ² at ℓ = 0.025, 0.05, …, 8.
- **Checks** (predict_out.txt section [A], validate_out.txt):
  - spline against direct evaluation: ≤ 2e-11;
  - Φoff with p < 10⁵ against p < 10⁶: ≤ 3e-6;
  - one off-diagonal integral against scipy.quad per half-period: equal to 1e-8;
  - Σ²_GUE(1) = 0.344163 and Σ²_GUE(8) = 0.556673, as by quad (and as r129.json's `gue_exact`).
- **Independent recomputation** (validate.py): the paper's literal Ψoff form, mpmath ζ at every node with no splines,
  (ζ'/ζ)' by mpmath numerical differentiation of log ζ, ∂² ln Φdiag by finite differences, scipy.quad per half-period.
  It covers 8 single-prime cells, 5 cells with other bad primes, the multi-prime consistency
  Σ²(B) − Σ²(B∖p) − diag(p) = off(p | B∖p), 3 diagonal cells, and ζ's full Σ² at 4 lengths. Largest difference 1.5e-8.

## 7. Results

### 7.1 ζ (k = 1), retrospective: the normalisation and the sign

Round 129 measured 26,812 zeros with 7.6 ≤ ln(t/2π) ≤ 8.6 (r129.json). The full BK Σ²(ℓ), L-averaged over that
range, against the data and the alternatives:

| ℓ | data | BK | BK, flipped sign | Berry (r129) | GUE |
|---|---|---|---|---|---|
| 0.5 | 0.2741 | 0.2740 | 0.2773 | 0.2709 | 0.2802 |
| 1 | 0.3204 | 0.3205 | 0.3229 | 0.3196 | 0.3442 |
| 2 | 0.3715 | 0.3712 | 0.3742 | 0.3695 | 0.4157 |
| 4.5 | 0.4027 | 0.4030 | 0.4075 | 0.4013 | 0.4985 |
| 6 | 0.3903 | 0.3905 | 0.3971 | 0.3865 | 0.5275 |
| 8 | 0.3633 | 0.3640 | 0.3704 | 0.3621 | 0.5567 |

Over all 16 lengths:
- BK: max |data − BK| = 0.0008, mean −0.0002;
- Berry: max 0.0037, mean +0.0016;
- flipped sign (e^{−iεL} with Φoff(ε)): max 0.0072.

So the BK lower-order terms account for the systematic +0.002 that the fortieth review found against Berry's
leading-order formula. This fixes the normalisation (R = R2/d̄², e^{2πir}) and the sign convention on real data.
**Computed, retrospective.** BK say their ζ formula agrees with Odlyzko's zeros, so this confirms it in the repo's own
data rather than adding anything new.

### 7.2 Per-prime changes at L = 9.38 (the value registered in round 123), single prime (S = ∅)

Rows: diagonal (Berry) / off-diagonal (exact) / off-diagonal to leading order in 1/p. Columns: ℓ = 0.5, 1, 1.5, …, 8.

```
p=2  diag   -0.0027 -0.0099 -0.0202 -0.0321 -0.0443 -0.0563 -0.0675 -0.0775 -0.0861 -0.0932 -0.0986 -0.1023 -0.1042 -0.1042 -0.1025 -0.0989
     off    +0.0012 -0.0001 +0.0014 -0.0003 +0.0016 -0.0004 +0.0015 -0.0001 +0.0011 -0.0005 +0.0030 -0.0036 +0.0066 -0.0064 +0.0076 -0.0059
     off1/p +0.0005 -0.0000 +0.0005 +0.0001 +0.0003 +0.0003 +0.0001 +0.0004 +0.0004 +0.0000 +0.0000 +0.0016 -0.0022 +0.0033 -0.0025 +0.0023
p=3  diag   -0.0033 -0.0121 -0.0240 -0.0366 -0.0482 -0.0577 -0.0644 -0.0679 -0.0681 -0.0648 -0.0583 -0.0490 -0.0375 -0.0250 -0.0130 -0.0038
     off    +0.0015 -0.0002 +0.0018 -0.0005 +0.0021 -0.0008 +0.0023 -0.0001 +0.0001 +0.0027 -0.0019 +0.0033 -0.0012 +0.0019 +0.0002 +0.0009
p=5  diag   -0.0035 -0.0124 -0.0233 -0.0329 -0.0390 -0.0406 -0.0374 -0.0299 -0.0196 -0.0090 -0.0016 -0.0004 -0.0061 -0.0161 -0.0269 -0.0355
     off    +0.0016 -0.0002 +0.0018 -0.0001 +0.0013 +0.0007 +0.0005 +0.0008 +0.0008 +0.0007 +0.0008 +0.0008 +0.0006 +0.0014 -0.0004 +0.0023
p=7  diag   -0.0034 -0.0117 -0.0208 -0.0272 -0.0289 -0.0254 -0.0177 -0.0084 -0.0014 -0.0005 -0.0060 -0.0151 -0.0236 -0.0284 -0.0282 -0.0230
     off    +0.0015 -0.0001 +0.0016 +0.0001 +0.0010 +0.0007 +0.0007 +0.0007 +0.0007 +0.0007 +0.0006 +0.0011 -0.0002 +0.0022 -0.0006 +0.0015
p=11 diag   -0.0030 -0.0100 -0.0163 -0.0184 -0.0154 -0.0087 -0.0021 -0.0001 -0.0041 -0.0112 -0.0170 -0.0182 -0.0144 -0.0074 -0.0013 -0.0004
     off    +0.0014 -0.0000 +0.0013 +0.0003 +0.0008 +0.0006 +0.0007 +0.0007 +0.0006 +0.0007 +0.0006 +0.0010 +0.0002 +0.0009 +0.0006 +0.0006
p=13 diag   -0.0029 -0.0092 -0.0145 -0.0153 -0.0112 -0.0047 -0.0003 -0.0014 -0.0071 -0.0132 -0.0156 -0.0129 -0.0068 -0.0012 -0.0004 -0.0051
     off    +0.0013 -0.0000 +0.0011 +0.0003 +0.0008 +0.0006 +0.0006 +0.0006 +0.0006 +0.0006 +0.0006 +0.0008 +0.0005 +0.0006 +0.0007 +0.0006
```

The same tables at L = 9.2246 (round 121's mean) and L = 10.6859 (round 122's mean) are in predict_out.txt [C].

### 7.3 The shape on the fine grid (step 0.025, L = 9.38), computed

**Period.** Extrema of the off-diagonal part sit about 1 apart. For p = 2:
- maxima at 0.525, 1.500, 2.475, 3.425, 4.700, 5.650, 6.550, 7.475;
- minima at 1.025, 1.975, 2.950, 3.925, 5.175, 6.100, 7.000, 7.950.

For larger p the spacing is slightly below 1 (0.875–0.975 for p = 3), because H's phase turns slowly.

**Signs.** The off-diagonal part is positive at every half ℓ for p = 2 (+0.0011 to +0.0076). At whole ℓ it is
−0.0001 to −0.0005 up to ℓ = 5, then −0.0036, −0.0064, −0.0059 at ℓ = 6, 7, 8.

**Growth.** The swing "half minus the mean of its two whole neighbours" for p = 2:
+0.0043 (at ℓ = 0.5), +0.0017, +0.0019, +0.0017, +0.0014, +0.0050, +0.0116, +0.0137 (at 7.5). It grows strongly for
ℓ ≳ 5.5. For p = 3 the swing reverses sign beyond ℓ ≈ 4.5: −0.0012, −0.0050, −0.0038, −0.0012. That is the lifted
zero of A_3 from §5.1. For p ≥ 7 the off part is a small positive offset, (ln p)²/((p−1)π²L²) ≈ 0.0006–0.0007, with
weak oscillation.

**ℓ = 0.5.** total/diagonal = 0.53–0.56 for every prime. **At ℓ = 1:** 1.00–1.02.

### 7.4 Retrospective comparison with round 123 (r123.json; post-repair fits, 181 real L-functions, 100 ≤ γ ≤ 600)

The fitted values are known, so this comparison is **retrospective**. No parameter was fitted; the predictions use only
the discriminants, the height range and the paper's formula. The χ² uses r123.json's OLS standard errors.

| prediction | χ² (96 cells) | per prime 2 / 3 / 5 / 7 / 11 / 13 | slope fit/pred |
|---|---|---|---|
| diagonal only (r123 `pred`) | 707.5 | 425.6 / 137.9 / 63.2 / 27.9 / 24.3 / 28.6 | 0.991 |
| diagonal only, regression-matched, L-averaged | 684.1 | 431.9 / 134.9 / 60.9 / 21.9 / 14.0 / 20.6 | 0.992 |
| diagonal + off, single prime (S = ∅), group mean at mean L | 308.9 | 114.7 / 112.3 / 43.1 / 12.3 / 13.4 / 13.2 | 1.002 |
| **diagonal + off, exact (all p \| D), regression-matched, L-averaged** | **15.7** | **2.0 / 2.8 / 1.8 / 2.3 / 1.1 / 5.8** | **0.999** |

χ² per ℓ (6 cells each), ℓ = 0.5 … 8:
- diagonal: 48 6 17 20 81 33 56 25 43 24 54 38 63 64 92 43;
- full: 3 1 1 2 3 1 1 0 1 0 0 0 1 0 0 1.

"Regression-matched" means: predict Σ²_i(ℓ) for each L-function i, with all its bad primes and averaged over its height
range, then apply round 123's OLS (constant plus indicators of 2, 3, 5, 7, 11, 13) to the predicted values. The base
Σ²(ℓ) without bad primes is in the predicted values but cancels except through L-variation.

- **The residual of the diagonal law against the predicted off-diagonal part** (96 cells): slope 0.995 (se-weighted
  0.973), correlation 0.942.
- **Review example** (p = 2): fit − diagonal = −0.0101 (se 0.0013) at ℓ = 7 and +0.0108 (se 0.0012) at ℓ = 7.5. The
  exact off-diagonal part predicts −0.0098 and +0.0107. The single-prime version predicts only −0.0063 and +0.0075.
  diag.py decomposes the difference: L-averaging moves it by 0.0001; the other primes dividing D move it to −0.0099
  and +0.0106.
- **ℓ = 0.5.**

  | | p = 2 | 3 | 5 | 7 | 11 | 13 |
  |---|---|---|---|---|---|---|
  | fits | −0.0015 | −0.0024 | −0.0021 | −0.0019 | −0.0014 | −0.0006 |
  | se | 0.0004 | 0.0004 | 0.0005 | 0.0005 | 0.0007 | 0.0007 |
  | diagonal | −0.0027 | −0.0033 | −0.0035 | −0.0034 | −0.0030 | −0.0029 |
  | residual vs diagonal | +0.0011 | +0.0010 | +0.0014 | +0.0015 | +0.0016 | +0.0024 |
  | full | −0.0016 | −0.0018 | −0.0020 | −0.0018 | −0.0015 | −0.0014 |
  | residual vs full | +0.0000 | −0.0006 | −0.0001 | −0.0001 | +0.0001 | +0.0008 |

  Slope against the diagonal 0.53, against the full law 1.01. The predicted ratio full/diagonal at ℓ = 0.5 for this
  design is 0.557; the measured ratio is 0.566.
- **Why χ² < 96 for the full law.** The OLS standard errors come from the residual SD s of the indicator regression.
  That SD includes between-L-function variation that BK predicts: L-dependence and primes above 13.
  - At ℓ = 1: s = 0.0046, but the per-D rms of measured − BK is 0.0028. That is a variance ratio of 2.7.
  - The 16 lengths are strongly correlated (the review estimated about 2.5 independent ones).

  So these χ² values are not calibrated goodness-of-fit numbers. The meaningful comparisons are the ratio 684/15.7 and
  the slope. The fresh test therefore uses bootstrap errors.

### 7.5 Absolute level at ℓ = 1 (no fitted constant), retrospective

| sample | mean measured Σ²(1) | mean BK | per-L-function |
|---|---|---|---|
| round 120 (181 real, 100–600) | 0.3085 | 0.3089 | correlation 0.973, rms difference 0.0028 |
| round 121 (88 complex, 100–500) | 0.3123 | 0.3120 | |
| round 122 (61 real, 900–1500) | 0.3165 | 0.3165 | |

β at ℓ = 1 (the regression on Berry's x):

| sample | measured | BK-predicted values | diagonal-predicted values |
|---|---|---|---|
| round 121 | 1.005 | 1.052 | 1.031 |
| round 122 | 0.944 ± 0.042 | 1.027 | 0.999 |

At ℓ = 1 the off-diagonal term changes little, as §5.4 says. Round 122's β is 2σ below BK, which is not new.

## 8. Doubts

1. **The paper's formula is heuristic.** It rests on Hardy–Littlewood for pairs in progressions, unjustified swaps of
   sums, and "≈" steps with h ≪ n, at claimed precision O(E^{−1/2}). The heights here are E = 100–600 (E^{−1/2} up to
   0.1 in some unit the paper does not fix). The absolute agreement at ℓ = 1 (≤ 0.0004 in three samples) and for ζ at
   all ℓ (≤ 0.0008) says the finite-height error is small here. That is evidence, not proof. For L-functions at
   ℓ > 1 the absolute level has not been checked (no stored data). The fresh test registers it as a secondary reading.
2. **Normalisation**: R = R2/d̄², and e^{iεL} = e^{2πir} in the repo's unfolding. This is checked three ways: the GUE
   limit, R(0) = 0 for every k, and ζ's data to 0.0008.
3. **Sign convention**: proved internally (§0) and preferred by ζ's data: max error 0.0008, against 0.0072 for the
   flipped convention. The two differ by 0.002–0.007 at every ℓ, while the data errors are ≤ 0.0008.
4. **Real against complex characters.** The formula does not distinguish them (§0). The repo's exploratory remark
   (STUDY round 121) that real characters might be more rigid because χ² is principal has no counterpart in the paper's
   single-L-function pair correlation. Twisted χ² terms belong to family averages (low-lying zeros), not to one
   L-function high up. Test A checks this.
5. **Complex characters at ℓ = 0.5.** The thirty-eighth review found a slope of 0.30 (fit against the diagonal law)
   before the data repair; BK predicts 0.55 for that sample. That is about 2–3σ, and the review's data still missed
   close pairs on 15 of 176 sides. Missing pairs bias exactly the small-ℓ statistics. Test A's reading (3) decides it.
6. **Non-additivity.** Per-prime regression coefficients depend on the design, that is on which other primes divide the
   conductors in the sample. Predictions must be regression-matched to the sample. The draft fixes this.
7. **Zero completeness** (missed close pairs) biases small-ℓ Σ². The draft requires the repo's completeness check.
8. **Finite-sample bias of np.var** on overlapping windows, and unfolding without the O(1) and O(1/t) terms of N(t),
   are not modelled. They affect mainly the constant, so they matter only for the absolute reading.
9. **Primitive characters only.** All samples are primitive.
10. **No cutoff.** The diagonal sum runs over all prime powers, as in the repo's formula. BK has no τ = 1 cutoff; the
    off-diagonal term does what Berry's saturation term did.
11. I did **not** search the citing literature. Whether someone tested Ψ and Ψoff numerically on actual zeros is
    unknown to me. The paper itself contains no numerics for L-functions.

## 9. What is whose

- **Bogomolny–Keating (2013)**: the two-point function with Ψ and Ψoff for primitive Dirichlet L-functions (heuristic).
  Its ζ version is Bogomolny–Keating 1996. The number variance as a double integral of R2 is standard. The diagonal
  per-prime term is Berry 1988 (as the repo already says).
- **Derived here** (simple consequences of the paper's formula, not found in the paper):
  - the per-prime off-diagonal number-variance term ΔΣ²_off(ℓ; p, L, S) and its evaluation;
  - the local-factor reading (a bad prime freezes its off-diagonal Euler factor at (1 − 1/p)^{−2} and drops out of the
    diagonal);
  - the small-ℓ form −(2(ln p)²/((p−1)L²))(ℓ² − sin²(πℓ)/π²), hence the ≈ 0.55–0.6 shortfall at ℓ = 0.5 and the exact
    agreement at ℓ = 1;
  - the non-additivity across bad primes.

  Whether any of these is in print elsewhere was not searched.
- **The repo's own** (as far as found, pending a search of citing papers): a numerical test of BK's bad-prime factors
  on actual zeros.
  - Retrospectively: round 123's 96 per-prime coefficients go from χ² 684 (diagonal) to 15.7 (BK), with slope 0.999.
  - The absolute Σ²(1) within 0.0004 in three samples with no fitted parameter.
  - Prospectively: the registration in round_draft.txt.
- **Not new**: that BK's ζ formula fits ζ's zeros (BK cite Odlyzko's data). Here it closes round 129's 0.002 residual.
