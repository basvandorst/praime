# The variance of log|L′(ρ, χ)| over the zeros, and what each prime dividing the conductor does to it

Folder: this directory (`varlp/`). Nothing in the repo was changed; the repo's `tplib.py`, `bklib.py` and `amp/grid.npz`
were copied into `lib/` (sha256 bd049701…, 9bfc1a3b…, 7af1eb14…) and used as libraries.

Labels: **proved** (algebra/calculus from stated formulas, or an identity checked to machine precision), **conjecture**
(the ratios conjecture / Bogomolny–Keating correlation functions with the bad primes deleted — the law registered in rounds
134–138), **computed** (numbers from the scripts here), **heuristic**, **retrospective** (a comparison with data already
seen), **literature**.

| script | output | what it does |
|---|---|---|
| `cue_check.py` | `cue_check_out.txt` | §2: the Palm formula against the exact CUE(N) variance, N = 2, 10, 20, 50 |
| `tab.py` | `tab_out.txt`, `lib/xgrid.npz` | ζ(1+iε), ζ′/ζ, Φ tables extended from ε ≤ 7.6 to ε ≤ 160 (same definitions as `bklib.Grid`, `gtab.py`) |
| `vlib.py`, `check_vlib.py` | `check_vlib_out.txt` | T₂, T₃ of the registered law on the long range; equal to bklib/tplib to 2·10⁻⁶ |
| `oned.py`, `oned2.py`, `mu.py` | `oned_out.txt`, `oned2.json`, `mu_out.txt`, `mu.json` | the one-dimensional pieces (§3) |
| `twod.py`, `twod2.py` | `t_*_out.txt`, `runs/t2_L*_V70*` | the three-point piece ΔV₃ (§4), quadrature checks |
| `assemble.py` | `assemble_out.txt` | ΔVar(S, L) from the pieces |
| `retro.py` | `retro_out.txt`, `retro.json` | §6: retrospective comparison with round 142 |
| `fresh.py` | `predictions_varlp.json` (+ sha256) | §7: frozen predictions for the fresh sets C, Bp |
| `var_test.py` | | the registered test (`round_draft.txt`) |
| `fresh_matched.py` | `fresh_matched_out.txt`, `fresh_matched.json` | matched frozen coefficients, expected se, power |
| `retro.py` (also) | `retro_pred_r142.json`, `var_test_retro_out.txt` | `var_test.py` run end to end on round 142's stored values (`--band 100 600`) |
| `farfield.py` | `farfield_out.txt` | a pointwise look at the far field of ΔT₃ (see §4.3: not a clean test) |

## 0. In short

1. **Exact representation** (proved, given the zeros on the line). log|L′(ρ_n)| − ln(L/2π) = F_n + h(γ_n), where
   F_n = Σ_{j≠n} ln|r_j| − ∫ ln|r| dr (r_j = unfolded distances of the other zeros from γ_n) and h is a smooth deterministic
   function that does not depend on which primes divide q. So the variance over the zeros at fixed L is Var F under the
   Palm measure (a zero at the origin).
2. **Palm formula** (proved): Var F = A + B + V₃ − μ², with f = ln|u|, T₂, T₃ the connected (unit-density) two- and three-point
   functions:
   A = ∫ f² T₂,  B = ∫ f² + ∫∫ f(u) f(v) T₂(u−v) = ∫ |f̂(τ)|² b(τ) dτ,  V₃ = ∫∫ f(u) f(v) T₃(0,u,v),  μ = ∫ f T₂.
   Checked on CUE(N) against Hughes–Keating–O'Connell's exact variance: agreement to 10⁻¹⁵ at N = 2, 10, 20, 50 (exact
   Fourier sums and 2-D quadrature).
3. **Per bad prime** (the law: p deleted from every prime sum/product of T₂, T₃): ΔVar = ΔB + ΔA + ΔV₃ − (2μΔμ + Δμ²) with
   closed forms ΔB_diag = −½Li₂(1/p), ΔA_diag = −(2/L) Σ_m ln p (γ_E + ln(2πm ln p/L))/(m p^m), ΔV₃,Q = (ln p/L) ln²(p/(p−1)),
   Δμ_diag = ln p ln(p/(p−1))/L, μ = γ_E + ln 2π − 1 + (γ_E + C₁)/L + μ_off(L) (C₁ = Σ_p Σ_{m≥2}(1−1/m) ln p/p^m = 0.434480),
   and the far field ΔT₃(0,u,v) → −ΔT₂,diag(v)·sin(2πu)/(πu), whose integral is known exactly.
4. **First order** (derived): ΔVar = −½Li₂(1/p) + (ln p/L)·ln(p/(p−1))·[2 + ln(p/(p−1))] + O(1/L²). This is exactly the
   first-order term of round 142's heuristic (LDM: 1.8667, 0.9753, 0.4961 per unit ln p/L for p = 2, 3, 5). The heuristic's
   26 % miss is therefore entirely second order (and beyond); it is supplied by the exact three-point integrals.
5. **Finite L** (computed): at L = 10.2, ΔVar = −0.191, −0.104, −0.050, −0.031, −0.017, −0.013 for p = 2, 3, 5, 7, 11, 13
   (leading order −0.291, −0.183, −0.106, −0.074, −0.047, −0.039; LDM −0.164, −0.078, −0.027, −0.011, +0.000, +0.003).
   Retrospective, matched to round 142's design: see §6 (the drops are reproduced to ≈ 1–2 se).
6. **Fresh test** (§7): frozen per-L-function predictions for C and Bp at heights 600–1000 (`predictions_varlp.json`),
   registration draft `round_draft.txt`.

(The numbers in items 5–6 and §6–7 are copied from the outputs named there.)

## 1. log|L′(ρ)| through the zeros (proved)

χ = χ_D real primitive, q = |D|, a = 0 (D > 0) or 1 (D < 0). Λ(s) = (q/π)^{(s+a)/2} Γ((s+a)/2) L(s, χ) is entire of order 1
with Λ(1−s) = Λ(s) (root number 1 for quadratic characters) and Λ(s̄) = conj Λ(s); so Ξ(t) := Λ(1/2 + it) is real and even,
and, under GRH, all its zeros are real: Ξ(t) = Ξ(0) Π_{γ>0} (1 − t²/γ²) (a zero at t = 0 adds a factor t^{2k}; none occurs in
the sets used). Let N(x) count the zeros of Ξ in [0, x] and N(−x) = −N(x) (the zeros are symmetric), and N̄ the smooth count,
N̄′(x) = (1/2π)[ln(q/π) + Re ψ((1/2 + a + ix)/2)] (it depends on q only through ln q). Taking logarithms of the product and of the
same product written for the smooth density (Stirling's formula for Γ is the Hadamard product of the Gamma factor, up to
smooth terms), for t not a zero

    ln|L(1/2 + it, χ)| = ∫ ln|t − x| d(N − N̄)(x) + h(t),                                                    (1.1)

the integral over the whole line (the symmetric pairing makes it converge; it is Abel-summable at ∞ because N − N̄ = S + O(1/x)
has mean 0), h smooth and O(1/t) deterministic, independent of the bad primes. At a simple zero γ_n,
ln|L′(ρ_n)| = lim_{t→γ_n} [ln|L(1/2+it)| − ln|t − γ_n|], so

    Y_n := ln|L′(ρ_n)| = ∫ ln|γ_n − x| d(N − N̄ − δ_{γ_n})(x) + h(γ_n).                                      (1.2)

Unfold around γ_n: x − γ_n = 2πr/L, L = ln(qγ_n/2π) (the mean density is L/2π + O(1/γ)); ln|γ_n − x| = ln|r| + ln(2π/L). The
measure d(N − N̄ − δ) has total mass −1 (Abel), so Y_n = ln(L/2π) + F_n + h̃(γ_n) with

    F_n = Σ_{j≠n} ln|r_j| − ∫ ln|r| dr,                                                                       (1.3)

where h̃ collects h and the (deterministic) deviation of the unfolding from linearity beyond |r| ~ L, both smooth in γ_n and
independent of which primes divide q. The round-142 statistic is the sample variance over the zeros in a band of
y = Y − ln L, i.e. (band average of) Var F at the local L, plus the band variance of the deterministic mean, which is common
to all L-functions to 10⁻⁵ (it changes by 2 Cov_band(μ, Δμ) ≈ 10⁻⁵ per bad prime: neglected).

## 2. The Palm formula; the check on CUE (proved; computed)

Condition on a zero at the origin. The other zeros then have one- and two-point densities R₂(0,u) and R₃(0,u,v) (unit
density; R₂ = 1 + T₂, R₃ = 1 + T₂(u) + T₂(v) + T₂(u−v) + T₃). With f(u) = ln|u|,
E F = ∫ f T₂ =: μ and E F² = ∫ f² R₂ + ∫∫ f f R₃ − 2∫f ∫f R₂ + (∫f)², hence

    Var F = ∫ f²(1 + T₂) + ∫∫ f(u) f(v) [T₂(u − v) + T₃(0, u, v) − T₂(u) T₂(v)]
          = A + B + V₃ − μ²,                                                                                  (2.1)
    A = ∫ f² T₂,   B = ∫ f² + ∫∫ f(u) f(v) T₂(u − v) = ∫ |f̂(τ)|² b(τ) dτ = ½ ∫₀^∞ b(τ)/τ² dτ,   V₃ = ∫∫ f f T₃(0,u,v),

using ∫∫ f(u) f(v) T₂(u) T₂(v) = μ². b = 1 + T̂₂ is the form factor and f̂(τ) = −1/(2|τ|) (τ ≠ 0; f̂ = ∫ f e^{−2πiτu}).
**Where ½ ln L comes from:** B is the only divergent-looking piece. Its diagonal part b(τ) = Σ_n Λ(n)²/(nL²) δ(τ − ln n/L)
(τ < 1) gives ½ Σ_{n<e^L} Λ(n)²/(n ln²n) = ½ Σ_p Σ_m 1/(m² p^m) ≈ ½ ln L + const; the GUE form factor |τ| would give a log
divergence at τ → 0, cut off by the primes (no orbit shorter than ln 2). The bad prime removes its own terms: §3.1.

**CUE check** (`cue_check.py`; the same construction is exact there): eigenangles of CUE(N), K(x) = (1/2π)Σ_{m<N} e^{imx},
f(u) = ln|2 sin(u/2)|, R₂, R₃ the 2×2, 3×3 determinants; ln|Λ′(θ_n)| = Σ_{j≠n} f(θ_j − θ_n). With f_k = −1/(2|k|):
Var = (N−1)π²/12 − (1/N)Σ_{k<N}(N−k)(H_{k−1}+H_k)/k − ½Σ_{k<N}(N−k)/k² + (1/2N)Σ_{b<N}(H_b+H_{N−1−b})² − (H_N − 1)²
against the exact K″(0) = Σ_{j=1}^{N−1}[ψ′(j+2) − ½ψ′(j+1)]:

| N | exact | Palm, Fourier sums | Palm, 2-D quadrature (graded Gauss panels) | mean H_N − 1 (quadrature) |
|---|---|---|---|---|
| 2 | 0.072467033424 | diff −1.7·10⁻¹⁶ | | |
| 10 | 0.568081040802 | diff 1.7·10⁻¹⁵ | diff −4.1·10⁻¹⁵ | 1.92896825 (1.92896825) |
| 20 | 0.867947780587 | diff −7.8·10⁻¹⁶ | diff 9.0·10⁻¹⁵ | 2.59773966 (2.59773966) |
| 50 | 1.297036232786 | diff 3.6·10⁻¹⁵ | diff −8.9·10⁻¹⁵ | 3.49920534 (3.49920534) |

So (2.1), its normalisation (densities divided by ρ₁ = N/2π, i.e. d̄) and the quadrature scheme with log singularities are
right. (The exact CUE variance is ½ ln N − 0.678785 + 0.956/N … ; at N = 10 the 1/N term is +0.096.)
GUE limit (computed and proved): μ_GUE = −∫ ln|u| sinc²u du = γ_E + ln 2π − 1 = 1.415093 (`oned_out.txt`: numerical
1.415093).

## 3. The bad prime: two-point pieces (conjecture-level input; closed forms proved; numbers computed)

The law (registered, rounds 134–138; `bklib`, `tplib`): adding p to the bad set S deletes p from every prime sum and
product. For two points (bklib `C_of_r`, `dsig_off`), with ε = 2πr/L,

    ΔT₂(r) = ΔT₂,diag + ΔT₂,off,   ΔT₂,diag = −(2/L²) Σ_m ln²p p^{−m} cos(m ε ln p),
    ΔT₂,off = (2/L²) Re[e^{2πir} ε²|ζ(1+iε)|² Φ(ε) Ψ_S(ε) (F_p(ε) − 1)/ε²],   F_p = (p − p^{−iε})/(p − 2 + p^{iε}).

With ∫ ln|r| cos(ωr) dr = −π/|ω|, ∫ ln²|r| cos(ωr) dr = 2π(γ_E + ln|ω|)/|ω| (from ∫₀^∞ r^{s−1}cos ωr dr = Γ(s)cos(πs/2)ω^{−s},
second derivative at s = 1) and ∫∫ ln|u| ln|v| cos(ω(u−v)) du dv = π²/ω² (i.e. |f̂(τ)|² = 1/(4τ²), ω = 2πτ):

    ΔB_diag = −½ Li₂(1/p)              (all L: the leading-order law),                                         (3.1)
    Δμ_diag = ln p · ln(p/(p−1))/L      (round 142's first-order mean shift),                                    (3.2)
    ΔA_diag = −(2/L) Σ_m [ln p/(m p^m)] (γ_E + ln(2π m ln p/L))   (contains ln L/L; cancelled by ΔV₃, §4),        (3.3)
    ΔB_off = −(π²/2) ∫ |r| ΔT₂,off dr  (uses the sum rule ∫ ΔT₂,off = 0, checked: |·| < 10⁻⁵),
    ΔA_off = ∫ ln²|r| ΔT₂,off,   Δμ_off = ∫ ln|r| ΔT₂,off   (numerical; Gaussian damping e^{−(ε/E0)²}, E0 = 15, 30, 45:
                                                             equal to < 10⁻⁵ for p ≤ 97 at L ≥ 7, `oned2.json`).

The closed forms (3.1)–(3.3) agree with the damped numerical integrals of ΔT₂,diag, converging like 1/E0² (`oned_out.txt`).
Δμ = Δμ_diag + Δμ_off at L = 10: +0.04301 (p = 2), +0.03831 (3), +0.02932 (5), +0.02362 (7) — the review's mean-shift
numbers with the swapped term (+0.0428, +0.0374, +0.0296, +0.0236 band- and design-averaged; §6 matches them).

**The base mean μ** (enters through −2μΔμ). μ = ∫ ln|r| T₂ = μ_GUE + ∫ ln|r| C(r) with C = T₂ − T₂,GUE. The diagonal part of
C has the Fourier transform Σ_n Λ(n)²/(nL²)[δ(τ−τ_n) + δ(τ+τ_n)] − |τ| (the −|τ| is the transform of the 1/ε² in bklib's
`regd`), so its ln-integral is −(1/L)[Σ_{n≤x} Λ(n)²/(n ln n) − ln x] → −c₀/L with c₀ = −γ_E − C₁ (Mertens:
Σ_{p≤x} ln p/p = ln x − γ_E − Σ_p ln p/(p(p−1)) + o(1)). Hence (proved given the law; `mu.py`)

    μ(L) = γ_E + ln 2π − 1 + (γ_E + C₁)/L + μ_off(L),   C₁ = 0.434480,                                       (3.4)

μ_off = ∫ ln|r| (2/L²) Re[e^{2πir} rego] (numerical, converged to 10⁻⁶): μ = 1.5328, 1.5025, 1.4896 at L = 7, 10, 12.
The direct damped integral of C converges slowly (like ln E0/E0) towards the same value (0.0618, 0.0724, 0.0764 → 0.0874
at L = 10).

## 4. The bad prime: the three-point piece

### 4.1 The diagonal (triple-orbit) term: closed form (proved)

The Q-part of T₃ (tplib; Conrey–Snaith's Q(x, y) = −Σ_p ln³p/((p^{1+x}−1)(p^{1+y}−1))) loses its p-summand:
ΔT₃,Q(0,u,v) = (2 ln³p/L³) Σ_{m,n≥1} p^{−m−n} Re Σ_{apex} e^{−iκ(m(b−a) + n(c−a))}, κ = 2π ln p/L. Each exponential integrates
against ln|u| ln|v| to π²/(κ² · two frequencies), and the three apexes give Σ p^{−m−n}[1/(mn) + 2/(n(m+n))] = 2 ln²(1−1/p):

    ΔV₃,Q = (ln p/L) · ln²(p/(p−1)).                                                                         (4.1)

(The 2-D quadrature reproduces it: 0.03318 vs 0.03330 for p = 2, L = 10, 0.01806 vs 0.01806 for p = 3, `runs/`.) This is the
same number as −E[A_p² ε_p] of the random-phase model (A_p = −ln|1 − χ(p)p^{−1/2−iγ}|, ε_p the Landau–Gonek density
modulation): the zero-weighting of the Euler factor.

### 4.2 The far field and the first order (derived; checked by L-scaling)

When one zero (v) is far (|v| ≫ 1) from the pair (0, u), the change of T₃ is governed by the long-wave density modulation
that the prime carries. In the Bogomolny–Keating picture (zeros = sine process in the coordinate N_Θ(x) that includes the
short orbits), expanding R₃(0,u,v) = E_Θ[R₂^Θ(0,u) R₁^Θ(v)] to first order in the modulation ε(x) (E ε(0)ε(v) = C(v),
C = T₂,diag) gives

    T₃(0, u, v) → −C(v) D(u),   D(u) = sinc²u + (u sinc²u)′ = sin(2πu)/(πu)   (D̂ = 1 on |τ| < 1),               (4.2)

consistent with the exact coincidence rule T₃(0,0,v) = −2T₂(v) (D(0) = 2). With the two other orientations (the far zero is u,
or the origin is far from the pair (u,v)) the model M(u,v) = −ΔC(v)D(u) − ΔC(u)D(v) − ΔC(u)D(v−u) has the exact integral

    ∫∫ ln|u| ln|v| M = 2(γ_E + ln 2π) Δμ_diag − ΔA_diag + J,   J = ∫ ln|u| Ci(2π|u|) ΔC(u) du = O(1/L²),       (4.3)

using ∫ ln|u| D = −(γ_E + ln 2π) and ∫ ln|v| D(v−u) dv = ln|u| − Ci(2π|u|). Adding (3.1)–(4.3) and −2μ_GUE Δμ:
the ln L/L terms of ΔA cancel, 2(γ_E + ln 2π) − 2μ_GUE = 2, and

    ΔVar = −½ Li₂(1/p) + (ln p/L) ln(p/(p−1)) [2 + ln(p/(p−1))] + O(1/L²).                                   (4.4)

**Check that (4.2) is the exact first order of the ratios T₃** (computed): the residual ΔV₃,P − ∫∫M from the full 2-D
integral scales like 1/L²: for p = 2 it is −0.0906, −0.0545, −0.0376 at L = 7.5, 10, 12, i.e. L²·residual = −5.10, −5.45,
−5.41. (A pointwise look, `farfield.py`, is not a clean test: ΔT₃,off(v) carries the non-decaying off-diagonal factor
e^{2πiv}|ζ(1+iε)|²Φ(F_p − 1), of size 0.01–0.1 for p = 2, which period-averaging does not remove because its envelope has
all the frequencies ln(m/n). It integrates to ≈ 0 against the slowly varying ln|v|, which is what (4.2) needs.)

(4.4) equals round 142's LDM first-order term (its coefficients −1.8667, −0.9753, −0.4961 per unit ln p/L are exactly
ln(p/(p−1))(2 + ln(p/(p−1)))). So the LDM heuristic was right at first order; its failure is in the higher orders.

### 4.3 Finite L: the 2-D integral (computed)

ΔV₃,P = ∫∫ ln|u| ln|v| ΔT₃,P(0,u,v), ΔT₃,P = (P-part of T₃ with S = {p}) − (S = {}) from tplib's formulas (`vlib.T3P`, equal to
`tplib.T3_parts` to 2·10⁻⁶), over the half plane u ≥ 0 (T₃ even), |u|, |v| ≤ V; Gauss panels of width 1 (6 nodes), graded
at u = 0 and v = 0 (log singularities); coincident nodes u = v shifted by 10⁻⁶ (T₃ continuous). The non-decaying far field
is handled by subtracting M (its exact integral (4.3) added back) and a product damping d(u)d(v) of widths V/5 … V/2.7, then
Richardson in 1/E0². Convergence (L = 10, p = 2): 6 vs 8 nodes 10⁻⁵; primes in B1 ≤ 500 vs ≤ 1500: 10⁻⁵; V = 30, 60, 70:
0.17134, 0.17239, 0.17225; spread over the four dampings at V = 70: 1·10⁻⁴. Values (V = 70):

| L | p = 2: ΔV₃,P (∫∫M; L²·resid) | p = 3: ΔV₃,P (∫∫M; L²·resid) | p = 5: ΔV₃,P (∫∫M; L²·resid) | p = 7: ΔV₃,P (∫∫M; L²·resid) | pair (2,3) non-additive, damping range |
|---|---|---|---|---|---|
| 7.5 | 0.2464 (0.3370; -5.09) | 0.2361 (0.3517; -6.50) | 0.1952 (0.3108; -6.50) | 0.1598 (0.2711; -6.26) | +0.0100 … +0.0318 |
| 8.5 | 0.2096 (0.2841; -5.39) | 0.2099 (0.2984; -6.39) | 0.1751 (0.2650; -6.49) | 0.1458 (0.2316; -6.20) | +0.0040 … +0.0060 |
| 9.0 | 0.1962 (0.2626; -5.38) | 0.1989 (0.2766; -6.29) | 0.1664 (0.2462; -6.47) | 0.1391 (0.2155; -6.18) | +0.0025 … +0.0040 |
| 10.0 | 0.1723 (0.2268; -5.45) | 0.1770 (0.2403; -6.33) | 0.1499 (0.2148; -6.49) | 0.1262 (0.1884; -6.22) | -0.0002 … +0.0014 |
| 11.0 | 0.1533 (0.1982; -5.44) | 0.1587 (0.2112; -6.35) | 0.1359 (0.1896; -6.49) | 0.1152 (0.1666; -6.21) | +0.0012 … +0.0012 |
| 12.0 | 0.1374 (0.1751; -5.42) | 0.1433 (0.1875; -6.37) | 0.1238 (0.1690; -6.51) | 0.1055 (0.1488; -6.23) | +0.0009 … +0.0012 |

L²·(ΔV₃,P − ∫∫M) is constant over 8.5 ≤ L ≤ 12 to ±0.06 for every p: the residual is O(1/L²), i.e. (4.2) is the exact first
order of the ratios-conjecture T₃ (§4.2), and the second-order three-point coefficients are ≈ −5.4, −6.35, −6.49, −6.21
(p = 2, 3, 5, 7; computed). At L = 7.5 the pair term is not converged (near-resonant terms, see below).

The second-order structure: ΔV₃,P − ∫∫M ≈ −5.4/L² (p = 2) is the large negative piece that the heuristic misses; the others
(ΔA_off ≈ +2.4/L² · …, Δμ_off, μ's 1/L part, Δμ²) partly cancel it (§5 table).

**Non-additivity** (several bad primes): products of F_p in P1 and the square in −μ² make ΔVar non-additive. The μ² cross
term −2Δμ_pΔμ_q is included exactly; the T₂,off and T₃,P cross terms are computed for the pairs among 2, 3, 5, 7, 11, 13 that
occur (`oned2.json`, `runs/`): ≤ 0.004 (|·|), with a damping spread up to ±0.002 at L ≥ 9 (near-resonant terms
p^a q^b ≈ e^L whose weight e^{−L} is below the conjecture's accuracy) — carried as theory uncertainty. At L = 7.5 the spread
is ±0.015 for (2, 3): the retrospective R120 set (L 7.4–9.5) is less sharply predicted.

## 5. The per-prime change at fixed L (computed; `assemble_out.txt`)

Columns: −½Li₂(1/p) | ΔB_off | ΔA (diag + off) | ΔV₃,Q | ΔV₃,P | −2μΔμ − Δμ² | **total ΔVar** | first order (4.4) = LDM.

| L | p | −½Li₂ | ΔB_off | ΔA | ΔV₃,Q | ΔV₃,P | −2μΔμ−Δμ² | **ΔVar** | first order |
|---|---|---|---|---|---|---|---|---|---|
| 9.0 | 2 | -0.29112 | +0.00322 | +0.01889 | +0.03700 | +0.19622 | -0.14465 | **-0.18043** | -0.14735 |
| 9.0 | 3 | -0.18311 | +0.00389 | -0.00824 | +0.02007 | +0.19889 | -0.12803 | **-0.09652** | -0.06405 |
| 9.0 | 5 | -0.10550 | +0.00406 | -0.02243 | +0.00890 | +0.16636 | -0.09699 | **-0.04560** | -0.01679 |
| 9.0 | 7 | -0.07416 | +0.00392 | -0.02446 | +0.00514 | +0.13911 | -0.07755 | **-0.02800** | -0.00236 |
| 9.0 | 11 | -0.04653 | +0.00355 | -0.02266 | +0.00242 | +0.10394 | -0.05548 | **-0.01477** | +0.00668 |
| 9.0 | 13 | -0.03923 | +0.00338 | -0.02131 | +0.00183 | +0.09221 | -0.04863 | **-0.01175** | +0.00822 |
| 10.2 | 2 | -0.29112 | +0.00250 | +0.02531 | +0.03265 | +0.16810 | -0.12865 | **-0.19122** | -0.16426 |
| 10.2 | 3 | -0.18311 | +0.00303 | -0.00022 | +0.01771 | +0.17302 | -0.11456 | **-0.10412** | -0.07806 |
| 10.2 | 5 | -0.10550 | +0.00317 | -0.01507 | +0.00786 | +0.14690 | -0.08752 | **-0.05016** | -0.02723 |
| 10.2 | 7 | -0.07416 | +0.00306 | -0.01818 | +0.00453 | +0.12387 | -0.07040 | **-0.03128** | -0.01081 |
| 10.2 | 11 | -0.04653 | +0.00277 | -0.01797 | +0.00214 | +0.09379 | -0.05085 | **-0.01665** | +0.00042 |
| 10.2 | 13 | -0.03923 | +0.00264 | -0.01718 | +0.00161 | +0.08364 | -0.04474 | **-0.01326** | +0.00264 |
| 11.4 | 2 | -0.29112 | +0.00199 | +0.02981 | +0.02921 | +0.14671 | -0.11578 | **-0.19918** | -0.17762 |
| 11.4 | 3 | -0.18311 | +0.00243 | +0.00576 | +0.01584 | +0.15227 | -0.10357 | **-0.11038** | -0.08911 |
| 11.4 | 5 | -0.10550 | +0.00254 | -0.00936 | +0.00703 | +0.13088 | -0.07963 | **-0.05405** | -0.03547 |
| 11.4 | 7 | -0.07416 | +0.00246 | -0.01321 | +0.00406 | +0.11123 | -0.06436 | **-0.03399** | -0.01747 |
| 11.4 | 11 | -0.04653 | +0.00222 | -0.01415 | +0.00191 | +0.08496 | -0.04680 | **-0.01840** | -0.00453 |
| 11.4 | 13 | -0.03923 | +0.00212 | -0.01378 | +0.00144 | +0.07611 | -0.04130 | **-0.01463** | -0.00177 |

Reading: the first-order sum (4.4) is +0.127 (p = 2, L = 10.2) against the leading −0.291; the full calculation adds a further
−0.027 (p = 2), −0.026 (3), −0.023 (5), −0.020 (7). The largest second-order pieces are ΔV₃,P − ∫∫M (−5.4/L² for p = 2:
−0.052), ΔA_off (+0.024), the swapped term in Δμ (via −2μ_GUE Δμ_off: +0.014), μ's arithmetic part (−2μ_arith Δμ: −0.007) and −Δμ² (−0.002).
A compact summary (`assemble_out.txt`): ΔVar ≈ −½Li₂(1/p) + a₁/L + a₂/L² + a₃/L³ with a₁ exact (4.4) and a₂, a₃ fitted on
8.5 ≤ L ≤ 12 (misfit ≤ 9·10⁻⁴): p = 2: a₁ 1.2939, a₂ −3.30, a₃ +5.4; p = 3: 1.0715, −3.06, +3.5; p = 5: 0.7984, −2.67, +2.9;
p = 7: 0.6462, −2.38, +2.6; p = 11: 0.4789, −2.04, +2.7; p = 13: 0.4270, −1.84, +1.9. The predictions use the full numbers.

## 6. Round 142, retrospectively (`retro_out.txt`)

All of round 142's L-functions (547; D = 1153 excluded as registered), zero-weighted band averages of ΔVar(S, L(t)) for each
L-function's full bad set, matched to the registered design. **Retrospective**: the measured numbers were known before this
calculation; the calculation has no free parameter and nothing was adjusted after the comparison (the only choices —
damping widths, V, the pair-term averaging — were fixed on convergence grounds at L = 10 before `retro.py` was run).

| p | 2 | 3 | 5 | 7 | 11 | 13 |
|---|---|---|---|---|---|---|
| measured (round 142, bootstrap se) | −0.19422 (0.0016) | −0.10706 (0.0017) | −0.05408 (0.0020) | −0.03256 (0.0024) | −0.01385 (0.0034) | −0.00886 (0.0034) |
| **three-point (this work)** | **−0.19329** | **−0.10559** | **−0.05194** | **−0.03263** | **−0.01720** | **−0.01365** |
| z (measured − three-point) | −0.6 | −0.9 | −1.1 | +0.0 | +1.0 | +1.4 |
| LDM heuristic (registered) | −0.16472 (z −18.4) | −0.07838 (−16.7) | −0.02721 (−13.6) | −0.01131 (−9.0) | +0.00068 | +0.00235 |
| leading order −½Li₂ | −0.29271 (z +61.6) | −0.18243 | −0.10710 | −0.07492 | −0.04776 | −0.03972 |

Amplitude of the beyond-leading-order part, (fit − lead)/(three-point − lead) over p = 2, 3, 5: **0.984** (bootstrap se 0.011
by `var_test.py` on the same data: 0.983 ± 0.011; against the LDM, γ_L = 1.055 ± 0.034). The registered heuristic's
amplitude 0.741 is replaced by 0.98: the three-point route gives the measured drops.

End-to-end check of the test script: `var_test.py` run on round 142's stored per-zero values (`valcheck/`, unpacked from
`r142/vals_*.npz`, `--band 100 600`, `var_test_retro_out.txt`) reproduces round 142's fitted coefficients digit for digit
(−0.19422 … −0.00886) and its matched LDM/lead predictions, so the fresh test reads exactly the registered statistic.

**Mean shift** (normalisation check), matched: +0.04276, +0.03738, +0.02955, +0.02355 (+0.01758, +0.01517) against the
review's +0.0428, +0.0374, +0.0296, +0.0236 — equal; measured +0.04327, +0.03885, +0.03007, +0.02356 (z +0.7, +2.0, +0.6, 0.0).

**L-dependence across the sets** (retrospective; per-set OLS on [1, ⟨ln L⟩, indicators] from the stored values; no se
computed per set, roughly √3 × the pooled se):

| set | ⟨L⟩ | p = 2: fit / 3pt / LDM | p = 3 | p = 5 | p = 7 |
|---|---|---|---|---|---|
| R120 | 9.37 | −0.1853 / −0.1855 / −0.1519 | −0.1003 / −0.0999 / −0.0682 | −0.0477 / −0.0478 / −0.0194 | −0.0260 / −0.0303 / −0.0059 |
| Dn | 10.66 | −0.1961 / −0.1963 / −0.1693 | −0.1080 / −0.1081 / −0.0823 | −0.0554 / −0.0531 / −0.0303 | −0.0339 / −0.0336 / −0.0137 |
| Fn | 10.97 | −0.2014 / −0.1980 / −0.1726 | −0.1135 / −0.1091 / −0.0852 | −0.0592 / −0.0550 / −0.0321 | −0.0382 / −0.0342 / −0.0148 |

The drop deepens with L as predicted (p = 2: −0.185 → −0.196 → −0.201 measured, −0.186 → −0.196 → −0.198 predicted).

## 7. Fresh test: frozen predictions (`fresh.py` → `predictions_varlp.json`)

Sets: **C** = the 181 discriminants of `r134_C.npz` (100 ≤ |D| ≤ 400), **Bp** = the 184 of `r138_Bp.npz` (400 < |D| ≤ 700),
zeros 600 ≤ γ ≤ 1000 (L′ never computed there). Only the keys of the zero files were read. `predictions_varlp.json`
(sha256 eb1f554352da5db63b02008899777df4fba042320ec8d0ac908f344578d3fe29): per L-function D, bad primes, n_zeros_expected,
mean_lnL, cum2_3pt (the frozen prediction), cum2_3pt_numunc (numerical uncertainty), cum2_ldm, cum2_lead, cum2_none,
cum1_3pt (mean shift, secondary). Regression-matched over all 365 (`fresh_matched_out.txt`):

```
L-functions 365, expected zeros 248828 (round 142: 449710) -> se scale 1.344
groups (2,3,5,7,11,13): [121, 92, 61, 44, 30, 26]
mean L range: 9.45 .. 11.39; zero-weighted <L> 10.72

p:                             2         3         5         7        11        13
matched 3pt             -0.19647  -0.10830  -0.05294  -0.03377  -0.01782  -0.01340
matched ldm             -0.16966  -0.08233  -0.03067  -0.01358  -0.00292  +0.00038
matched lead            -0.29328  -0.18433  -0.10626  -0.07629  -0.04528  -0.03762
matched none            +0.00000  +0.00000  +0.00000  +0.00000  +0.00000  +0.00000
numerical unc. of 3pt    0.00071   0.00073   0.00070   0.00070   0.00073   0.00073
expected se              0.00215   0.00231   0.00265   0.00319   0.00462   0.00457
power |3pt - ldm |/se       12.5      11.2       8.4       6.3       3.2       3.0
power |3pt - lead|/se       45.0      32.9      20.1      13.3       5.9       5.3
power |3pt - none|/se       91.3      46.8      20.0      10.6       3.9       2.9
amplitude (0 = lead, 1 = 3pt) over p = 2, 3, 5: expected se 0.0169, separation 59.3 se
amplitude (0 = ldm, 1 = 3pt) over p = 2, 3, 5: expected se 0.0533, separation 18.8 se
```

(Expected se: round 142's bootstrap se × √(449,710/248,828) = 1.344; the group sizes scale the same way, 365/547 L-functions
with ≈ 0.83 of the zeros each.) The theory's own numerical uncertainty (matched, conservative: Σ|H|·unc) is 0.0007 per
coefficient — a third of the measurement se; it is *not* added to the pass bands (registered separately).

## 8. The third cumulant and the mean

**Mean** (check of normalisations): Δμ = Δμ_diag + Δμ_off from the same T₂, band-averaged and matched to round 142's design
(§6), against the review's +0.0428, +0.0374, +0.0296, +0.0236.

**Third cumulant** (no prediction): κ₃ of F over the zeros needs E F³ = ∫ f³R₂ + 3∫∫ f² f R₃ + ∫∫∫ f f f R₄ — the
four-point function with the bad primes deleted, which the project does not have (the ratios recipe gives it, but it is not
coded). The two- and three-point terms alone are not a meaningful part (the R₄ term is of the same order and carries the
leading −κ₃(A_p) = −0.142, −0.054, −0.017). What the three-point function does give is the *first-order structure by
analogy*: in the far-field picture (4.2) the conditional variance of the zeros part depends on the local density, so
κ₃ ⊃ 3 Cov(A_p + ε_p, Var(F|Θ)), the LDM mechanism; its coefficient needs R₄ (the LDM value used ½ where the linear-response
estimate of §4.2 gives 1). Measured −0.037, +0.008, +0.022 for p = 2, 3, 5 against leading −0.142, −0.054, −0.017 and LDM
+0.030, +0.060, +0.057: no reading is predicted here.

## 9. Prior art (searched; abstracts/pages opened as cited)

- Hejhal (1989), and Çiçek, "On the logarithm of the Riemann zeta-function near the nontrivial zeros" (arXiv:2007.10533,
  abstract read): ln|ζ′(ρ)| is approximately Gaussian with variance ½ ln ln T (leading order; RH + pair correlation or a
  spacing hypothesis), with a rate of convergence; no O(1) constant, no arithmetic lower-order term stated.
- Hiary–Odlyzko, "Numerical study of the derivative of the Riemann zeta function at zeros" (arXiv:1105.4312, read via
  ar5iv): near the 10²³rd zero, mean 3.4907, sd 1.0977 against the leading-order 2.0 and 1.4 and the CUE values 3.47, 1.12;
  no lower-order theory for these two statistics; differences "probably not significant".
- Lugar–Milinovich–Quesada-Herrera (arXiv:2211.14918, abstract read): lower-order terms in the second moment of
  ln|ζ(1/2+it)| on the line, related to the pair correlation — the two-point analogue on the critical line, not at zeros.
- Conrey–Snaith, triple correlation with lower-order terms (arXiv:math/0610495; 0803.2795): the input, not applied to ln|ζ′|.
- Hughes–Keating–O'Connell (2000), Bui–Gonek–Milinovich (hybrid model at zeros): leading-order arithmetic factors of the
  moments; the leading variance change ½Li₂(1/p) follows (KNOWN IN SUBSTANCE).
- Searches for log|L′(ρ, χ)| with conductor dependence found only CLTs for ln|L(ρ, χ)| at zeros of ζ (arXiv:2109.09097,
  2504.08322) and low-lying-zero statistics.

Verdict: the leading-order bad-prime factor of the variance (−½Li₂(1/p)) — KNOWN IN SUBSTANCE (HKO/hybrid model). The exact
first-order formula (4.4) and the finite-L three-point evaluation — NOT FOUND (nothing on lower-order terms of
Var ln|ζ′(ρ)| or on its computation through the three-point function was found).

## 10. Doubts

1. **Conjecture.** Everything finite-L rests on the ratios-conjecture T₂, T₃ with the bad primes deleted. Rounds 134–138
   confirmed them in windows of ≤ 8 mean spacings (raw ε ≲ 5); here they are integrated against ln|u| ln|v| over all
   separations (the relevant range is ε ≲ 30–50 raw). The diagonal far field is Berry/Landau–Gonek (safe); the off-diagonal
   far field enters ΔA_off, ΔB_off (+0.024, +0.003 for p = 2) with no independent test.
2. **Near-resonant terms.** Off-diagonal Euler-product terms with frequency L − ln n ≈ 0 (n = p^a q^b ≈ e^L) carry weight
   ~e^{−L}, below the conjecture's accuracy, but are amplified by the 1/τ² of the log kernel. They make the pair terms and
   very large single primes (ln p close to L) converge poorly under damping: carried as uncertainty (±0.002 at L ≥ 9; ±0.015
   for (2,3) at L = 7.5, which affects only the retrospective R120 set; single p ≥ 151 at L ≤ 8 likewise).
3. **Fixed L.** The correlation functions are taken at the height of the zero; zeros 30–50 raw units away sit at heights
   differing by 3–8 % (fresh sets) or up to 50 % (γ ≈ 100). The resulting error is second order in 1/L times ΔL/L: estimated
   ≲ 0.002 at γ ≥ 600, larger for R120 at low heights (not computed).
4. **Additivity.** Base S = {} plus pair corrections among 2…13 (the μ² cross term exact); triple and large-prime cross
   terms neglected (each pair ≤ 0.004).
5. **First-order law (4.2)** is derived from the warped-sine (BK) picture and confirmed only through the 1/L² scaling of the
   residual (six L values); the pointwise far field is masked by the oscillating off-diagonal part (`farfield_out.txt`).
6. **Retrospective.** The agreement in §6 is with data seen before; it is parameter-free, but the convergence choices were
   made knowing round 142's numbers existed (not their comparison). Hence the fresh registration.
7. **Bootstrap se** of var_test with 500 resamples (validation) differ from round 142's 2000-resample se by ≤ 12 %; the
   registered run uses 2000.
8. **Third cumulant**: no prediction (needs R₄).
