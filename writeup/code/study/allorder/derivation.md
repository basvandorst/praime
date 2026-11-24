# The long-range Hardy–Littlewood term of the within part, per direction and without truncation

> **In short.** (1) The within marks are doubly centred: E[z | class of p_n] = E[z | class of p_(n+1)] = 0 (*proved*). So every four-pair term of the Hardy–Littlewood expansion at the level of residue classes mod q vanishes in all (C − 1)² directions, including the cross pairings that carry Jacobi sums. The ln H term comes only from finer residues, through the insertion-response power (*heuristic*: local-density form of HL). (2) The response power that one more prime p adds is isotropic: p·inc_d(p) is the same in every direction to 1–2 %, and the same for q = 3, 5, 7 when p ≥ 29 (*computed*). So the calculation predicts **no mixed-direction excess**, as round 139 found. (3) Without truncation the first-order (pair-correlation) sum does **not** converge: inc(p) ≈ κ/p for every p up to 499, and the power of p switches on in the measured range of H. Over [10⁴, 10⁶] it is a closed, convergent sum, s = LR₁₇ − 0.625 κ/λ, and it overshoots the fresh slopes by 3–5 errors (§5.2). (4) The all-order local model (every local factor up to P = 10⁵, random profinite start) is converged in P to ±0.0005 for H ≤ 10⁶. *Retrospectively*, with nothing fitted, it gives the fresh round-139 slopes in ten (q, k) cells with χ² = 14.5/10: slope/LR = 1.6, 1.7, 1.4, 1.3 for q = 7, 9, 11, 13 at k = 1 (measured 2.0, 1.6, 1.6, 1.4) and 1.4–2.4 at k = 2 (measured 0.9 ± 0.9, 1.3 ± 0.9, 1.8 ± 0.3, 1.5 ± 0.3), with W(10⁴) to 0.0004. So the common k = 1 deficit is a Hardy–Littlewood effect of all primes at all orders, beyond the pair-correlation tail; k = 2 is predicted to show it too. A registration draft for a third height is in round_draft.txt.

Status labels: **proved** (algebra), **HL** (follows from the Hardy–Littlewood conjecture), **heuristic**, **computed** (a script here), **retrospective** (compared with numbers known before the computation). All scripts and outputs are in this folder. Real primes were listed only in A = [1.5·10¹¹, 1.5·10¹¹ + 10⁸) and B = [5·10⁹, 5·10⁹ + 10⁸) (`sieve.py`, checked by `check_primes.py` → `check_primes_output.txt`). The model runs sieve only rough numbers at random 120-digit starting points; they list no primes.

## 1. Setting

Blocks of length H near x, λ = ln x, N = H/λ primes per block. c_n = class of p_n mod q (C = φ(q) classes), k = 1. Within directions as in review44/dirs.py: for non-principal characters χ_i, χ_j the function χ_i(a)·conj χ_j(b), weighted by √(π_a P_ab), projected by Π_W (row and column parts removed), Löwdin-orthonormalised; `common.py` re-types that construction. The per-step mark of direction d is z_d(a, b) = B_(d,ab)/√(π_a P_ab), and W_d(H) = E[z-block²] = 1 + 2 Σ_j (1 − j/N) τ_j(d) (repo derivation §2, proved up to O(1/N)). τ₀(d) = 1 exactly on both stretches (`lagcov.py`).

## 2. Exact algebra: the class-level four-pair terms vanish (proved)

**2.1 Double centring.** B_d lies in the range of R = ⊕_a (I − √P_a √P_aᵀ), so Σ_b √P_ab B_(d,ab) = 0 for every a; and Π_W annihilates RV, so Σ_a √(π_a P_ab) B_(d,ab) = 0 for every b. Hence

  Σ_b P_ab z_d(a, b) = 0 for every a,   Σ_a π_a P_ab z_d(a, b) = 0 for every b.   (2.1)

Check (`algebra.py` → `algebra_output.txt`): both sides are below 3·10⁻¹⁶ for q = 5, 7 on both stretches.

**2.2 Consequence for the four-point expansion.** Take two consecutive pairs far apart, (p, p′) with classes (a, a′) and (r, r′) with classes (b, b′). In the Hardy–Littlewood four-point density with empty interiors, the first-order connected part is a sum of pair terms. At the level of classes mod q a pair term couples one class of the first event with one class of the second: it is a weight w(x, y) with x ∈ {a, a′}, y ∈ {b, b′}, multiplied by the product of the two event laws. Summing the marks against it, the mark of event 1 is averaged over its other class with the transition law, and (2.1) gives zero. This holds for all four pairings (p–r, p–r′, p′–r, p′–r′), for any weight, and so for any gap-shifted version that still depends only on the two coupled classes. Check: the largest class-level pairing term over random weights is 2.5·10⁻¹⁷ (`algebra_output.txt`).

So **no direction gets a ln H term from same-class pairings or from the cross pairings at the level of classes**. The Jacobi sums Σ_a χ_i(a) conj χ_j(a + d) (|J| = √5 for i ≠ j mod 5, 1 for i = j; `algebra_output.txt`) appear only if the coupling weight also depends on the gap g or on the interior, that is, on residues finer than mod q. That dependence is the insertion response of §3.

**2.3 Parities.** Under (a, b) → (−a, −b) a direction picks up χ_i(−1)χ_j(−1). The odd directions are exactly the four mixed ones for q = 5; none for q = 3; 12 of 25 for q = 7; 40 of 81 for q = 11; 60 of 121 for q = 13 (`algebra_output.txt`). The law is invariant only under the reflection (a, b) → (−b, −a), which maps the group (1,2) to (2,1); the per-direction expectations of the two mixed groups are therefore equal, and nothing singles out the odd sector.

## 3. The long-range term is the coherent response power (heuristic: local-density form of HL)

Condition on the residues of all sites modulo a primorial Q; in the local-density model the sites are then independent, and the block sum of the marks has conditional mean Σ_x f(x), f = local mean of the event at x. Its variance over block positions is Σ_θ |f̂(θ)|² |D_H(θ)|². Only frequencies within 1/H of 0 grow with H. Writing f to first order in the site densities, the frequencies near 0 combine the contributions of all events that see one site u coherently, which gives the insertion response r_u (review44, repo §4.4):

  dW_d/d ln H = −(1/λ) × [power of r_d at the denominators resolved at H].   (3.1)

With the residues mod Q = 510510 exact and the pair factors of the primes > 17 this is review44's tail, slope −V_d/λ at large H. Within the y-truncated model, (3.1) per direction (`tailM.py` → `tailM_B.json`): every direction's finite-H slope is 0.89–0.93 of −V_d/λ; no direction type turns on differently.

## 4. Per direction: the large-prime power is isotropic (computed), hence no mixed excess

Moving one more prime p from the pair factors into the exact residues adds response power inc_d(p) = V_d(210p) − V_d(210). On both stretches (`respdir.py` → `respdir.json`, `respsum.py`):

| | λ | p·inc(p), p = 19–61, per direction, by group |
|---|---|---|
| q = 5, A | 25.7 | ratio 0.269, mixed 0.270, product 0.265, (2,2) 0.263 |
| q = 5, B | 22.35 | ratio 0.330, mixed 0.332, product 0.331, (2,2) 0.329 |
| q = 7, A | 25.7 | all 13 groups 0.254–0.264 |
| q = 7, B | 22.35 | all 13 groups 0.312–0.321 |
| q = 3, A / B | | 0.244 / 0.317 |

Checks (`inccheck.py`): a cross-stretch estimate (no noise correction) gives p·inc = 0.303, 0.298, 0.297, 0.288 for p = 19, 47, 97, 199; replacing the residue mod p by a random one gives |p·inc| ≤ 0.0014.

*Why isotropic (heuristic).* A large prime removes one site at offset s, and the class of that site, c_u + s mod q, is uniform over the classes as u varies. The induced change of the response is a difference of rows of F, and its power averaged over a uniform class is proportional to ‖B_d‖² = 1 in every direction.

Per group mod 5 at λ = 22.69, slope per unit ln H (`pergroup.py` → `pergroup_output.txt`; measured values retrospective):

| group | all-order model CG (P = 10³ and 10⁴, 102 windows) | registered 24 windows | old height, 9 windows (λ 27.0) | y = 17 LR −V_d/λ |
|---|---|---|---|---|
| (1,1) ratio | −0.0093 ± 0.0019 | −0.0070 ± 0.0029 | −0.0143 ± 0.0066 | −0.0066 |
| (1,2) mixed | −0.0113 ± 0.0021 | −0.0193 ± 0.0035 | −0.0117 ± 0.0058 | −0.0062 |
| (1,3) product | −0.0077 ± 0.0019 | −0.0021 ± 0.0030 | +0.0007 ± 0.0050 | −0.0038 |
| (2,1) mixed | −0.0078 ± 0.0019 | −0.0091 ± 0.0040 | −0.0148 ± 0.0045 | −0.0062 |
| (2,2) ratio | −0.0088 ± 0.0031 | −0.0070 ± 0.0067 | +0.0108 ± 0.0106 | −0.0048 |

The model has no group pattern beyond its errors: (1,2) and (2,1) are reflection partners with equal expectations, and their model difference is 1.2σ. Against it, the registered (1,2) group is 2.0σ low and the product group 1.5σ high. The model was not run at λ = 27.0.

So the directions differ only through the small primes (V_d at y = 17: mixed 0.140, ratio (1,1) 0.150, product 0.087, (2,2) 0.108 per direction, review44), and **the Hardy–Littlewood calculation predicts no excess in the mixed directions mod 5**. The all-order model of §5.3 agrees: odd minus even slope mod 5 is +0.0005 ± 0.0021 (CG(1000), 62 windows) and −0.0002 ± 0.0039 (CG(10⁴)). Round 139 found −0.0003 ± 0.0030 on fresh windows.

## 5. The slope without truncation

### 5.1 First order: the sum over primes does not converge (computed)

The increments do not fall faster than 1/p: p·inc(p) is flat from p = 19 to 199 (cross-stretch 0.30 at λ ≈ 24, about 7.2/λ per direction). In the exact finite-H formula with the residues mod 210p exact, the power of p switches on inside H ∈ [10⁴, 10⁶] (turn-on fraction f_p of the slope over that range, `tailbig.py`, `tailx.py` → `turnon_output.txt`, `tailx1.log`): see the table in §5.2. Therefore

  s_∞(H₁, H₂) = s_y(H₁, H₂) − (1/λ) Σ_(p>y) inc(p) f_p(H₁, H₂),   inc(p) = κ/p,   (5.1)

is a convergent sum for a fixed range (f_p → 0 for p ≫ H₂/c), but its limit as the range grows diverges like (κ/λ) ln ln H. **There is no convergent constant slope at first order; the comparator's growth with y (6–33 % per step) is these increments, and it would continue for every prime up to about H₂/c.** The comparator's own steps give κ directly: p·[LR(y = p) − LR(y = previous)] at λ = 24.16 is 0.011 (q = 7), 0.007–0.010 (q = 11), 0.0065 (q = 13) per unit ln H (`fresh139_cells.txt`).

### 5.2 Turn-on and the closed form

Turn-on of one prime's power, measured without any noise correction (responses on stretch A times responses on stretch B, `tailx.py` → `tailx1.log`, q = 5, all directions, λ = 22.35):

| p | 19 | 47 | 97 | 199 | 499 |
|---|---|---|---|---|---|
| p·inc(p) per direction | 0.349 | 0.343 | 0.342 | 0.331 | 0.336 |
| fraction of −inc(p)/λ in the slope over [10⁴, 10⁶] | 0.95 | 0.91 | 0.74 | 0.74 | 0.43 |

A one-constant model f_p = clip(ln(H₂/(c p))/ln(H₂/H₁), 0, 1) fits with c = 266 (heuristic: the power of p sits at denominators p·e₀ with e₀ | 210, and a denominator is resolved at H ≈ c·p). Then, with Mertens' density of primes,

  Σ_(p>17) f_p/p = 0.625 (c/2: 0.736; 2c: 0.502)  for [H₁, H₂] = [10⁴, 10⁶]   (`firstorder.py` → `firstorder_output.txt`),

and the closed form is s_∞ = LR₁₇ − (κ/λ)·0.625, with κ/λ per direction read off the comparator's own steps (p·ΔLR/f_p at p = 13, 17). At λ = 24.16 this gives slope/LR = 2.31 (q = 7), 2.33 (9), 1.93 (11), 1.70 (13) for k = 1 and 2.73, 2.61, 2.22, 1.86 for k = 2, against the fresh 2.03, 1.59, 1.63, 1.37 and 0.92, 1.27, 1.78, 1.51 (retrospective). **The first-order sum without truncation overshoots**: for q = 9, 11, 13 at k = 1 by 4, 3 and 5 measurement errors, and the measured excess over the y = 17 comparator is 0.4–0.8 of the first-order excess. The ordering in q is right (the excess shrinks as q grows because V_17 grows with q while κ does not).

So the first-order expansion is not the answer: the increments are real (§4), they do switch on in the measured range, and yet the slope of the primes is smaller than their sum. The missing piece is of higher order in the local densities. §5.3 computes it.

### 5.3 All orders: the profinite Cramér–Granville model (computed)

CG(P): a window of 10⁸ starts at a random 120-digit integer (so the residues of all sites are a random profinite position: no sieve-limit effects); a site with no prime factor ≤ P is a pseudo-prime with probability ∏_(p≤P)(1 − 1/p)⁻¹/λ_ref; same statistic as review44/dirs.py (`cgsim3.py`, `cgsim4.py`, summaries `cgPsum.py`). It contains every order of the local factors up to P and nothing beyond.

**Convergence in P at λ = 22.69** (`cgPsum_output.txt`; slope per unit ln H, all directions; error from the window scatter):

| P | windows | q = 3 | q = 5 | q = 7 | W(10⁴), W(10⁵), W(10⁶) for q = 5 |
|---|---|---|---|---|---|
| 10³ | 62 | −0.0023 ± 0.0039 | −0.0084 ± 0.0012 | −0.0084 ± 0.0008 | 0.9760, 0.9618, 0.9373 |
| 10⁴ | 40 | −0.0038 ± 0.0042 | −0.0099 ± 0.0013 | −0.0089 ± 0.0008 | 0.9756, 0.9591, 0.9300 |
| real primes, 24 registered windows (review44 dirs.json) | | +0.0088 | −0.0091 | −0.0094 | 0.9759, 0.9547, 0.9339 |
| comparator y = 17 (round 139 map, new height) | | −0.0028 | −0.0051 | −0.0064 | |

(An earlier set of runs started at x₀ ≈ 3·10⁹ instead of a random profinite position (`cgsim.py`, `cgcompare_output.txt`) is not used: for P ≥ 10³ at that height more than half of the rough numbers are real primes, so those runs are partly real data. They gave the same picture: P = 17, 100, 10³ give slopes −0.0012, −0.0048, −0.0086 for q = 5.)

**At λ = 24.16, q = 5–13, k = 1 and 2** (`cg4sum.py` → `cg4sum_output_24.txt`), P = 10³ / 10⁴ / 10⁵ with 70 / 71 / 105 windows. Between P = 10⁴ and 10⁵ the slopes move by at most 0.0008, within the Monte Carlo errors (0.0003–0.0010). From 10³ to 10⁴ they move by up to 0.0013. So for H ≤ 10⁶ the model is converged in P to about ±0.0005; that is the law's own uncertainty below.

**What the model contains that the comparator lacks.** Every order of the local factors of every prime up to P, at once. The first-order sum of §5.1 keeps only one prime's modulation at a time, on top of pair factors for all the others. The two differ by negative higher-order terms that are about as large as the first-order excess beyond y = 17, so no closed formula is given here (doubt 3).

**Per direction.** Odd minus even slope mod 5: +0.0005 ± 0.0021 (P = 10³), −0.0033 ± 0.0028 (P = 10⁴). Mod 7: +0.0017 ± 0.0016, −0.0028 ± 0.0019. No sector is special, as §4 says.

## 6. Comparison with round 139 (retrospective)

*Retrospective*: the fresh numbers (round 139, 24 windows at λ = 24.16, `writeup/data/study/r139/rep_fresh.json`, read only) were known when the model was run. Nothing was fitted to them: P is the largest value the model allows comfortably (ρ = 0.85) and the λ is the windows' mean.

| cell (q, k) | fresh slope ± err | comparator y = 17 | first-order, no truncation (§5.2) | all-order model CG(10⁵) ± MC | z (fresh − CG) | slope/LR fresh | slope/LR CG |
|---|---|---|---|---|---|---|---|
| 5, 1 | −0.0053 ± 0.0019 | −0.0043 | −0.0115 | −0.0088 ± 0.0009 | +1.7 | 1.24 | 2.06 |
| 7, 1 | −0.0108 ± 0.0012 | −0.0053 | −0.0123 | −0.0086 ± 0.0006 | −1.6 | 2.03 | 1.61 |
| 9, 1 | −0.0080 ± 0.0009 | −0.0050 | −0.0117 | −0.0084 ± 0.0006 | +0.3 | 1.59 | 1.67 |
| 11, 1 | −0.0094 ± 0.0006 | −0.0058 | −0.0112 | −0.0079 ± 0.0003 | −2.3 | 1.63 | 1.37 |
| 13, 1 | −0.0079 ± 0.0004 | −0.0058 | −0.0098 | −0.0074 ± 0.0002 | −1.0 | 1.37 | 1.29 |
| 5, 2 | −0.0044 ± 0.0020 | −0.0013 | −0.0038 | −0.0038 ± 0.0010 | −0.2 | 3.53 | 3.08 |
| 7, 2 | −0.0015 ± 0.0015 | −0.0017 | −0.0045 | −0.0040 ± 0.0006 | +1.6 | 0.91 | 2.43 |
| 9, 2 | −0.0021 ± 0.0015 | −0.0017 | −0.0044 | −0.0024 ± 0.0005 | +0.1 | 1.27 | 1.40 |
| 11, 2 | −0.0037 ± 0.0007 | −0.0021 | −0.0046 | −0.0038 ± 0.0003 | +0.1 | 1.79 | 1.84 |
| 13, 2 | −0.0034 ± 0.0006 | −0.0022 | −0.0041 | −0.0035 ± 0.0003 | +0.2 | 1.51 | 1.57 |

z uses the fresh error and the model's MC error. Over the ten cells χ² = 14.5 for the model (p ≈ 0.15). For the comparator alone the k = 1 cells are at z = −0.5, −4.6, −3.3, −6.0, −5.4. The levels agree too: W(10⁴) of the model against the fresh windows is 0.9791/0.9777, 0.9823/0.9823, 0.9820/0.9824, 0.9839/0.9837, 0.9853/0.9856 for k = 1 (q = 5, 7, 9, 11, 13).

**Answer to the question.**
- The k = 1 deficit **is** a Hardy–Littlewood effect: the local model with all local factors gives slope/LR = 1.6, 1.7, 1.4, 1.3 for q = 7, 9, 11, 13, against the measured 2.0, 1.6, 1.6, 1.4.
- It is **not** given by the first-order (pair-correlation) sum without truncation. That sum does not converge as the range of H grows. Over [10⁴, 10⁶] it is a convergent sum, and it overshoots (1.9–2.3).
- For k = 2 the model does **not** predict slope/LR ≈ 1. It predicts 1.4–2.4, about as large as for k = 1. The fresh k = 2 cells agree with it (q = 9, 11, 13 within 0.2 errors; q = 7 at +1.6). "k = 2 shows nothing" in round 139 is low power, because the k = 2 comparator is 3 times smaller.
- No direction type is singled out (§4).

## 7. The other candidates of the brief

- **(a) Three-point connected terms.** In the local-density form they are the second-order terms of f in the site densities (two sites of one event, one site of the other). They are inside the all-order model of §5.3 and inside the increments inc(p) of §5.1 (a prime that removes the partner site of an event). They are not a separate ln H term beyond those.
- **(b) Estimated P̂ and the factor NB/(NB − 1).** With block correlation ρ̄ inside a window, the statistic estimates σ²(1 − ρ̄), and (NB − 1)ρ̄ = W(window)/W(H) − 1 ≈ −0.1 at H = 10⁶. That is a bias of +0.001 at 10⁶ and about 0 at 10⁴: at most +0.0002 per unit ln H, of the opposite sign to the excess.
- **(c) Products of rigid class contrasts.** They enter through the ratio estimator (actual N_a in the denominators) at relative order 1/N: −0.001 in level at 10⁴ (round 139's chain null) and a slope consistent with 0 (chain null, 55 cells, rms z 1.00). Upper bound 0.0003 per unit ln H. They are not direction-specific beyond the column projection, which the double centring (2.1) makes exact.

## 8. Prior art

Mathematical queries only.
- **Lemke Oliver–Soundararajan, arXiv 1603.03720** (read through the ar5iv HTML; the PDF could not be text-extracted here). Main Conjecture, Prop. 2.1 and (2.23): c₂(q;(a, b)) = (log 2π)/2 + q·S₀ᶜ(q, b − a) + q·B_q(b − a) − (q₀/φ(q₀))Σ_(d|q₀)(μ(d)/φ(d))Σ_χ C_(q₀,χ)(χ̄(b) − χ̄(a)), with B_q(v) = 1/2 − v/q and S₀ᶜ(q, v) = −(φ(q)/(2q))Λ(q/d)/φ(q/d) − B_q(v) + (1/φ(q/d))Σ_(χ≠χ₀) χ̄(v/d)L(0, χ)L(1, χ)A_(q,χ). The sawtooth cancels; the dependence on the difference b − a goes through odd characters only (L(0, χ) = 0 for even χ). This is the class-level, order-1/λ part of the transition law; it enters the within part only through the class-level response power (part of V), and by (2.1) it gives no class-level ln H term. LOS treat no block variances.
- **David–Devin–Nam–Schlitt, arXiv 2105.05048** (abstract): the LOS heuristic for consecutive sums of two squares; averages of the Hardy–Littlewood constant; no variances or long-range correlations.
- **Montgomery–Soundararajan, math/0409258** (abstract and search text): the variance of prime counts in short intervals from sums of singular series. This is the d ↦ μ²(d)/φ(d) log density used in §3.
- **Freiberg, arXiv 2609.33692** (abstract): under a uniform HL hypothesis the arithmetic correction to the distribution of prime counts in intervals of the mean-gap scale differs from the binomial correction of Cramér's model; finite-sieve argument; a refined random model. This is a relative of §5.3's observation that the all-order local model and its first-order expansion differ.
- **Granville–Lumley** and later sieve-based random models (search snippets only): sieve by small primes, keep survivors independently. §5.3's CG(P) is of this type, started at a random profinite position.

Verdict. (2.1)–(2.2) (the class-level four-pair terms vanish by double centring): elementary, **not found** stated. The isotropy of the large-prime response power and the non-convergence of the first-order sum: **not found**. Using a sieve random model for the block-length dependence of transition-table statistics: **not found**; random models of this kind are **known** (Granville–Lumley; refined models in Freiberg 2026).

## 9. Doubts

1. **Is the local model the Hardy–Littlewood conjecture?** CG(P) reproduces the singular series of every finite pattern for the primes ≤ P exactly, and none beyond. For H ≤ 10⁶ the primes that matter are those whose power switches on in the range (p ≲ H₂/c ≈ 4·10³, §5.2), all below P = 10⁵. Primes above P enter only through denominators that are not resolved. The model fails as P approaches exp(λ e^(−γ)) (7.8·10⁵ at λ = 24.16), where its density reaches 1. This is a heuristic identification, not a theorem; Freiberg (2026) shows that such models and HL can differ in second-order terms.
2. **Retrospective.** The 22.69 and 24.16 comparisons were made with the real numbers known. No parameter was fitted. But presenting the model as the answer was decided after its 22.69 agreement, and the choice P = 10⁵ after seeing that 10³ was too small at 24.16. The registration below is the honest test.
3. **No closed form at all orders.** The first-order formula (5.1) is closed and convergent over a fixed range, and it is rejected (§6). The all-order slope is a model computation, and why the higher-order terms nearly halve the first-order excess was not derived. A plausible cause is saturation: with many primes exact, the remaining rough sites are fewer and each is more likely prime, so removing one more site matters less. The increments were measured relative to the mod-210 classes only (ρ = 0.2), where this cannot show.
4. **The turn-on constant** c = 266 rests on five primes. Their fractions (0.95, 0.91, 0.74, 0.74, 0.43) are not monotone at 97 and 199. c/2 and 2c change Σ f_p/p from 0.625 to 0.74 and 0.50, and the first-order excess by ±20 %. The rejection of the first-order sum holds for all three.
5. **κ from the comparator's steps.** Only two steps (y = 11 → 13 → 17) were available, and q = 13 has only one. The cross-stretch increments (p·inc = 0.34 at λ = 22.35, flat to p = 499) agree with them within 10 %.
6. **Isotropy** is measured for p ≤ 61 in every direction (both stretches) and in the totals to p = 499. The argument for it (uniform class of the excluded site) is heuristic.
7. **q = 3.** At λ = 22.69 the model gives −0.0038 ± 0.0042 and the registered windows +0.0088 ± 0.0064 (1.6σ). Fresh: −0.0041 ± 0.0050. The single within direction of q = 3 is the noisiest cell. It was not run at 24.16.
8. **The k = 2 comparator** (round 139's response with second neighbours) was not re-derived here. respk2.py confirms the scale on stretch B: p·inc ≈ 0.13 for k = 2 against 0.33 for k = 1, with V(210) 0.011 against 0.051. So relative to its own comparator k = 2 has the larger untruncated excess.
9. **My lag-covariance bins** on the two stretches (`lagbins.py`) are too noisy (±0.005–0.01 per sector) to test anything at the level of the excess. They are reported only as a check that τ₀ = 1 and that the lag-1 term sits in the product and (2,2) directions (−0.020, −0.015 for q = 5, as the leading-order theory says).

## 10. Files (all in this folder)

- **Sieve and checks.** `sieve.py` (stretches A, B only; the primes_*.npy were deleted afterwards and can be rebuilt), `check_primes.py` → `check_primes_output.txt`, `common.py` (basis as review44).
- **Exact algebra.** `algebra.py` → `algebra_output.txt`.
- **Lag covariances per direction.** `lagcov.py` → `lagcov.json`; `lagbins.py` → `lagbins.json`.
- **Response power per direction.** `respdir.py` → `respdir.json`; `respsum.py`; checks `inccheck.py`; k = 2 `respk2.py` → `respk2_B5.log`.
- **Turn-on.** `tailM.py` → `tailM_B.json`, `tailbig.py` → `tailbig_B.log`, `turnon.py` → `turnon_output.txt`, `tailx.py` → `tailx1.log`, `tailx_19-47-97-199-499.json`.
- **First-order closed form.** `firstorder.py` → `firstorder_output.txt`, `firstorder.json`.
- **All-order model.** `cgsim3.py`, `cgsim4.py` → `cgP_*.json`, `cg4_*.json`; summaries `cgPsum.py` → `cgPsum_output.txt`, `cg4sum.py` → `cg4sum_output_24.txt`, `pergroup.py` → `pergroup_output.txt`, `oddeven.py`. Superseded runs from x₀ ≈ 3·10⁹: `cgsim.py`, `cgsim2.py`, `cgcompare.py` → `cgcompare_output.txt`.
- **Round 139 numbers used** (read only): `fresh139_cells.txt`.
- **Registration draft.** `round_draft.txt`, numbers `reg_numbers.py` → `reg_numbers_output.txt`, `reg_numbers.json`.
