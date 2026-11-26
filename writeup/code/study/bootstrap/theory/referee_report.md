# Referee report on theory/proofs.txt (sections 1-5) and on section 5b of bootstrap.md (26 November 2026)

Brief: assume each proof wrong. Verdict per result:

| result | verdict |
|---|---|
| Lemma 1 (the comb determines the system) | no failing step found |
| Lemma 2 (D_primes tempered iff RH) | no failing step found |
| Theorem A (antichain) | no failing step found |
| Theorem C (thin changes are rigid) | gaps but fixable |
| Proposition 5 | gaps but fixable |
| 5b (price of an edit) | main statement holds; one consequence false as worded, one with a fixable gap |

1. **Lemma 2.** The cutoff argument is valid; signs are right (L = 1/(w - 1/2) + zeta'/zeta(w + 1/2) + E); a zero with real part above 1/2 is a pole with residue m >= 1 that 1/(s - 1) cannot cancel. Notation: "chi*T" is a product.
2. **Theorem A.** The sign nu = C_Q - C_P = D_P - D_Q is correct; complex test functions are fine (D_P is real and even); sigma >= 0 with a vanishing integral sits on t0 + Z, Z discrete; a shift t1 outside Z - Z gives sigma = 0.
3. **Theorem C.** (i) holds, sign slip (D_primes = D_P - nu). (ii) holds: |sigma| <= mu_P + mu_0 is locally finite and the indicator of [gamma* - 1, gamma* + 1] dominates phi_T. (iii) holds. Case (a) holds. Case (b): the per-term bound is false uniformly in k for eta < 0 (p = 2, eta = -0.5: ratio 3.8 at k = 20, 2900 at k = 50); fix with min(p, p e^eta)^(-k/2) and a constant depending on sup|eta|; require p e^eta > 1. Additive moves o(sqrt p) qualify (Mertens). "Sharpness" is only the contrapositive of (ii).
4. **Proposition 5.** (i), (ii), (iv) hold. (iii): Wiener-Ikehara in Stieltjes form applies. (v): a_C undefined in the text; with the Gamma_C term the conclusion holds.
5. **Section 5b.** Steps (3)-(5) hold; the identity of step (3) was checked for complex f to 6e-8. RH is used in step (5) as well as in step (2). Step (6): the extremes are approached, not attained. "Moving a prime by any amount costs both prices" is false when log x / log p = m/n with m odd (2 -> 2^(11/10): 2 inf h = -3.85771 against the sum -3.86009; 2 -> 8: -2.811 against -2.849); true for an irrational ratio or m even. "Every finite edit needs q > 1" was not proved for dependent logarithms; fix: h_P is almost periodic with mean 0 and h_P = 0 forces equal combs (Lemma 1), or cite Theorem C (a).

Not checked: the normalisation of Weil's formula (the 27-digit check was accepted); the theorem numbers cited in Schwartz and Korevaar; the table of edit_cost_out.txt; sections 5c, 6, 7 of the note; the Herglotz remark.
