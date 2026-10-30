# Facts about primes

A long list of facts about prime numbers, each with its formula. The emphasis is on what matters for puzzles about runs of consecutive primes: residues, gaps, run lengths, patterns, formulas, and the open questions ("missing parts").

**How to read the tags.** Every fact ends with its status:

- *theorem* (who, year): proved; we did not reprove it.
- *checked*: the statement is a theorem or a claim from the literature, and our script confirmed it on the stated range.
- *computed*: a number we found ourselves with the C programs listed at the end.
- *literature*: a number we did not recompute. The year says when it was true.
- *conjecture*, *heuristic*, *conditional*: not proved.
- *site*: from a primepuzzles.net page. *ours*: from a write-up in [results/](results/) or [writeup/](writeup/).

**Notation.** p, q are primes, p_n is the n-th prime (p₁ = 2), g_n = p_(n+1) − p_n is the gap after p_n. π(x) counts primes ≤ x and π(x; m, a) those ≡ a (mod m). li(x) = ∫₀ˣ dt/ln t (principal value), (a/p) is the Legendre symbol, k# is the product of the primes ≤ k, φ is Euler's function, λ is Carmichael's function, v_p is the p-adic valuation, γ = 0.5772156649… is Euler's constant.

---

## 1. Residues: what a prime is modulo m

### 1A. "A prime mod m is always r"

1. Every prime p > 3 is 6k ± 1: `p ≡ 1 or 5 (mod 6)`, and so `p ≡ 1, 5, 7, 11 (mod 12)`. The other residues are multiples of 2 or 3. *(theorem)*
2. Primes p > 5 lie in the 8 classes coprime to 30: `p mod 30 ∈ {1, 7, 11, 13, 17, 19, 23, 29}`. Primes p > q lie in the φ(q#) classes coprime to q#, a fraction `∏_(p≤q) (1 − 1/p) ~ e^(−γ)/ln q` of all residues (Mertens). *(theorem)*
3. Primes p > 5 end in 1, 3, 7 or 9, and for p > 3 the digit sum is not divisible by 3. *(theorem)*
4. `p² ≡ 1 (mod 24)` for every p > 3, and 24 is the largest such modulus: `gcd(p² − 1 : p > 3) = 24`. *(theorem; checked p < 10⁵)*
5. Higher powers: `p⁴ ≡ 1 (mod 240)` for p > 5, `p⁶ ≡ 1 (mod 504)` for p ∉ {2, 3, 7}, `p⁸ ≡ 1 (mod 480)` for p > 5, each the largest possible modulus. Rule: the largest m with `p^k ≡ 1 (mod m)` for all p ∤ m is the largest m with λ(m) | k. *(theorem; the four gcds checked p < 10⁵)*
6. `p² + 2 ≡ 0 (mod 3)` for p > 3, so p² + 2 is prime only for p = 3. Likewise p² + k is divisible by 3 whenever k ≡ 2 (mod 3). *(theorem)*
7. `2^p + 1 ≡ 0 (mod 3)` for odd p, which is why Wagstaff numbers are written (2^p + 1)/3. *(theorem)*
8. Gap residues: for p > 3 the gap to the next prime is ≡ 0 (mod 6), or ≡ 2 (mod 6) exactly when p ≡ 5 (mod 6), or ≡ 4 (mod 6) exactly when p ≡ 1 (mod 6). *(theorem; checked below 2·10⁸)*
9. Twin primes (p, p + 2) with p > 3 have `p ≡ 5 (mod 6)`, so `p + 1 ≡ 0 (mod 6)` and `p + (p + 2) ≡ 0 (mod 12)`. (3, 5, 7) is the only prime triple p, p + 2, p + 4. *(theorem)*
10. Sophie Germain primes p > 3 (2p + 1 prime) have `p ≡ 5 (mod 6)`; their safe primes satisfy `2p + 1 ≡ 11 (mod 12)`. *(theorem)*
11. Every prime other than 2, 3, 5 divides the repunit `R_(p−1) = (10^(p−1) − 1)/9`. The decimal period of 1/p is the order of 10 mod p, a divisor of p − 1. *(theorem)*
12. If n² + 1 > 5 is prime then n is even and n² + 1 ends in 1 or 7. *(theorem)*
13. A palindromic prime with an even number of digits can only be 11: every even-length palindrome is divisible by 11. *(theorem)*

### 1B. Fermat, Wilson and friends

14. **Fermat's little theorem.** `a^p ≡ a (mod p)`, and `a^(p−1) ≡ 1 (mod p)` if p ∤ a. *(theorem: Fermat 1640, Euler 1736)*
15. **Wilson's theorem.** `(p − 1)! ≡ −1 (mod p)`; conversely n > 1 with (n − 1)! ≡ −1 (mod n) is prime. *(theorem: Wilson, Lagrange 1771)*
16. Half-factorial: `((p−1)/2)!² ≡ (−1)^((p+1)/2) (mod p)`. So for p ≡ 1 (mod 4), x = ((p−1)/2)! is an explicit square root of −1 mod p. *(theorem; checked p < 10⁴)*
17. **Wolstenholme.** `C(2p−1, p−1) ≡ 1 (mod p³)` and `H_(p−1) = 1 + 1/2 + … + 1/(p−1) ≡ 0 (mod p²)` for p ≥ 5. *(theorem: Wolstenholme 1862; first checked p < 2000)* Only 16843 and 2124679 satisfy it mod p⁴. *(literature)*
18. **Lucas's theorem.** `C(n, k) ≡ ∏ C(n_i, k_i) (mod p)` with n_i, k_i the base-p digits. **Kummer.** v_p(C(m+n, m)) is the number of carries when adding m and n in base p. **Legendre.** `v_p(n!) = Σ_i ⌊n/p^i⌋ = (n − s_p(n))/(p − 1)`, s_p = digit sum in base p. *(theorems)*
19. **Fibonacci.** `F_p ≡ (p/5) (mod p)` and `p | F_(p − (p/5))` for p ≠ 2, 5. *(theorem; checked p < 10⁴)* No Wall–Sun–Sun prime (p² | F_(p−(p/5))) is known. *(literature)*
20. **Wieferich primes** (`2^(p−1) ≡ 1 (mod p²)`): only 1093 and 3511 below 10⁷, and no others are known. *(checked below 10⁷; literature)*
21. **Wilson primes** (`(p−1)! ≡ −1 (mod p²)`): 5, 13, 563 below 1000, and none other below 2·10¹³. *(checked below 1000; literature: Costa–Gerbicz–Harvey 2014)*
22. **Clement (1949).** (n, n + 2) are twin primes iff `4((n−1)! + 1) + n ≡ 0 (mod n(n + 2))`. *(theorem; checked odd n < 1500)*
23. **Rowland (2008).** a(1) = 7, `a(n) = a(n−1) + gcd(n, a(n−1))`. Every difference a(n) − a(n−1) is 1 or a prime. *(theorem; checked n < 2·10⁵)*

### 1C. Squares mod p, reciprocity, quadratic forms

24. **Euler's criterion.** `a^((p−1)/2) ≡ (a/p) (mod p)`. There are (p − 1)/2 squares and (p − 1)/2 non-squares mod p. *(theorem)*
25. **Quadratic reciprocity.** `(p/q)(q/p) = (−1)^(((p−1)/2)((q−1)/2))` for odd p ≠ q. Supplements: `(−1/p) = (−1)^((p−1)/2)`, `(2/p) = (−1)^((p²−1)/8)`. *(theorem: Gauss 1796)*
26. Consequences: −1 is a square mod p iff p ≡ 1 (mod 4). 2 is a square iff p ≡ ±1 (mod 8). 3 is a square iff p ≡ ±1 (mod 12). −3 is a square iff p ≡ 1 (mod 3). 5 is a square iff p ≡ ±1 (mod 5). *(theorem; checked p < 10⁵)*
27. Least non-square n₀: `n₀(n₀ − 1) < p`. So n₀ < √p + 1. *(theorem)* Under GRH, n₀ < 2 ln² p. *(conditional: Bach 1990)*
28. **Sums of two squares.** An odd prime is a² + b² iff p ≡ 1 (mod 4), uniquely up to order and signs. *(theorem: Fermat, Euler 1749; checked p < 10⁴)* Gauss's formula finds a: for odd a ≡ 1 (mod 4), `C((p−1)/2, (p−1)/4) ≡ 2a (mod p)`. *(theorem: Gauss 1828; checked p < 10⁴)*
29. **Primes x² + ny².** *(theorems: Fermat, Euler, Gauss; see Cox, "Primes of the form x² + ny²"; all checked p < 10⁴)*
    ```
    p = x² + 2y²   ⟺  p = 2 or p ≡ 1, 3 (mod 8)
    p = x² + 3y²   ⟺  p = 3 or p ≡ 1 (mod 3)
    p = x² + 5y²   ⟺  p ≡ 1, 9 (mod 20)                         (p ≠ 2, 5)
    p = x² + 7y²   ⟺  p ≡ 1, 2, 4 (mod 7)                       (p odd, p ≠ 7)
    p = x² + 27y²  ⟺  p ≡ 1 (mod 3) and 2 is a cube mod p
    p = x² + 64y²  ⟺  p ≡ 1 (mod 8) and 2 is a fourth power mod p
    ```
30. **The non-residue trick.** An odd prime p ∤ D that divides x² − Dy² with gcd(x, y) = 1 has (D/p) = 1. So if (D/p) = −1, then p never divides x² − D. Examples: odd prime factors of n² + 1 are ≡ 1 (mod 4); of n² − 2 are ≡ ±1 (mod 8); of n² + n + 1 are 3 or ≡ 1 (mod 3). *(theorem)* This is the engine behind good prime-producing polynomials (facts 98–100).

### 1D. Primitive roots

31. The units mod p form a cyclic group, so primitive roots exist; there are `φ(p − 1)` of them. g is a primitive root iff `g^((p−1)/q) ≢ 1` for every prime q | p − 1. Primitive roots exist mod n iff n = 1, 2, 4, p^k or 2p^k. *(theorem: Gauss 1801)*
32. **Artin's conjecture.** 2 is a primitive root mod infinitely many p, with density `A = ∏_p (1 − 1/(p(p−1))) = 0.3739558136…`. *(conditional on GRH: Hooley 1967; the density checked for p < 10⁵, within 0.01)* At most two primes a fail Artin's conjecture. *(theorem: Heath-Brown 1986)*
33. **Consecutive primitive roots.** Every odd q > 169 has three consecutive primitive elements in F_q; for a prime q, three consecutive integers that are all primitive roots. *(theorem: Cohen–Oliveira e Silva–Trudgian 2015)*

### 1E. Special forms and their prime factors

34. If `2^n − 1` is prime then n is prime; more generally, if a^n − 1 is prime then a = 2 and n is prime. *(theorem)*
35. If `2^n + 1` is prime then n = 2^k, a Fermat number F_k = 2^(2^k) + 1. F₀ … F₄ = 3, 5, 17, 257, 65537 are prime; F₅ … F₁₄ are composite by Pépin's test. *(checked)* F₅ = 641 · 6700417 (Euler 1732). F₅ … F₃₂ are known to be composite. *(literature)*
36. Prime factors q of `2^p − 1` (p an odd prime) satisfy `q ≡ 1 (mod 2p)` and `q ≡ ±1 (mod 8)`. *(theorem; checked p < 60)*
37. Prime factors of F_n (n ≥ 2) satisfy `q ≡ 1 (mod 2^(n+2))`. *(theorem: Lucas 1878; checked n ≤ 6)*
38. If p ≡ 3 (mod 4) and q = 2p + 1 is prime, then `q | 2^p − 1`, so 2^p − 1 is composite for p > 3 (example: 23 | 2¹¹ − 1). *(theorem: Euler 1750; checked p < 10⁵)*
39. A repunit R_n = (10^n − 1)/9 can only be prime when n is prime. The proven repunit primes are R₂, R₁₉, R₂₃, R₃₁₇, R₁₀₃₁. *(theorem; literature)*
40. Mersenne exponents up to 1279: 2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279. *(checked with the Lucas–Lehmer test)* 52 Mersenne primes are known; the largest known prime is `2^136279841 − 1` (41,024,320 digits, October 2024). *(literature, as of 2025)*
41. A regular n-gon is constructible with ruler and compass iff n = 2^k times distinct Fermat primes. *(theorem: Gauss–Wantzel)*

### 1F. How primes spread over residue classes

42. **Dirichlet.** For gcd(a, m) = 1 there are infinitely many p ≡ a (mod m), and `π(x; m, a) ~ li(x)/φ(m)`. *(theorem: Dirichlet 1837, de la Vallée Poussin 1896)* Below 1.4·10¹¹ the last digits 1, 3, 7, 9 occur 1,421,581,391 / 1,421,587,551 / 1,421,582,979 / 1,421,574,235 times, equal to within 10⁻⁵. *(computed)*
43. **Chebyshev's bias.** π(x; 4, 3) ≥ π(x; 4, 1) for most x. The first exception is x = 26861 (Leech 1957). Up to 1.4·10¹¹, π(p; 4, 1) > π(p; 4, 3) holds at only 25,124,707 of the 5,686,326,158 primes (0.44 %), in short stretches. *(computed)* The logarithmic density of {x : π(x; 4, 3) > π(x; 4, 1)} is 0.9959. *(conditional on GRH and LI: Rubinstein–Sarnak 1994)*
44. Mod 3 the bias is stronger: `π(x; 3, 2) − π(x; 3, 1) ≥ 2` for every x in (10⁵, 1.4·10¹¹]. *(computed)* The first x with π(x; 3, 1) > π(x; 3, 2) is 608,981,813,029. *(literature: Bays–Hudson 1978; cited on [puzzle 83](https://www.primepuzzles.net/puzzles/puzz_083.htm))*
45. **Consecutive primes avoid repeating their last digit** (Lemke Oliver–Soundararajan 2016). Among the first 10⁸ primes, the last digits (row: p_n, column: p_(n+1)) are: *(computed)*
    ```
              1           3           7           9
    1   4,623,041   7,429,438   7,504,612   5,442,344      1 → 1 only 18.5 % instead of 25 %
    3   6,010,981   4,442,561   7,043,695   7,502,896
    7   6,373,982   6,755,195   4,439,355   7,431,870
    9   7,991,431   6,372,940   6,012,739   4,622,916
    ```
    The diagonal is depleted by about a quarter. The explanation (from the Hardy–Littlewood conjectures) is a second-order term of size about `ln ln x / ln x`, so the bias fades very slowly. *(heuristic)*

---

## 2. How many primes

46. **Euclid.** There are infinitely many primes. **Euler (1737).** `Σ 1/p` diverges. *(theorems)*
47. **Mertens (1874).** `Σ_(p≤x) 1/p = ln ln x + M + O(1/ln x)` with M = 0.2614972128…, `∏_(p≤x) (1 − 1/p) ~ e^(−γ)/ln x`, `Σ_(p≤x) ln p / p = ln x + O(1)`. *(theorem; M checked to 10⁻⁴ with primes < 2·10⁸)*
48. **Prime number theorem.** `π(x) ~ x/ln x ~ li(x)`, equivalently `θ(x) = Σ_(p≤x) ln p ~ x` and `ψ(x) = Σ_(p^k≤x) ln p ~ x`. Error: `π(x) − li(x) = O(x·exp(−c√(ln x)))`. *(theorem: Hadamard, de la Vallée Poussin 1896)* Consequences: `p_n ~ n ln n`, `Σ_(p≤x) p ~ x²/(2 ln x)`, `x# = e^((1+o(1))x)`.
49. **Riemann's explicit formula.** For x not a prime power, `ψ(x) = x − Σ_ρ x^ρ/ρ − ln 2π − ½·ln(1 − x⁻²)`. The sum runs over the non-trivial zeros ρ of ζ; the first is ½ + 14.134725…i. *(theorem: Riemann 1859, von Mangoldt 1895)* Every fluctuation of the primes around li(x) is a sum of waves x^ρ.
50. **Riemann hypothesis ⟺** `π(x) = li(x) + O(√x ln x)`. Under RH, `|π(x) − li(x)| < √x ln x / (8π)` for x ≥ 2657. *(conditional: Schoenfeld 1976)* RH has been verified up to height 3·10¹². *(literature: Platt–Trudgian 2021)*
51. **Riemann's R(x)** `= Σ_(n≥1) μ(n) li(x^(1/n))/n = 1 + Σ_(k≥1) (ln x)^k / (k·k!·ζ(k+1))` approximates π(x) much better than li(x): *(computed)*
    ```
    x        π(x)              li(x) − π(x)   R(x) − π(x)   d-digit primes
    10       4                      2.2           0.6        4
    10²      25                     5.1           0.7        21
    10³      168                    9.6           0.4        143
    10⁴      1,229                 17.1          −2.1        1,061
    10⁵      9,592                 37.8          −4.6        8,363
    10⁶      78,498               129.5          29.4        68,906
    10⁷      664,579              339.4          88.4        586,081
    10⁸      5,761,455            754.4          96.9        5,096,876
    10⁹      50,847,534         1,701.0         −78.6        45,086,079
    10¹⁰     455,052,511        3,103.6      −1,827.7        404,204,977
    10¹¹     4,118,054,813     11,587.6      −2,318.4        3,663,002,302
    1.4·10¹¹ 5,686,326,158
    ```
    A random d-digit number is prime with probability about 1/((d − ½)·ln 10): 4.07 % of the 11-digit numbers are prime.
52. **Explicit bounds.** `x/ln x < π(x)` for x ≥ 17 and `π(x) < 1.25506·x/ln x` for x > 1; `p_n > n ln n` for n ≥ 1; `p_n < n(ln n + ln ln n)` for n ≥ 6 *(theorem: Rosser–Schoenfeld 1962)*. `p_n > n(ln n + ln ln n − 1)` for n ≥ 2, and `p_n < n(ln n + ln ln n − 0.9484)` for n ≥ 39017 *(theorem: Dusart 1999)*.
53. **Littlewood (1914).** π(x) − li(x) changes sign infinitely often. *(theorem)* No sign change happens for 2 ≤ x ≤ 10¹⁹, and the first one is below 1.4·10³¹⁶. *(literature: Büthe 2018; Bays–Hudson 2000, Saouter–Demichel 2010)*
54. **Brun–Titchmarsh** (sharp form). `π(x + y) − π(x) ≤ 2y/ln y`, and `π(x; q, a) ≤ 2x/(φ(q) ln(x/q))` for x > q. *(theorem: Montgomery–Vaughan 1973)*
55. **Explicit primes in progressions.** `|θ(x; q, a) − x/φ(q)| < x/(160 ln x)` for 3 ≤ q ≤ 10⁵ and x ≥ 8·10⁹. *(theorem: Bennett–Martin–O'Bryant–Rechnitzer 2018)*
56. **Bombieri–Vinogradov.** On average over moduli q ≤ √x/(ln x)^B, primes are as well distributed in progressions as GRH predicts. The Elliott–Halberstam conjecture extends this to q ≤ x^(1−ε). *(theorem 1965; conjecture)*
57. **Short intervals.** `π(x + h) − π(x) ~ h/ln x` holds for h ≥ x^(7/12+ε) *(theorem: Huxley 1972)*. It fails for h = (ln x)^λ with λ > 1: the count oscillates by a constant factor *(theorem: Maier 1985)*.
58. **Hardy–Littlewood's second conjecture** `π(x + y) ≤ π(x) + π(y)` is incompatible with the prime k-tuple conjecture (fact 83). *(theorem: Hensley–Richards 1974)*
59. **Counting algorithms.** Legendre: `π(x) − π(√x) + 1 = Σ_(d | ∏_(p≤√x) p) μ(d)·⌊x/d⌋`. Combinatorial methods compute π(x) in about x^(2/3) steps (Lagarias–Miller–Odlyzko 1985, Deléglise–Rivat 1996), analytic ones in x^(1/2+ε) (Lagarias–Odlyzko 1987). No method polynomial in ln x is known. *(theorems)*
60. **Goldbach-type results.** Every odd n > 5 is a sum of three primes *(theorem: Helfgott 2013)*. Every large even n is p + q or p + q₁q₂ *(theorem: Chen 1973)*. Every integer > 6 is a sum of distinct primes *(theorem: Richert 1949)*. Even Goldbach has been verified up to 4·10¹⁸ *(literature: Oliveira e Silva, Herzog, Pardi 2014)*.
61. **Sylvester (1892).** A product of k consecutive integers, all greater than k, has a prime factor > k. *(theorem)*

---

## 3. Gaps between consecutive primes

62. g₁ = 1 (from 2 to 3); every later gap is even. The average gap near x is ln x, since `Σ_(p_n≤x) g_n ~ x`. *(theorem)*
63. **Arbitrarily long gaps.** `k# + 2, k# + 3, …, k# + k` (or n! + 2, …, n! + n) are all composite. How long such a run can be made with small primes is the Jacobsthal function j(m), the largest gap between integers coprime to m. *(theorem)*
64. **Bertrand's postulate.** There is a prime in (n, 2n] for every n ≥ 1, so `p_(n+1) < 2p_n`. *(theorem: Chebyshev 1852)* Sharper: a prime in (n, 6n/5] for n ≥ 25 *(Nagura 1952)*, and in `(x, x(1 + 1/(2 ln² x))]` for x ≥ 3275 *(Dusart 2010)*. Every x from 3275 to 1.4·10¹¹ respects Dusart's bound, and 3275 is sharp: from x = 3274 the next prime, 3299, is just out of reach *(checked)*.
65. **Baker–Harman–Pintz (2001).** For large x there is a prime in `[x − x^0.525, x]`. *(theorem)* **Ingham (1937).** There is a prime between n³ and (n + 1)³ for large n. *(theorem)*
66. Under RH, `g_n = O(√p_n · ln p_n)`. *(conditional: Cramér 1920)*
67. **Legendre's conjecture** (a prime between n² and (n + 1)²) and **Andrica's conjecture** (`√p_(n+1) − √p_n < 1`) are open. Below 1.4·10¹¹ every gap after p ≥ 11 satisfies `g < 1.32·√p`. The largest ratio, 1.317, is at 113 → 127. Since g < 2√p implies both conjectures, they hold in this range. Andrica's maximum is at 7 → 11: √11 − √7 = 0.6708. *(computed)*
68. **Cramér's conjecture.** `g_n = O(ln² p_n)` *(conjecture: Cramér 1936)*. Granville's refinement of the model suggests `limsup g_n / ln² p_n ≥ 2e^(−γ) ≈ 1.1229` *(heuristic: Granville 1995)*. Firoozbakht's conjecture would give `g_n < ln² p_n − ln p_n` for n > 4 *(Kourbatov)*. The largest known ratio g/ln² p (ignoring p = 2, 3, 7) is 0.9206386, at p = 1693182318746371 *(literature)*. Below 1.4·10¹¹ the largest is 0.7953, for the gap 456 after 25056082087 *(computed)*.
69. **Small gaps.** `liminf g_n / ln p_n = 0` *(Goldston–Pintz–Yıldırım 2005)*. `liminf g_n < 7·10⁷` *(Zhang 2013)*, ≤ 600 *(Maynard 2013)*, ≤ 246 *(Polymath 2014)*, and ≤ 6 under the generalised Elliott–Halberstam conjecture. For every m, `liminf (p_(n+m) − p_n) < ∞`, with a bound of order m³e^(4m) *(Maynard–Tao)*. *(theorems)*
70. **Large gaps.** `limsup g_n / ln p_n = ∞` *(Westzynthius 1931)*. Infinitely often `g_n > c · ln p ln ln p ln ln ln ln p / (ln ln ln p)` *(Ford–Green–Konyagin–Maynard–Tao 2018, improving Rankin 1938)*. *(theorems)*
71. **Erdős–Turán (1948).** The gaps rise (g_(n+1) > g_n) infinitely often and fall infinitely often. *(theorem)*
72. **Maximal gaps** (a gap larger than every earlier one), as gap:prime-before-it, below 1.4·10¹¹: *(computed; equal to OEIS A005250 / A002386)*
    ```
    1:2  2:3  4:7  6:23  8:89  14:113  18:523  20:887  22:1129  34:1327  36:9551  44:15683  52:19609
    72:31397  86:155921  96:360653  112:370261  114:492113  118:1349533  132:1357201  148:2010733
    154:4652353  180:17051707  210:20831323  220:47326693  222:122164747  234:189695659  248:191912783
    250:387096133  282:436273009  288:1294268491  292:1453168141  320:2300942549  336:3842610773
    354:4302407359  382:10726904659  384:20678048297  394:22367084959  456:25056082087
    464:42652618343  468:127976334671
    ```
    85 maximal gaps are known. The largest is 1854, after 101412319996363309069 (R. Smith, 2026). The largest known merit g/ln p is 42.52, for a gap of 59724 after a 611-digit prime (2026). *(literature: Wikipedia "Prime gap", archived September 2026)* The number of maximal gaps below p_n is conjectured to be about 2 ln n. *(heuristic)*
73. **Merit.** Below 1.4·10¹¹ the largest merit g/ln p is 19.04, for the gap 456 after 25056082087. *(computed)*
74. **Which gaps are common.** 6 is the most common gap at every checkpoint from 10⁶ to 1.4·10¹¹. Below 1.4·10¹¹: *(computed)*
    ```
    gap      2           4           6           8           10          12          14          16          18          20
    count    305,560,420 305,552,719 557,127,352 253,193,700 329,666,173 437,953,710 246,234,893 188,935,421 350,765,217 200,151,778
    gap      22          24          26          28          30          36          42          48          60
    count    170,746,147 267,867,928 131,044,274 143,130,077 268,847,861 155,996,782 139,311,351  90,840,234  69,310,603
    ```
    Multiples of 6 are favoured, and 30 already beats 24. Gaps 2 and 4 are almost equally common: a gap of 4 after p > 3 forces p ≡ 1 (mod 6), so p + 2 is divisible by 3 and nothing is lost. 218 different gap sizes occur, the largest 468. *(computed)*
75. **Why.** By Hardy–Littlewood, the number of prime pairs (p, p + 2k), consecutive or not, is about `2C₂ · ∏_(q | k, q > 2) (q − 1)/(q − 2) · ∫₂ˣ dt/ln² t`. A factor 3 in the difference doubles the count, a factor 5 multiplies it by 4/3. *(conjecture)* Twin primes against this formula: *(computed)*
    ```
    x         twin pairs ≤ x    2C₂ ∫₂ˣ dt/ln²t    difference
    10⁶       8,169             8,248              +79
    10⁸       440,312           440,368            +56
    10¹⁰      27,412,679        27,411,417         −1,263
    10¹¹      224,376,048       224,368,865        −7,183
    1.4·10¹¹  305,560,420
    ```
76. **Jumping champions.** The most common gap below x is conjectured to be 6 until about 1.7·10³⁵, then 30, then 210, …: the primorials. *(heuristic: Odlyzko–Rubinstein–Wolf 1999)*
77. **Poisson law.** If the Hardy–Littlewood conjectures hold uniformly, the number of primes in (x, x + λ ln x] for random x is asymptotically Poisson with mean λ. *(conditional: Gallagher 1976)*
78. **Capelle's floor identity.** For consecutive primes r < q < p, `⌊pq/r⌋ = p + q − r` iff `(p − r)(q − r) < r`. It fails exactly for r ∈ {2, 3, 5, 7, 11, 13, 19, 23, 31, 47, 83, 89, 113, 199, 1327} below 2·10⁸. *(site: [conjecture 53](https://www.primepuzzles.net/conjectures/conj_053.htm); checked below 2·10⁸)*
79. **Gaps look exponential.** In Cramér's model (each n > 2 is prime with probability 1/ln n, independently), the gap after x exceeds t·ln x with probability about e^(−t). Real gaps follow this law with the Hardy–Littlewood corrections of facts 74–75. *(heuristic: Cramér 1936; conditional form: fact 77)*

---

## 4. Patterns in consecutive primes

80. **Admissible tuples.** A set H = {h₁ < … < h_k} is admissible if, for every prime q, the numbers h_i miss at least one residue class mod q. Only q ≤ k need checking. Only admissible patterns can be prime infinitely often (example: {0, 2, 4} is not admissible mod 3). *(theorem)*
81. **Least width s(k) of an admissible k-tuple** (h_k − h₁): *(checked by exhaustive search)*
    ```
    k     2  3  4  5   6   7   8   9   10  11  12  13  14  15  16  17  18
    s(k)  2  6  8  12  16  20  26  30  32  36  42  48  50  56  60  66  70
    ```
82. **Tightest clusters actually found** (least span of k consecutive primes, first prime > k): *(computed below 1.4·10¹¹)*
    ```
    k              2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19
    least span     2   6   8   12  16  20  26  30  32  40  44  50  54  56  62  68  72  76
    first prime    3   5   5   7   7   11  11  11  11  13  17  17  17  17  17  41  37  37
    ```
    For k = 11–14 and 16–18 the minimum s(k) is not reached below 1.4·10¹¹. The small primes are exceptionally dense: the 15 primes 17 … 73 have the least possible width 56, and the 19 primes 37 … 113 the least width 76.
83. **Hardy–Littlewood k-tuple conjecture.** For admissible H, `#{n ≤ x : n + h_i all prime} ~ 𝔖(H) ∫₂ˣ dt/ln^k t` with `𝔖(H) = ∏_q (1 − ν_H(q)/q)(1 − 1/q)^(−k)`, where ν_H(q) is the number of classes mod q that H occupies. For twins, `𝔖 = 2C₂`, `C₂ = ∏_(q>2) (1 − 1/(q − 1)²) = 0.6601618158…`. *(conjecture: Hardy–Littlewood 1923; C₂ checked to 10⁻⁷)* 𝔖 is unchanged when the pattern is mirrored (h → h_k − h), so a pattern and its mirror image are equally common. Fact 88 confirms this.
84. **Brun (1919).** `Σ 1/p` over twin primes converges, to Brun's constant B₂ = 1.902160583… (an estimate based on the Hardy–Littlewood tail, not a proof). *(theorem; literature)* The largest known twin primes are `2996863034895 · 2^1290000 ± 1` (388,342 digits, 2016). *(literature, as of 2025)*
85. **Consecutive primes in arithmetic progression (CPAP).** Primes in arithmetic progression, consecutive or not: if p, p + d, …, p + (k − 1)d are all prime and p > k, then every prime q ≤ k divides d, so `k# | d`. (If q ∤ d, the k terms cover every class mod q.) *(theorem)* First CPAP-k: k = 3 at 3 (d = 2), 4 at 251 (d = 6), 5 at 9843019 (d = 30), 6 at 121174811 (d = 30). *(computed; these are the minimal ones, and the minimal CPAP-k is known only for k ≤ 6)* Only two CPAP-10 are known. The first (1998, M. Toplic, in the CP10 project of Dubner, Forbes, Lygeros, Mizony, Zimmermann) starts at the 93-digit number `100996972469714247637786655587969840329509324689190041803603417758904341703348882159067229719` with d = 210. *(checked: the 10 terms are PRP and consecutive)* A CPAP-11 needs d to be a multiple of 2310 and at least 23090 composites between its terms; Dubner and Zimmermann estimate it at least 10¹² times harder than a CPAP-10. *(literature)*
86. **Green–Tao (2008).** The primes contain arbitrarily long arithmetic progressions (not consecutive). *(theorem)* The longest known has 27 terms: `224584605939537911 + 81292139 · 23# · n`, n = 0 … 26 (PrimeGrid, R. Gahan, 2019). *(checked: all 27 are prime, below 2⁶⁴)*
87. **Euler–Rodríguez chains** ([puzzle 151](https://www.primepuzzles.net/puzzles/puzz_151.htm)): k consecutive primes `A + x(x + 1)`, x = 0 … k − 1, so the gaps are 2, 4, 6, …, 2(k − 1). Least A for a chain of exactly k primes:
    ```
    k   2  3  4   5    6      7       8       9          10           11              12               13                 14
    A   3  5  17  347  13901  665111  128981  113575727  2426256797   137168442221    4656625081181    101951758179851    484511389338941
    ```
    For k ≤ 11 these values are minimal by our exhaustive pass to 1.4·10¹¹ *(computed)*. k = 12–14 are J. K. Andersen's (2003); we verified the chains *(site; checked)*. **No chain of 15 is known.** Andersen searched to 10¹⁷ *(site)*. A necessary condition when A > 2k − 1: `1 − 4A` must be a quadratic non-residue mod every odd prime q ≤ 2k − 1. (Otherwise x² + x + A ≡ 0 (mod q) has a root x₀ ≤ (q − 1)/2 < k, and q divides the value at x₀.) Once that holds, q never divides any value of x² + x + A. *(theorem, short proof above)*
88. **Mirror chains** (gaps 2(k − 1), …, 4, 2, i.e. primes `L − x(x + 1)`): first chains of exactly k primes are at 3, 7, 31, 1979, 41203, 752251, 5647457, 32465047, 2636566333 (k = 2 … 10), 245333213 (k = 11) and 27797667517 (k = 12). *(computed)* Forward and mirrored chains below 10¹⁰ are equally common, as fact 83 predicts: *(computed)*
    ```
    k          5        6       7      8      9     10    11
    forward    43,938   5,154   767    109    8     3     0
    mirrored   44,418   5,364   729    118    12    1     1
    ```
    So the mirrored 12-chain at 2.8·10¹⁰, which comes 170 times earlier than the first forward one, is most likely luck and not a bias.
89. **Runs of consecutive primes on a polynomial.** p_n, …, p_(n+R−1) lie on a polynomial of degree d in the index exactly when their (d + 1)-th differences vanish. Any d + 1 primes fit a degree-d polynomial (Lagrange), so a run longer than d + 1 is a coincidence. First run of length R below 1.4·10¹¹: *(computed)*
    ```
    degree 1 (CPAP)      R=3: 3   4: 251   5: 9843019   6: 121174811                    (no run of 7)
    degree 2             R=4: 17  5: 347   6: 2903   7: 15373   8: 128981   9: 19641263
                         10 and 11: 245333213          12: 27797667517                   (longest: 12)
    degree 3             the same first runs as degree 2 from R = 5 on
    degree 4             R=6: 41   7: 8081, then the same first runs as degree 2 from R = 8 on
    ```
    The long quadratic runs have gaps falling by 2: 19641263 has gaps 20, 18, …, 6, 245333213 has 20, …, 2, and 27797667517 has 22, 20, …, 2. So the 12 consecutive primes from 27797667517 to 27797667649 are exactly `27797667649 − x(x + 1)`, x = 0 … 11. Below 1.4·10¹¹, cubics and quartics never fit longer runs than quadratics.
90. **Strings in one residue class.** For any m, any a coprime to m and any k, there are k consecutive primes all ≡ a (mod m). *(theorem: Shiu 2000)* Longest strings below 1.4·10¹¹: *(computed)*
    ```
    mod 3:  25 primes ≡ 1 from 43553959717;    25 primes ≡ 2 from 21549657539
    mod 4:  26 primes ≡ 1 from 84676452781;    32 primes ≡ 3 from 113391385603
    mod 10: 11 ending in 1 from 14293856441;   12 ending in 3 from 131145172583;
            12 ending in 7 from 73712513057;   12 ending in 9 from 130426565719
    ```

---

## 5. Formulas that produce primes

91. **No polynomial gives only primes.** If f ∈ ℤ[x] is non-constant and f(a) = p, then `p | f(a + kp)` for every k. *(theorem: Goldbach 1752)*
92. **Euler's polynomial.** `n² + n + 41` is prime for n = 0 … 39. `n² − 79n + 1601` gives 80 primes for n = 0 … 79, each twice. No prime below 41 divides any value: n² + n is even, and the discriminant −163 is a non-square mod every odd prime < 41. *(theorem; checked)*
93. **Rabinowitsch (1913).** `n² + n + A` is prime for n = 0 … A − 2 iff `4A − 1 ∈ {3, 7, 11, 19, 43, 67, 163}`, i.e. Q(√(1 − 4A)) has class number 1. By Heegner, Baker and Stark there are no others, so A ∈ {1, 2, 3, 5, 11, 17, 41} and 40 is the best possible for this family. *(theorem; checked A ≤ 10⁴)*
94. **Bunyakovsky conjecture.** An irreducible f ∈ ℤ[x] with positive leading coefficient and no fixed prime divisor is prime infinitely often. Dirichlet's theorem is the linear case. *(conjecture)*
95. **Bateman–Horn conjecture.** `#{n ≤ x : f₁(n), …, f_m(n) all prime} ~ C/(∏ deg f_i) · ∫₂ˣ dt/ln^m t`, with `C = ∏_q (1 − N(q)/q)(1 − 1/q)^(−m)` and N(q) the number of roots of f₁⋯f_m mod q. A polynomial whose values avoid small primes (N(q) = 0) has a large C. *(conjecture: Bateman–Horn 1962)*
96. **Partial results.** Infinitely many n have n² + 1 with at most two prime factors *(Iwaniec 1978)*. Infinitely many primes have the form a² + b⁴ *(Friedlander–Iwaniec 1998)* and x³ + 2y³ *(Heath-Brown 2001)*. Infinitely many primes have no digit 7 *(Maynard 2019)*. *(theorems)*
97. **Polynomial records from the site.** For |f(n)| prime at consecutive n: `36n² − 810n + 2753` gives 45 distinct primes (Fung–Ruby). Among cubics, `3n³ − 183n² + 3318n − 18757` gives 47 primes, 43 distinct (S. M. Ruiz). *(site: [puzzle 232](https://www.primepuzzles.net/puzzles/puzz_232.htm))* Our own polynomial records are in section 9.
98. **Why good polynomials work.** Choose the constant so that small primes cannot divide any value: the discriminant must be a non-square mod as many small primes as possible (facts 30, 92, 95). This is the trick behind Euler's 41 and several of our records. *(ours: [IDEAS.txt](IDEAS.txt))*
99. **A list as a formula.** Any k consecutive primes are the values of a unique polynomial of degree ≤ k − 1 at 0 … k − 1 (Lagrange). And `Σ_j e_j·((m + j) mod k)` can produce any list of k values, since the circulant system is invertible. So "a formula for k primes" only means something when the formula is much shorter than the list. *(theorem; ours: [1097](results/1097.md))*
100. **Exact formulas** (all of them circular: they encode the primes they produce). *(checked as stated)*
    ```
    Willans (1964)   p_n = 1 + Σ_(i=1)^(2^n) ⌊ (n / Σ_(j=1)^i ⌊cos²(π((j−1)! + 1)/j)⌋)^(1/n) ⌋     (Wilson inside; checked n ≤ 8)
                     equivalently p_n = 1 + #{ m ≤ 2^n : π(m) < n }
    Mills (1947)     ⌊A^(3^n)⌋ is prime for all n; least A = 1.30637788386308069046…  (computed from 2, 11, 1361, 2521008887,
                     16022236204009818131831320183, …; the 29-digit term is PRP; minimality assumes RH)
    Fridman et al.   f₁ = Σ_n (p_n − 1)/∏_(i<n) p_i = 2.92005097731613471209…,  f_(n+1) = ⌊f_n⌋·(f_n − ⌊f_n⌋ + 1),  ⌊f_n⌋ = p_n
    (2019)           (all primes in order; it works because of Bertrand's postulate; checked for the first 60 primes)
    Wright (1951)    ⌊2^2^…^2^ω⌋ (n twos) is prime, ω = 1.9287800…                         (literature)
    JSWW (1976)      a polynomial of degree 25 in 26 variables whose positive values are exactly the primes
    ```

---

## 6. Primality tests and special primes

101. **Lucas–Lehmer.** For an odd prime p, M_p = 2^p − 1 is prime iff `s_(p−2) ≡ 0 (mod M_p)`, where s₀ = 4 and s_(k+1) = s_k² − 2. **Pépin (1877).** F_n is prime iff `3^((F_n−1)/2) ≡ −1 (mod F_n)`. **Proth (1878).** N = k·2^n + 1 (k odd, k < 2^n) is prime iff `a^((N−1)/2) ≡ −1 (mod N)` for some a. *(theorems)*
102. **AKS (2002).** Primality can be decided in deterministic polynomial time. *(theorem)* BPSW has no known counterexample and is correct below 2⁶⁴, which is why this repo calls sympy.isprime a proof only below 2⁶⁴. *(literature)*
103. **Pseudoprimes.** The least base-2 Fermat pseudoprime is 341 = 11 · 31. Korselt: a composite n satisfies a^n ≡ a for all a iff n is squarefree and p − 1 | n − 1 for every p | n (Carmichael numbers). 561 = 3 · 11 · 17 is the least, and there are infinitely many. *(theorems: Korselt 1899; Alford–Granville–Pomerance 1994)*

---

## 7. Constants

104. Constants used above. *(the first three checked with our sieves; the others literature or computed as stated)*
    ```
    twin prime constant   C₂ = ∏_(p>2) (1 − 1/(p−1)²)        = 0.6601618158…
    Artin's constant      A  = ∏_p (1 − 1/(p(p−1)))           = 0.3739558136…
    Meissel–Mertens       M  = lim (Σ_(p≤x) 1/p − ln ln x)     = 0.2614972128…
    Brun's constant       B₂ = Σ over twin pairs (1/p + 1/(p+2)) ≈ 1.902160583   (estimate)
    Mills' constant       A  = 1.30637788386308069046…          (computed; least, assuming RH)
    Fridman's constant    f₁ = 2.92005097731613471209…          (computed)
    e^γ = 1.7810724179…,  2e^(−γ) = 1.1229189671…
    ```

---

## 8. Record lengths at a glance

| Structure | Longest known | Where |
|---|---|---|
| Consecutive primes in arithmetic progression | 10 (two known, 1998 and 2008) | fact 85 |
| Primes in arithmetic progression (not consecutive) | 27 (2019) | fact 86 |
| Euler–Rodríguez chain, consecutive primes A + x(x + 1) | 14 (2003); 15 open, none below 10¹⁷ | fact 87 |
| Consecutive primes on a quadratic, below 1.4·10¹¹ | 12 (27797667649 − x(x + 1)) | facts 88, 89 |
| Euler-type polynomial n² + n + A, prime for n = 0 … A − 2 | 40 (A = 41), provably best | fact 93 |
| Tightest prime cluster below 1.4·10¹¹ with least possible width | 15 primes (17 … 73), 19 primes (37 … 113) | fact 82 |
| Consecutive primes ≡ 3 (mod 4), below 1.4·10¹¹ | 32 | fact 90 |
| Run of composites (maximal gap), below 1.4·10¹¹ | 467, after 127976334671 | fact 72 |
| Largest known maximal gap | 1854, after 101412319996363309069 (2026) | fact 72 |
| Largest known prime | 2^136279841 − 1, 41,024,320 digits (2024) | fact 40 |
| Largest known twin primes | 2996863034895 · 2^1290000 ± 1 (2016) | fact 84 |

---

## 9. Facts from this repository

(In progress: the facts proved or used in our write-ups are being collected and checked.)

## 10. Facts from the site and its references

(In progress: the facts stated on the archived pages and in the reference documents are being collected and checked.)

---

## 11. Missing parts: what nobody knows yet

Each line gives the open claim, its formula, and the best known result.

- **Twin primes.** `g_n = 2` infinitely often? Known: g_n ≤ 246 infinitely often (fact 69).
- **Polignac.** Every even gap occurs infinitely often? Not known for any single gap.
- **Prime k-tuples and Bateman–Horn** (facts 83, 95). Not one pattern with k ≥ 2 is known to occur infinitely often. Polynomials: Dirichlet (degree 1) is the only proved case.
- **Landau's four problems (1912):** Goldbach, twin primes, Legendre, and primes n² + 1. All open. Partial: three primes (Helfgott), Chen, Iwaniec (facts 60, 96).
- **Gaps from above.** Cramér `g = O(ln² p)`, Legendre, Andrica, Oppermann, Firoozbakht. Best proved: g ≤ p^0.525 (fact 65), O(√p ln p) under RH.
- **Riemann hypothesis**, which would give the best possible error in the prime number theorem (fact 50).
- **Artin's conjecture**, proved only under GRH (fact 32).
- **Special primes.** Are there infinitely many Mersenne, Sophie Germain, Wieferich, Wilson or Wall–Sun–Sun primes, or primes n! ± 1? Is F₄ the last Fermat prime?
- **Consecutive primes in arithmetic progression.** Not even "infinitely many CPAP-3" is proved, though it follows from the k-tuple conjecture. CPAP-11: none known. The least CPAP-7: unknown (fact 85).
- **Euler–Rodríguez chains.** No chain of 15 is known (fact 87).
- **A formula for p_n.** No formula that is not circular (fact 100). No algorithm for π(x) polynomial in ln x (fact 59). No deterministic polynomial-time way to *construct* a prime of d digits is known (the Polymath4 project did not find one).
- **Bias questions.** The explanations for Chebyshev's bias (fact 43) and the last-digit bias (fact 45) are conditional on GRH or on Hardy–Littlewood.
- **Our own open attempts** are listed in [writeup/INDEX.md](writeup/INDEX.md).

---

## 12. Verify it yourself

- [results/code/facts_check.py](results/code/facts_check.py) checks every *checked* fact and the small *computed* ones in about 15 seconds. It uses its own numpy sieve to 2·10⁸, sympy (exact below 2⁶⁴) and gmpy2. Run: `.venv/bin/python results/code/facts_check.py`.
- [results/code/facts.c](results/code/facts.c) is one pass over all primes below 1.4·10¹¹, about 7 minutes. It produces the counts, maximal gaps, gap frequencies, polynomial runs, chains, strings, clusters and biases of facts 42–45, 51, 67–75, 82, 85, 87, 89 and 90. Run: `gcc -O3 -o facts results/code/facts.c -lm && ./facts 140000000000`.
- [results/code/facts_chains.c](results/code/facts_chains.c) counts forward and mirrored Euler–Rodríguez chains (fact 88). Run: `./facts_chains 10000000000`, about 13 seconds.
