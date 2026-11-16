# Literature check: ten statements (Mersenne heuristics, wheel patterns, biases)

Checked 2026-11-16, about 60 minutes of web search. Rules: only sources I opened and read, or
search-engine snippets clearly marked "(snippet)". Access limits: the AMS PDF of Wagstaff 1983
returned HTTP 403, so its content is known only from snippets and from a later paper by Wagstaff
that quotes it. oeis.org returned 403, so I used the local OEIS copy in `references/oeis/`.
PDFs that WebFetch could not read were saved by the tool and converted with pypdf in the scratchpad.

---

## 1. Wagstaff's constant refined by p mod 3 / p mod 12; counts of Mersenne exponents mod 12

**Verdict: PARTLY KNOWN.** The mod-4 refinement is classical. I found no refinement by p mod 3
or p mod 12, and no published counts of Mersenne exponents mod 12.

What is known:
- Gillies (1964), sharpened by J. R. Ehrman, "The number of prime divisors of certain Mersenne
  numbers", Math. Comp. 21 (1967) 700-704. In S. S. Wagstaff Jr., "Divisors of Mersenne numbers",
  Math. Comp. 40 (1983) 385-397, the sharpened third consequence has a = 2 for p = 3 mod 4 and
  a = 6 for p = 1 mod 4 (snippet). Wagstaff gives the allowed k as 3,4,7,8,11,12,... for
  p = 1 mod 4 and 1,4,5,8,9,12,... for p = 3 mod 4, and says M_p "has a better chance of being prime
  if p = 1 (mod 4)". His count was 16 against 10 among the then-known exponents, which he called
  "suggestive but not statistically significant" (snippets of the AMS PDF; not read in full).
- Wagstaff's heuristic averages away every finer local factor. S. S. Wagstaff Jr.,
  P. L. Montgomery and S. Nahm, "The period of the Bell numbers modulo a prime", Math. Comp. 79
  (2010) 1793-1800, §3, p. 1799 (read). It quotes Bateman-Horn for p and 2kp+1 both prime:
  2·C2·f(2k)/log(2kp), with f(n) = prod over odd q | n of (q-1)/(q-2). It then says "Just as in the
  heuristic argument on page 386 of [Wagstaff 1983] we may replace C2 f(2k) by 1". So the
  3 | k structure, and with it any p mod 3 effect, is explicitly averaged out in Wagstaff 1983.
- Shanks and Kravitz (item 2) give per-k densities that include the factor over odd q | k (q = 3
  included). These are averaged over all p, not split by p mod 3.
- The Prime Pages heuristic page (t5k.org/mersenne/heuristic.html, read) uses only a = 2/6 by
  p mod 4. Its prose explanation of which k is excluded looks garbled. The k lists above are the
  reliable version: the smallest admissible divisor is 6p+1 when p = 1 mod 4, hence a = 6.
- Counts mod 4: M. Wolf, "Computer experiments with Mersenne primes", arXiv:1112.2412 (2011), p. 1
  (read): "Of these 47 known M_n ... 27 with p mod 4 = 1 and 19 with p mod 4 = 3". OEIS A112648 and
  A112649 split the Mersenne exponents by p mod 4. OEIS A000043 has no residue statistics.
- J. Domínguez, "Divisor structure of p-1 in Mersenne prime exponents", arXiv:2603.08994 (2026).
  It studies a p-1 divisor statistic, not residue classes, and says (§5) "No modification of the
  classical Wagstaff asymptotic law is proposed."

Searches with no hit: "Wagstaff heuristic Mersenne mod 12"; "Mersenne prime exponents p mod 3
probability 6kp+1 divisible by 3"; mersenneforum site: queries on "mod 12", "mod 3", "1 mod 4 /
3 mod 4 ratio"; "Mersenne exponents modulo 3 / mod 6 / mod 12 residue table"; local OEIS names
grep "Mersenne exponent ... mod/congruent".

## 2. Per-exponent local factor from which 2kp+1 (small k) are prime

**Verdict: PARTLY KNOWN** (the per-k ingredients); **NOT FOUND** (a per-p singular-series
correction to P(M_p prime)).

Known ingredients:
- D. Shanks and S. Kravitz, "On the distribution of Mersenne divisors", Math. Comp. 21 (1967) 97-101
  (text read on t5k.org/mersenne/LukeMirror/lit/lit_048s.htm; the formulas there are images). They
  assume "1/k is the probability that 2 is a (2k)ic residue of 2kn+1". They define f_k(x) = number
  of p < x such that M_p has the prime divisor 2kp+1. Conjecture (3) carries a factor that depends
  on k mod 4 and a product over the odd primes q | k. They check the k = 3 case by Landau's prime
  ideal theorem. Table 1 compares counts for k < 60 up to 10^6.
- Wagstaff 1983, p. 386, as quoted in Wagstaff-Montgomery-Nahm 2010, p. 1799: for a fixed p, the
  expected number of prime factors in (A, B] is the sum over k of
  [2·C2·f(2k)/log(2kp)] · [prob(divides)], and C2·f(2k) is then replaced by 1. The k-dependence is
  written down and then averaged; no per-p conditioning on which 2kp+1 are actually prime.
- Wagstaff 1983 (snippet): among the 4783 primes p = 1 mod 4 below 10^5, 6p+1 is prime for 1037,
  and 6p+1 divides M_p for 350 of them (0.34, against 1/3).
- GIMPS (mersenne.org/various/math.php, read): the chance of primality after trial factoring is
  based only on the factoring depth. There is no adjustment by p's residue or by which 2kp+1 are
  prime.
- Domínguez 2026 (item 1) is the only per-exponent correction I found. It uses a different
  mechanism (the divisors of p-1) and is empirical.

Searches with no hit: "probability Mersenne prime given 2kp+1 prime small k conditional heuristic
per exponent"; "Mersenne heuristic singular series exponent local factor"; mersenneforum
"exponent weight / sieve weight"; "Mersenne exponent 2p+1 composite more likely prime".

## 3. Frequencies of the 256 primality patterns of the wheel-24 "prime byte", from HL + inclusion-exclusion

**Verdict: PARTLY KNOWN** (the method is classical); **NOT FOUND** (the wheel-24 byte or any
2^8-pattern table).

- R. P. Brent, "The distribution of small gaps between successive primes", Math. Comp. 28 (1974)
  315-324. It applies inclusion-exclusion to Hardy-Littlewood k-tuple constants to predict exact
  prime/composite patterns: "p prime, p+2r the next prime", with constants A_{r,k} (abstract and
  secondary description read).
- R. J. Lemke Oliver and K. Soundararajan, "Unexpected biases in the distribution of consecutive
  primes", PNAS 113 (2016) E4446-E4454 (arXiv:1603.03720). It predicts the frequencies of
  consecutive-prime residue patterns mod q from Hardy-Littlewood, with slowly decaying secondary
  terms.
- L. Lacasa, B. Luque, I. Gómez and O. Miramontes, "On a dynamical approach to some prime number
  sequences", Entropy 20 (2018) 131 (arXiv:1802.08349). They "estimate the densities of gap block
  residues via Hardy-Littlewood k-tuple conjecture" (abstract).
- A. Granville and A. Lumley, "Primes in short intervals: heuristics and calculations",
  arXiv:2009.05000. A modified Cramér model built on small-prime sieving, the same local-factor
  idea; I read only secondary descriptions.
- The pzktupel.de page "Patterns of prime k-tuplets & the Hardy-Littlewood constants" (read) covers
  all-prime patterns only, not mixed prime/composite patterns.

Searches with no hit: "frequency of prime patterns in wheel residues mod 30 byte 256 patterns
Hardy-Littlewood inclusion-exclusion"; "prime constellations distribution of patterns".

## 4. Entropy / mutual information of the prime indicator predicted by HL; what ML can extract

**Verdict: PARTLY KNOWN.** Entropies of prime-related symbol sequences have been studied, but I
found no HL-based Shannon entropy or mutual information of the prime indicator, wheel-reduced or
not.

- Lacasa et al. 2018 (item 3): Rényi block entropies of residues of primes mod k ("maximally
  chaotic") and of gap residues, with HL used for the gap-residue block probabilities.
- M. Blank, "Are prime numbers and quadratic residues random?", arXiv:2403.04490 (v2 2025, "to
  appear in Discrete and Continuous Dynamical Systems"), PDF read.
  - It defines an "information entropy" h_info = lim (1/n) log2 L(n, N), where L counts distinct
    words. This is a word-count (topological-type) quantity, not a frequency-weighted Shannon
    entropy.
  - Thm 5.1: h_info(prime indicator) <= log2((1+√5)/2) ≈ 0.694 from "no 11" (Lemma 5.3). Part (b)
    claims equality under Hardy-Littlewood.
  - **Skeptical note (my own check, not literature):** part (b) cannot be right.
    - Two 1s an odd distance apart would need an even prime, so the bound is already <= 1/2 bit.
    - Under HL the words that occur at large heights are exactly the admissible patterns. An
      admissible set in a window of length n has O(n/log n) elements (Hensley-Richards scale). So
      L(n) = 2^{o(n)} and the HL value of h_info is 0.
- ML:
  - A. Kolpakov and A. Rocke, "On the impossibility of discovering a formula for primes using AI",
    arXiv:2308.10817 (2023), and "Machine learning of the prime distribution", arXiv:2403.12588,
    PLOS ONE 19(9) e0301240 (2024). They use maximum-entropy arguments (Erdős-Kac, Hardy-Ramanujan),
    not HL, to explain Yang-Hui He's experiments on the learnability of primes.
  - F. Charton, "Learning the greatest common divisor", ICLR 2024 (arXiv:2308.15594): transformers
    learn divisibility by the base's divisors and by small primes.
  - Wu et al. 2023, "On the prime number divisibility by deep learning".
  - None quantifies the extractable information by HL.

Searches with no hit: "entropy of prime indicator sequence Hardy-Littlewood mutual information";
"transformer learns primality only small prime divisibility wheel".

## 5. Bias between twin pairs (12j-1, 12j+1) and (12j+5, 12j+7)

**Verdict: PARTLY KNOWN.** The race has been tabulated numerically (in the equivalent mod-4 form)
with no explanation. I found no conjecture or theory.

- For twins (p, p+2) with p > 3, p = 5 mod 6. So p = 1 mod 4 ⇔ p = 5 mod 12 ⇔ the centre is 6 mod 12,
  and p = 3 mod 4 ⇔ p = 11 mod 12 ⇔ the centre is 0 mod 12. The mod-4 split of the first member is
  exactly the requested race.
- S. Sahoo, "On twin prime distribution and associated biases", arXiv:2111.09053 (v3 2023), §9.1
  "Type-I bias", Table 1 (read via ar5iv).
  - Δ2 = π2(x;4,3) − π2(x;4,1) = +16 at 10^7 and −505 at 10^10.
  - The lead changes near 5·10^4, near 10^7 and near 1.1·10^10.
  - The paper notes twins lean toward a = 1, unlike primes. No mechanism is given.
  - §9.2 adds an Oliver-Soundararajan-type bias for consecutive twin pairs, with Conjecture 9.1
    that it tends to 1/2.
- W. Puszkarz, arXiv:1807.00406 (2018): twin centres favour non-squarefree multiples of 6 (2.427
  against a naive 2.290).
  - My remark, not in the paper: this looks like the HL local factor. For q >= 5, under HL a twin
    centre n is spread evenly over the q − 2 residues n ≢ ±1 mod q. So q | n with probability
    1/(q−2), not 1/q. It says nothing about 0 versus 6 mod 12, which HL treats equally.

Searches with no hit: "twin primes Chebyshev bias residue classes mod 12 prime races"; "twin prime
race mod 4 / mod 12 bias explanation second order term Rubinstein-Sarnak".

## 6. Number of primes q <= x with ord_q(2) prime

**Verdict: NOT FOUND** (no asymptotic or heuristic formula). The sequence itself is catalogued.

- OEIS A122094 (local copy): "Prime divisors of Mersenne numbers. Primes p such that the
  multiplicative order of 2 modulo p is prime" (Alekseyev 2006, b-file 10,000 terms). It has no
  counting function or asymptotic. Related: A137332 (the prime orders) and A089162 (factors of
  2^prime − 1). No OEIS name matches "number of ... up to 10^n" for this set.
- Closest bounds: C. Pomerance, "Cyclotomic primes", talk, 28 May 2025 (slides read).
  - #{p <= x : ℓ(p) < x^{1/2}/log x} < x/(log x)^2.
  - Under prime k-tuples, #{p <= x : p = 3 mod 4, 2p+1 prime} ~ c·x/(log x)^2. Each such q = 2p+1
    has ord_q(2) = p, so this gives a conditional lower bound of order x/(log x)^2 for our count.
- Moree's papers on the order of g mod p in residue classes (J. Number Theory 2005; arXiv
  math/0211259, math/0212220, math/0608468) and Kurlberg-Pomerance on Arnold's problem treat
  divisibility conditions on the order, not primality of the order.
- Note (my own heuristic, not literature): Shanks-Kravitz/Wagstaff with Bateman-Horn give a count
  of order x/(log x)^2, dominated by small k in q = 2kℓ+1. I found no paper that writes the
  constant.

Searches with no hit: "Moree order of 2 prime order density count"; "number of primes p <= x such
that multiplicative order of 2 is prime asymptotic"; "ord_p(2) is prime counting function
Sophie Germain analogue"; "multiplicative order of 2 is prime x/(log x)^2 Kurlberg Pomerance
Shparlinski".

## 7. p = 1 mod 4, q = 6p+1 prime, q = x^2 + 27y^2 ⇒ q | 2^p − 1 (and 8p+1)

**Verdict: KNOWN.** It is even known as an "if and only if".

- Zhi-Wei Sun, message to the NMBRTHRY list, 1 March 2008, proving a conjecture of O. Latinne
  (29 Feb 2008, tested to p = 10^11). Text read as reposted in the Physics Forums thread "Is there a
  proof for the new Mersenne numbers conjecture?" (physicsforums.com/threads/...218948).
  - Statement: 6p+1 | 2^p − 1 iff p = 1 mod 4 and 6p+1 is a prime x^2 + 27y^2.
  - Proof of "if": 2 is a cubic residue (Ireland-Rosen Cor. 9.6.2), so q | 2^{2p} − 1 =
    (2^p − 1)(2^p + 1); then 2^p = −1 is ruled out by the quadratic character.
  - Proof of "only if": a composite d has a factor < 2p+1.
- Same message, the 8p+1 analogue: for prime p, 8p+1 | 2^p − 1 iff 8p+1 is a prime x^2 + 64y^2
  with y odd. It cites Berndt-Evans-Williams, *Gauss and Jacobi Sums* (1998), Cor. 7.5.8: 2 must be
  an octic residue.
  - Consistency check (mine): q = 8p+1 = 9 mod 16, so (−1/q)_8 = (−1)^p = −1. With BEW's
    (−2/q)_8 = (−1)^y this gives (2/q)_8 = 1 iff y is odd.
- Earlier, in probabilistic form: Shanks-Kravitz 1967 (item 2), k = 3 case. "any p = 6n+1 can be
  written" a^2 + 3b^2, "but only those p where 3 | b have 2 as a cubic residue", so 1/3 of such
  6p+1 divide M_p.
- OEIS (local copy):
  - A188130 (Hasler 2011): primes p with 6p+1 | M_p; "p = 1 (mod 4) and 6p+1 is prime, but there
    are other primes with these properties (13, 17, ...) not in this sequence".
  - A122095: 8p+1 | 2^p − 1.
  - A188133: 10p+1 | 2^p − 1.
  - A350702: 14p+1 | 2^p − 1.
  - A014752: primes x^2 + 27y^2 = primes 1 mod 3 with 2 a cubic residue; it cites Ireland-Rosen
    Prop. 9.6.2, p. 119.
- Possibly earlier, not verified: N. Y. Wilson, "Conjectures as to a factor of 2^p ± 1", Edinburgh
  Math. Notes 39 (1954) 6-9. It states seven conjectured "theorems" on factors of 2^p ± 1 and is
  cited by Alderson 1964 (quintic character of 2) and by Lemmermeyer's *Reciprocity Laws* (2000).
  I could not read the body.

## 8. Chebyshev bias in Goldbach representations by residue class of the summands

**Verdict: PARTLY KNOWN.** The analytic framework exists: averages of Goldbach counts with summands
in arithmetic progressions, with secondary terms from zeros of Dirichlet L-functions. I found no
paper that states a bias between class pairs such as (1, 11) and (5, 7) mod 12, or a prime-square
correction to Goldbach counts.

- F. Rüppel, "Convolutions of the von Mangoldt function over residue classes", Šiauliai Math.
  Semin. 7(15) (2012) 135-156. Studies Σ Λ(m)Λ(n−m) with m, n−m in residue classes, as cited by
  the next two.
- Y. Suzuki, arXiv:1504.01967 (published 2017): mean value of the representation function for
  sums of two primes in progressions, improving Rüppel.
- G. Bhowmik, K. Halupczok, K. Matsumoto and Y. Suzuki, "Goldbach representations in arithmetic
  progressions and zeros of Dirichlet L-functions", Mathematika 65 (2019) 57-97 (arXiv:1704.06103;
  abstract read). It gives asymptotics for the average number of such representations under a
  distinct-zeros conjecture.
- Further work on the same averages: Cantarini-Gambini-Zaccagnini arXiv:1911.07220 (Cesàro
  averages with summands in progressions) and T. T. Nguyen arXiv:2407.18266 (2024).
- I did not open the full texts, so I cannot say whether any isolates the −x^{1/2}-type term from
  prime squares as a "bias".
- Not peer-reviewed:
  - A DEV Community post (2026, "subhansh"): an n mod 3 bias of 0.26% in Goldbach counts, no
    citations.
  - C. Spackman's 2018 blog: only the ordinary mod-4 prime race.

Searches with no hit: "Chebyshev bias Goldbach representations residue classes prime squares";
"Goldbach prime number race Chebyshev's bias additive problem"; "Goldbach representations residue
classes bias p = 1 mod 4 vs 3 mod 4".

## 9. Residues mod 3 (or mod 12/24) of the prime factors of Mersenne numbers

**Verdict: NOT FOUND** (as an aggregated statistic or a bias claim).

- For p > 3, q = 2kp+1 = 1 mod 3 exactly when 3 | k. So the mod-3 split of factors is a statement
  about k mod 3.
- Shanks-Kravitz 1967 tabulate per-k counts (Table 1, k < 60, x <= 10^6; an image I could not
  read). Their conjecture (3) has the factor over odd q | k, so the 3 | k classes get their own
  density. No aggregated mod-3 or mod-24 statistic is drawn there or anywhere else I found.
- Wagstaff 1983 averages C2·f(2k) to 1 (item 1). Those are exactly the factors that would create
  the mod-3 effect.

Searches with no hit: "distribution of prime factors of Mersenne numbers modulo 3 / mod 8 bias";
mersenneforum "factors mod 3 / 1 mod 8 / 7 mod 8 distribution k values".

## 10. One dominant low zero in prime races; the zero of L(s, χ_{-19}) at 1.5161

**Verdict: KNOWN** (the substance; the phrase "one dominant zero" is not used).

- C. Bays, K. Ford, R. H. Hudson and M. Rubinstein, "Zeros of Dirichlet L-functions near the real
  axis and Chebyshev's bias", J. Number Theory 87 (2001) 54-76, doi:10.1006/jnth.2000.2601 (preprint
  PDF read).
  - §1, point 2: "when q has a primitive root the bias depends heavily on the location of the first
    few zeros (closest to the real axis) of L(s, χ_q) ... the bias is heavily influenced by the size
    of the first zero". By the Chowla-Selberg formula, class number 1 for Q(√−q) implies a
    relatively low first zero.
  - Table I: the non-residue vs residue race mod 19 has bias ≈ 0.8042, against ≈ 0.99 for
    q = 3, 4, 5.
  - Table II: q = 19 with T = 10000, 1000 and 100.
  - §2: "for q = 19 we obtained with just 10 zeros (T = 20) a bias value which differs by less than
    1% from the bias value using the 3184 zeros to height T = 10000". Fig. 3 plots P*_{19;N,R}(x;T).
  - Table III: actual zeros of L(s, χ_d) against the Chowla-Selberg approximation for h(d) = 1. For
    d = −19 the first zero is 1.51608375 (approximation 1.51607047); the next are 5.4766 and 7.1607.
  - The smallest is for d = −163, at 0.2029.
- LMFDB, L-function 1-19-19.18-r1-0-0 (read): first zero 1.516083753160053883697692705674.
- D. Fiorilli, "Highly biased prime number races", Algebra & Number Theory 8 (2014) no. 7
  (arXiv:1210.6946): "low-lying zeros (excluding real zeros) have a significant effect on
  decreasing the bias", citing BFHR and Fiorilli-Martin.
- D. Fiorilli and G. Martin, "Inequities in the Shanks-Rényi prime number race", J. reine angew.
  Math. 676 (2013) (arXiv:0912.4908): asymptotic formula for the densities. I found no mod-19
  specifics there.

---

## Summary table

| # | Verdict | Best reference |
|---|---------|----------------|
| 1 | PARTLY KNOWN (mod 4 only) | Wagstaff, Math. Comp. 40 (1983) 385-397 (a = 2/6; C2 f(2k) replaced by 1 on p. 386, per Wagstaff-Montgomery-Nahm, Math. Comp. 79 (2010) p. 1799); Ehrman, Math. Comp. 21 (1967) 700-704 |
| 2 | PARTLY KNOWN (per-k) / NOT FOUND (per-p) | Shanks-Kravitz, Math. Comp. 21 (1967) 97-101 (1/k law, conjecture (3)) |
| 3 | PARTLY KNOWN (method) / NOT FOUND (byte table) | Brent, Math. Comp. 28 (1974) 315-324; Lemke Oliver-Soundararajan, PNAS 113 (2016) |
| 4 | PARTLY KNOWN | Lacasa et al., Entropy 20 (2018) 131; Blank arXiv:2403.04490 (Thm 5.1(b) looks wrong); Kolpakov-Rocke PLOS ONE 2024 |
| 5 | PARTLY KNOWN (data, no theory) | Sahoo arXiv:2111.09053 §9.1, Table 1 (mod-4 split of the first member = the mod-12 race) |
| 6 | NOT FOUND (formula) | OEIS A122094 (sequence only); Pomerance "Cyclotomic primes" (2025) slides 17, 22 (bounds) |
| 7 | KNOWN (iff) | Z.-W. Sun, NMBRTHRY 1 Mar 2008 (Latinne's conjecture); Shanks-Kravitz 1967 (k = 3); OEIS A188130, A122095 |
| 8 | PARTLY KNOWN (framework only) | Bhowmik-Halupczok-Matsumoto-Suzuki, Mathematika 65 (2019) 57-97; Rüppel 2012 |
| 9 | NOT FOUND | (implicit in Shanks-Kravitz per-k densities) |
| 10 | KNOWN | Bays-Ford-Hudson-Rubinstein, J. Number Theory 87 (2001) 54-76, Tables I-III, §2 |
