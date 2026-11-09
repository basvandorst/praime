# Conjecture 36: How many primes must be subtracted from an even number before the difference is prime

> **In short:** T(N) counts the odd primes p = 3, 5, 7, … tried until N − p is prime. The conjecture says T(N)/ln²N stays below a constant (largest found 0.7006) and asks for the constant and the next records. (1) Each trial succeeds with probability about 2C₂·S(N)/ln N, so the record index tends to 1/(2C₂) = 0.757. (2) The Goldbach verification to 4·10¹⁸ gives the next records: 121005022304007026 needs 1055 trials, and no index above 0.7006 occurs below 4·10¹⁸.

**Conjecture:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_036.htm)

**Status:** explained (Q1: the constant is 1/(2C₂) = 0.7574, heuristic; Q2 and Q3 answered from the Goldbach verification to 4·10¹⁸)

## The question

Subtract the successive odd primes 3, 5, 7, … from an even number N until the difference is prime. The number of trials divided by (ln N)² is the "quality index". A table of 13 even numbers with index above 0.62 is given, the largest being 0.7006 at N = 1047610575836828 (838 trials).

- **Q1.** An argument for the conjecture (the index stays bounded), and the constant it tends to.
- **Q2.** The next even number needing more than 982 trials.
- **Q3.** The next even number with an index above 0.7006.

## What the page had vs ours

| | On the page | Ours |
|---|---|---|
| Q1 | no argument; the table suggests a bound near 0.7 | the index of the record-holders tends to **1/(2C₂) = 0.7574** (heuristic, with the extreme-value argument); the 69 known records are at 82–94 % of it |
| Q2 | records to 4.3·10¹⁶ (2004) | **121005022304007026**, least prime 8443, 1055 trials (Oliveira e Silva's verification, 2012; OEIS A025018) |
| Q3 | — | none below 4·10¹⁸: the largest index among the eight further records is 0.6905 (N = 258549426916149682, 1110 trials) |

## How we found it

T(N) is the position of the least prime in a Goldbach partition of N, so the records are exactly the records of the least Goldbach prime, which the Goldbach verification project recorded up to 4·10¹⁸. For the constant: the trials are nearly independent coin flips with a success chance that depends on N only through ln N and the odd prime factors of N, so T(N) is a geometric variable, and the maximum of X/2 geometric variables is ln(X/2) times the mean. Everything then follows from the Hardy–Littlewood probability that N − p is prime.

## The math

**Step 1 (heuristic, Hardy–Littlewood). One trial.** For an odd prime p ≪ N, N − p is an odd number of size N. It is divisible by an odd prime q iff p ≡ N (mod q). If q ∤ N, p (a prime, so p ≢ 0 mod q) hits the forbidden class with probability 1/(q−1) instead of the 1/q of a random odd number, which multiplies the chance of primality by (1 − 1/(q−1))/(1 − 1/q) = q(q−2)/(q−1)²; if q | N, the class p ≡ N ≡ 0 is impossible, factor (1)/(1 − 1/q) = q/(q−1). With the factor 2 for oddness, P(N − p prime) ≈ (2/ln N) · ∏_{q ∤ N} q(q−2)/(q−1)² · ∏_{q | N, q>2} q/(q−1) = 2C₂·S(N)/ln N, where C₂ = ∏_{q>2} (1 − 1/(q−1)²) = 0.6602 and S(N) = ∏_{q | N, q>2} (q−1)/(q−2) ≥ 1.

**Step 2 (heuristic). T(N) is geometric.** The trials for consecutive primes p are nearly independent (p mod q runs over the nonzero classes), so P(T(N) > t) ≈ (1 − π_N)^t ≈ exp(−π_N t), π_N = 2C₂S(N)/ln N. Mean ln N/(2C₂S(N)).

**Step 3 (heuristic, extreme values). The records.** Among the X/2 even numbers up to X, the largest T is reached where exp(−π_N T) ≈ 2/X, i.e. T ≈ (ln(X/2) + γ)·ln N/(2C₂S(N)). Record-holders have S(N) close to 1 (no small odd prime factor: in the table S(N) ≤ 1.067, and 1.000 for five of the thirteen numbers on the page), so the record index T/ln²N is about (ln(X/2) + γ)/(2C₂ ln X), which is 0.755 already at X = 4·10¹⁸ and tends to 1/(2C₂) = 0.7574. This is the constant of Q1, and it answers the conjecture in the sense of Cramér's model: the index is bounded, by a limit that the records approach from below.

**Step 4 (computed). Check against the data.** The normalised index T·2C₂·S(N)/ln²N, which the model puts near 1 for records, is 0.82–0.94 for all 69 known records (0.934 for the page's best, 0.6905·1.013·1.32 = 0.924 for the largest beyond it). The approach from below is slow, as it is for maximal prime gaps, where the data sit at 80–90 % of Cramér's ln²p for the same reason: the extreme-value step treats the trials as fully independent, and the deficit is the usual finite-range correction of sieve models.

**Step 5 (literature). Q2 and Q3.** The least primes of minimal Goldbach partitions set records at 69 even numbers up to 4·10¹⁸ (OEIS A025018, A025019, from the verification by Oliveira e Silva, Herzog and Pardi). The page's table ends at the 61st record (31284177910528922, least prime 7753, 982 trials). The 62nd is 121005022304007026 (8443, 1055 trials, index 0.6819): the answer to Q2. The remaining seven have indices 0.6592, 0.6905, 0.6655, 0.6646, 0.6756, 0.6585, 0.6625, so no even number below 4·10¹⁸ has an index above 0.7006 (Q3). All 69 records are rechecked: N − p is prime, N − q is composite for every odd prime q < p, and T = π(p) − 1.

**What this does not settle.** The limit 0.7574 is heuristic (it assumes the Hardy–Littlewood probabilities and independence); no bound on T(N) is proved for any constant, as that would imply Goldbach's conjecture.

## Validation

[code/conj36_verify.py](code/conj36_verify.py) recomputes the trial counts of the page's table, rechecks the 69 records of [data/conj36_records.txt](data/conj36_records.txt), prints their indices and normalised indices, answers Q2 and Q3 from them, and prints the model's record index at several scales. Run after `pip install sympy`:

```
python3 conj36_verify.py      # a few seconds
```
