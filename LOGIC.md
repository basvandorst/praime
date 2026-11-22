# The logic of the primes: one identity, three layers

> **In short:** Every regularity this repo has measured in the primes belongs to one of three layers of a single exact identity, the explicit formula: a local layer (each prime p on its own: residues, the sieve, the singular series), an archimedean layer (size: the density 1/ln x, the pole, the gamma factor) and a spectral layer (the zeros). The layers are not independent: the small primes fix the lowest zeros, the lowest zeros set the races, and the missing prime below 2 makes the zeros rigid. Tested end to end on 3,040 L-functions: rank correlations 0.87, 0.99 and 0.98. The anomalies left in the repo sort into lower-order terms of one layer, or a noise model that ignores the rigidity.

The documents this map draws on: [STUDY.md](STUDY.md) (seven lines, rounds 1–105), [RIEMANN.md](RIEMANN.md) (sections 1–23), [MERSENNE.md](MERSENNE.md) (sections 2–55), [FACTS.md](FACTS.md), [FINDINGS.txt](FINDINGS.txt). Each statement carries its status: *theorem*, *computed*, *heuristic* or *open*; *ours* marks what this repo derived or measured first, as far as we found.

## 1. The identity

For a test function h with Fourier transform ĥ, Weil's explicit formula (*theorem*) reads, schematically,

  Σ_ρ h(γ) = [pole terms, only for ζ] + (1/2π)∫ h(r)·[ln(q/π) + Re Γ′/Γ-term] dr − 2 Σ_n Λ(n) Re χ(n) ĥ(ln n)/√n,

with the zeros ρ = ½ + iγ on the left, the archimedean place (the pole of ζ at s = 1 and the gamma factor of the completed L-function) in the middle, and every prime power on the right. Riemann's formula ψ(x) = x − Σ x^ρ/ρ − … is the case of a sharp cut-off; the variance laws, the tube and the races of this repo are other choices of h. The identity is exact, so no statistic of the primes can be "free": each one is a local sum, plus a size term, plus a sum over zeros.

## 2. The three layers

| layer | what it is | what it produces in this repo |
|---|---|---|
| **local** | each prime p separately: which residues mod p are allowed | the sieve and wheels (MERSENNE §17–20); the singular series S(H) and the Hardy–Littlewood k-tuple densities (STUDY formula box); gaps mod 6, the exact gap law, the conditioned law (STUDY rounds 22–55); the last-digit bias and the k-th neighbour law (rounds 79–90; substance Lemke Oliver–Soundararajan); Ulam's diagonals and Bateman–Horn (rounds 9, 18, 103); the record prime-producing polynomials (first prime that can divide a value: 41, 59, 47, 67, 59, 79, 89 for degrees 2–6, STUDY polynomial section); the prime squares behind Chebyshev's bias (facts 1602, 1605; MERSENNE §9–10) |
| **archimedean** | the size of n: the place p = ∞ | the density 1/ln x (the pole of ζ at 1); Cramér's ln X − 3/2; the 7/8 in N(T) (the gamma factor), which gives the (7/4)δ in the exact variance D(δ) (RIEMANN §23, round 95); the pole term 2 ln(1 + δ) that the race between classes lacks (round 102); the tube constant C = 2 + γ − ln 4π (fact 2048) |
| **spectral** | the zeros | ψ(x) − x and every fluctuation: the variance R(X, δ) = ln(1/δ) + B (Goldston–Montgomery, rounds 91–93), its swing and the zeros' rigidity (rounds 96–100), the races (round 105), the tube and its fill fraction 0.53–0.88 (§22), decoding primes from zeros (RIEMANN §1–16) |

## 3. The links, tested

Each link joins two layers. All are measured in this repo; the theory behind each is classical unless marked ours.

**L1. Local → spectral, ζ: the pair constants are the zeros' pair correlation.** The Hardy–Littlewood pair constants summed over a window give Montgomery–Soundararajan's variance; the same variance read from the zeros (Goldston–Montgomery) is ln(1/δ) + B. *Computed*: primes alone give R within 0.065 of ln(1/δ) + B for δ = 10⁻³ … 10⁻⁶ at 10⁹ (rounds 91–92, pre-registered); the zeros up to 2,949 reproduce the primes' R(X, δ) within 0.003 at X = 10⁸, 10⁹, 10¹⁰ and δ = 0.1, 0.01 (round 93, a numerical identity).

**L2. Archimedean → variance.** With a test function whose Fourier support lies below ln 2, no prime enters the explicit formula, so the mean variance is fixed by the pole and the gamma factor alone: D(δ) = ln(1/δ) − γ − ln 2π + 2 ln(1 + δ) + ((2 + δ)/δ) ln(1 + δ/2) (*ours*, via Guinand–Weil, under RH and simple zeros; round 95, numerical check 2·10⁻⁸). The race between residue classes has no pole: D_χ(δ) lacks 2 ln(1 + δ) and gains ln q (*ours*; round 102, three of four readings pass). So the "+1" of the counting function, the pole, shows in the primes as an extra variance that races do not have.

**L3. Local → the lowest zero.** By Hadamard's product, Σ_ρ 1/|ρ|² = 2 Re L′/L(1, χ) + ln(q/π) + ψ((1 + a)/2) (*theorem*; fact 2049). L′/L(1, χ) weighs every prime by about ln p/p, so characters that are −1 on many small primes have a large right side and must have a zero near the central point. *Computed* (round 104, 3,040 fundamental discriminants with |D| ≤ 5,000): the predicted first zero, from L′/L(1, χ) alone with the other zeros taken as typical, has rank correlation 0.987 (D < 0) and 0.986 (D > 0) with the true one. For complex characters, the family of RIEMANN §18.1 (where small primes alone fixed the lowest zero to R² = 0.45 mod 395,137), round 107 tested all 1,007 characters mod 1009: rank correlation 0.937 between predicted and true lowest zero (passes), but its registered straight-line measure failed (R² 0.22 against 0.745 for a fit on Re χ(p), p ≤ 199); in the identity's own variable, 1/(¼ + γ²) against S, Pearson is 0.926 (chosen after the fact). So the weighted sum over all primes is the natural variable for §18.1's open point, with that caveat.

**L4. Local → prime-rich polynomials → low zeros.** A quadratic n² + bn + c is prime-rich when no small prime divides its values, i.e. χ_D(p) = −1 for small p; its Bateman–Horn constant is ∏(1 − χ(p)/(p − 1)), essentially 1/L(1, χ_D) (*theorem* for the identity of products; the density is *conjecture*). With L3 the prime-rich quadratics are the ones with low zeros: Spearman 0.87 between ln L(1, χ_D) and the scaled first zero, for both signs (round 104, pre-registered). Euler's n² + n + 41 sits on L(s, χ₋₁₆₃), whose first zero is at 0.2029 (11th lowest of 1,524 negative discriminants up to 5,000; χ₋₄: 6.02).

**L5. The lowest zero → the race.** The race between the −1 and the +1 primes of χ_D is E(x) = 1 + 2 Σ Re(x^(iγ)/ρ) (Rubinstein–Sarnak; the 1 is the prime squares, a local term). A low zero is a slow wave in ln x of period 2π/γ₁: with γ₁ ≈ 0.07 the −1 side leads at every prime below 10⁵, with γ₁ ≈ 0.3 the wave has turned near 10³ and the +1 side leads. *Computed* (round 105): the zeros up to height 30 predict the share of primes at which the −1 side leads with Pearson 0.981 (mean error 0.039); the lowest zero alone gives 0.889. Our pre-registered guess, that a lower zero means a larger bias, failed with the wrong sign: the zero sets the phase, not the size.

**L6. The missing prime below 2 → rigid zeros.** In the zeros' form factor the primes enter at τ = ln p^k/ln(T/2π); since no prime power is below 2, the form factor has no weight below ln 2 (Berry, Bogomolny–Keating; *heuristic*). *Computed*: the swing of R(X, δ) over ln X is below stationary GUE by 3.2–4.6σ at δ = 0.01–0.003, and the form-factor prediction accounts for it (round 100; its overall normalisation was not re-derived independently). So the zeros are more rigid than random matrices because of a local fact at the smallest prime.

**L7. Local → the records.** Every record prime-producing polynomial of degree 2 to 6 has a late first prime P that can divide a value (41 to 89); values below P² are then automatically prime (*computed*, STUDY polynomial section). The cubic record is a cubic analogue of Euler's: inert at 5, 7, 13, …, 43, with 2, 3, 11 sent to infinity by its leading coefficient 66. By L3–L4 the quadratic records sit on low zeros; for degrees ≥ 3 the zeros of ζ_K/ζ were beyond our machine (PARI ran out of memory at discriminant 8·10¹²), so that link is *open*.

The chain L3 → L4 → L5 is measured end to end on the same 3,040 L-functions: the small primes set L′/L(1, χ); that one number sets how prime-rich the quadratic is (0.87), where the lowest zero lies (0.99) and the phase of the race (0.89, 0.98 with more zeros). Three topics that the repo treated apart, Euler's polynomial, low zeros and Chebyshev's bias, are one number read three ways.

## 4. The anomalies, sorted

Two catalogues of everything in the repo that did not fit (one agent per side, both read-only) give the list below; each item names the layer it belongs to. *Explained* means a measured explanation in the repo; *proposed* means a reading that has not been tested.

**Lower-order terms of the local layer (proposed: truncation and normalisation).**
- The conditioned gap law's residuals are all positive and grow with g (0.065 % at g = 700; STUDY round 53); the law's share of gaps ≥ 40 is 0.26 % high at 10⁸–10⁹ (round 41) and the two-gap law 0.27 % high in long-gap cells at 10⁹ but not at 10¹² (stress test); the twin-gap tail beyond 1,200 is high by 0.24 % (rounds 59, 64). The two-gap law's row sums are 1.0006–1.0012 and the k-th neighbour law's mass is off by 2·10⁻⁴: one order-2 normalisation error would produce all of these. Order 3 would decide.
- Runs of prime values of n² + n + 41: observed/expected falls from 1.000 (r = 1) to 0.856 (r = 8) (round 103); long patterns at small x overweight the small-t end of the Bateman–Horn integral.

**The archimedean layer at finite height (explained).** The (7/4)δ of the variance (round 95); the 7/8 of N(T); the race without the pole (round 102).

**The spectral layer: rigidity, and a noise model that ignores it (proposed: one cause).** The repo has a long list of χ² values below their degrees of freedom: last digits 20.0–25.4 on 36 (rounds 85–86), k-th neighbour distances 8.2 on 19 (round 83), block errors 0.75–0.87 of binomial (round 85), the residue tables mod 210 at 0.691 against 0.847 (rounds 13, 39), prime bytes, the lighthouse, pair clocks (MERSENNE §10, §45, §53), the Goldbach comet 13 % narrower than sieved fakes (§28). Their sizes, 0.6–0.85 of the binomial variance, match the measured rigidity of gap-event counts across residue classes, 0.60–0.80 of Poisson (round 40). Counts per class are governed by the zeros of the L-functions mod q, which are rigid (L6), so binomial or multinomial nulls overstate the noise everywhere. The core of it, twin counts spread across residue classes more evenly than the Hardy–Littlewood four-point sums predict (u² = 0.13 against 0.45 at 10⁸–10⁹; Fisher p ≈ 0.002 over three windows; rounds 41–44), stays *open*: twins have no L-function, so there is no explicit formula to supply the missing long-range term.

**Variance at long range (explained, round 106).** At N = 10⁷ the variance of ψ(n + h) − ψ(n) − h over n in [N/2, N] is 0.963, 0.924, 0.879 of Montgomery–Soundararajan at h = 10³, 3·10³, 10⁴ (MERSENNE §34). The explicit formula with 10,000 zeros gives 0.969, 0.936, 0.881: the window, only ln 2 long in ln x, sits in a trough of the zeros' slow swing (rounds 96–98). The values below the law at H = 10⁵–10⁶ in STUDY rounds 10 and 34 lie within their noise (t₉ from −0.35 to −1.6) and were not recomputed.

**Polynomial tails (open).** The heuristic chance per candidate overestimates the run-length tail of sextic searches by about 10³ (STUDY polynomial section).

**Explained by this map.** RIEMANN §18.1 (the lowest zero is fixed by L′/L(1, χ), L3); the low-zero discriminants with the "wrong" race sign (L5); the long-range variance shortfall of MERSENNE §34 (round 106).

## 5. What "not random" means here

The primes are not random in any layer. The local layer is a sieve and is fully described by the singular series; the archimedean layer is a smooth function of size. What looks random is the spectral layer, the sum over zeros, and that is not free either: the zeros are fixed by the primes (backwards decoding, RIEMANN §6 and §13; the lowest zero by one weighted prime sum, L3), and the zeros are more rigid than random matrices because of the primes (L6). The only ingredient that behaves like chance is that the heights γ are, as far as anyone knows, linearly independent over the rationals, so sums of x^(iγ) behave like sums of independent phases. That is the logic in one sentence: the primes and the zeros determine each other through one exact identity, and what remains random-like is only the incommensurability of the zeros' heights.

## 6. Open tests

- **L7 beyond degree 2:** the first zero of ζ_K/ζ for the cubic, quartic and quintic records against fields of similar discriminant (needs a larger machine or an explicit-formula estimate of the first zero).
- **Twins across classes:** a long-range correction to the four-point sums of the same type as L2, derived for pairs of primes.
- **The low-χ² family:** recompute the χ² of rounds 83–86 with errors from the rigidity of round 40 instead of binomial errors; the prediction is χ²/dof near 1.
- **Order 3 of the gap laws:** whether one normalisation term removes the 0.26–0.27 % long-gap excess and the twin-gap tail together.

## Validation

The numbers in section 3 come from [lowzero/r104.py](writeup/code/study/lowzero/r104.py), [r104b.py](writeup/code/study/lowzero/r104b.py), [r105.py](writeup/code/study/lowzero/r105.py) and [r105b.py](writeup/code/study/lowzero/r105b.py) (their JSON outputs next to them; PARI for zeros and L-values); the rest cite the rounds and sections where they were computed and checked.
