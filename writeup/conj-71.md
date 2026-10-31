# Conjecture 71: The twin-prime product is 3.19618…, not π

> **In short:** Werner Sand conjectured that the product of 1/(1 − 1/t) over the smaller members t of all twin prime pairs equals π. The product is fixed by Brun's constant, the famous sum of reciprocals of twin primes, and with its known value the product comes out at **3.1961782**, about 1.7% above π. Sand's own figure, 3.02437 at t ≈ 2.4·10¹⁰, fits this value: the product climbs so slowly that it would only pass π near 10³³. So the conjecture is false, assuming the accepted value of Brun's constant (which rests on the Hardy–Littlewood heuristic).

**Conjecture:** [primepuzzles.net](https://www.primepuzzles.net/conjectures/conj_071.htm)

**Status:** explained (heuristic disproof: the limit is 3.1961782 given Brun's constant; a rigorous disproof is out of reach)

## The question

With t running over 3, 5, 11, 17, 29, 41, … (the smaller member of each twin pair), Sand conjectured

```
∏ 1/(1 − 1/t) = π
```

He observed the partial products rising monotonically from 3/2, reaching 3.024374573 at t = 23,887,065,569. The page asks for comments.

## What the page had vs ours

| | On the page | Ours |
|---|---|---|
| The limit | conjectured to be π | **3.1961782**, computed from Brun's constant |
| Partial product at 2.4·10¹⁰ | 3.02437 | consistent: the predicted gap to the limit is a factor e^(1.32/ln t) ≈ 1.057 |

## How we found it

Take logarithms: ln P = Σ −ln(1 − 1/t) = Σ 1/t + Σ (1/(2t²) + 1/(3t³) + …). The second sum converges fast, so everything hinges on Σ 1/t over the smaller twins. That is almost Brun's constant

```
B₂ = (1/3 + 1/5) + (1/5 + 1/7) + (1/11 + 1/13) + … = 1.902160583…
```

which adds both members of each pair. The two are linked by

```
Σ 1/t = (B₂ + Σ (1/t − 1/(t+2)))/2 = (B₂ + Σ 2/(t(t+2)))/2,
```

and the last sum also converges fast. So the product is a simple function of B₂.

## The math

With twin primes up to 10⁸ for the fast-converging sums (the tails beyond 10⁸ are below 10⁻¹⁰):

```
Σ 2/(t(t+2))            = D
Σ (−ln(1−1/t) − 1/t)    = E
ln P = (B₂ + D)/2 + E
B₂ = 1.902160583104 (Sebah–Gourdon)   →  P = 3.196178168
B₂ = 1.902160583209 (Nicely 2010)     →  P = 3.196178169
```

The value of B₂ comes from summing twin primes up to about 10¹⁶ and estimating the tail with the Hardy–Littlewood twin-prime formula, so it is not a proven constant. But its uncertainty is around 10⁻⁹, and π is 0.055 away.

**Why the partial products look like they approach π.** By Hardy–Littlewood, the twin primes beyond x contribute about Σ_{t > x} 1/t ≈ 2C₂/ln x, with 2C₂ = 1.3203. So

```
P(x) ≈ 3.19618 · e^(−1.3203/ln x)
x = 10⁸:       2.9751 (computed exactly: 2.975117)
x = 2.4·10¹⁰:  3.0244 (Sand: 3.024375)
x ≈ 10³³:      π
```

The product does pass π eventually, but only near 10³³, so every computation will show values below π for a long time.

**Can it be proved?** P = π would need B₂ = 1.8677 instead of 1.9022. Ruling that out needs a proven *lower* bound for B₂ (or for the product), and the only ones available are partial sums. They stay below the needed value until about 10³³. So an unconditional disproof is not possible today; the disproof is as solid as the Hardy–Littlewood estimate of B₂.

## Validation

This script computes D and E from the twin primes below 10⁷, combines them with B₂, and checks the partial product against the tail estimate. It runs in a few seconds.

```python
import math
import numpy as np

X = 10**7
s = np.ones(X + 3, dtype=bool); s[:2] = False
for i in range(2, int((X + 2)**0.5) + 1):
    if s[i]: s[i*i::i] = False
t = np.nonzero(s[:X + 1] & s[2:X + 3])[0].astype(np.float64)       # smaller twin primes t <= 10^7

D = float(np.sum(2 / (t * (t + 2))))
E = float(np.sum(-np.log1p(-1 / t) - 1 / t))
B2 = 1.902160583104                                                # Brun's constant (Sebah-Gourdon)
P = math.exp((B2 + D) / 2 + E)
assert abs(P - 3.1961782) < 1e-6 and P - math.pi > 0.05

partial = float(np.exp(np.sum(-np.log1p(-1 / t))))                 # rigorous partial product up to 10^7
C2 = 0.6601618158
assert abs(partial - P * math.exp(-2 * C2 / math.log(X))) < 0.003   # matches the Hardy-Littlewood tail
print("limit P = %.7f (pi = %.7f); partial product to 1e7 = %.5f" % (P, math.pi, partial))
```
