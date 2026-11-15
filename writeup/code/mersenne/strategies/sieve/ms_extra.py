"""Extra mersieve numbers (numpy; reads ms_first.txt and ms_analyse.out).  Rechecked by ms_extra_check.py.
1. Neighbours (previous and next prime) of the 21 exponents in [100, N]: survivors at K = 64.
2. Totals of the 1/k table in ms_analyse.out: all k, and 3 <= k <= 64.
3. Poisson(1) tail P(X >= 6) times the number of Gallagher windows (31,848)."""
import numpy as np, re, math
P, F = np.loadtxt("ms_first.txt", dtype=np.int64, unpack=True)
EX = [107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091]
pos = np.searchsorted(P, EX)
nb = np.concatenate([pos - 1, pos + 1])
assert len(set(nb.tolist())) == 42 and not set(P[nb].tolist()) & set(EX)
surv = (F[nb] == 0) | (F[nb] > 64)
rate = ((F == 0) | (F > 64))
# expected: each neighbour's band survivor rate (bands as in Q4)
bands = [(2, 10**4), (10**4, 10**5), (10**5, 500001)]
exp = 0.0
for i in nb:
    for lo, hi in bands:
        if lo <= P[i] < hi:
            b = (P >= lo) & (P < hi); exp += rate[b].mean()
print(f"neighbours of exponents >= 100: {int(surv.sum())}/42 survivors at K=64; expected {42*rate.mean():.1f} at the overall rate, {exp:.1f} at band rates")
print("  non-survivor neighbours (prime, row):", [(int(P[i]), int(F[i])) for i in nb if not ((F[i] == 0) | (F[i] > 64))])
rows = [l for l in open("ms_analyse.out") if l.strip().startswith("k ") and "admissible" in l]
tot = [0, 0.0, 0.0]; sub = [0, 0.0, 0.0]
for l in rows:
    a, b = map(int, re.search(r"k (\d+)-(\d+)", l).groups())
    h = int(re.search(r"divides (\d+)", l).group(1)); e = float(re.search(r"= ([\d.]+) \(sd", l).group(1)); sd = float(re.search(r"sd ([\d.]+)", l).group(1))
    for t in ([tot] + ([sub] if 3 <= a and b <= 64 else [])):
        t[0] += h; t[1] += e; t[2] += sd * sd
for name, t in (("all k <= 1000", tot), ("3 <= k <= 64", sub)):
    print(f"1/k table, {name}: found {t[0]}, expected {t[1]:.1f}, sd {math.sqrt(t[2]):.1f}, z = {(t[0]-t[1])/math.sqrt(t[2]):+.2f}")
print("Poisson(1) expected windows with >= 6 primes:", round(31848 * (1 - sum(math.exp(-1) / math.factorial(j) for j in range(6))), 1))
