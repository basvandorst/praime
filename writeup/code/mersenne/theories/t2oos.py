"""Theory 2 out of sample: the local factor of no. 52, 136,279,841 (above the verified bound, used by no fit), against the primes
around it in the same class mod 4; and the lane shares that theory 1 forecasts for no. 53 (median forecast 184 million)."""
import numpy as np, math, itertools
SMALL = [3, 5, 7, 11, 13]
K = 2000
k = np.arange(1, K + 1, dtype=np.int64)
corr = np.ones(K)
for i in range(K):
    n = int(k[i]); r = 17
    while r <= n:
        if n % r == 0 and all(r % d for d in range(2, int(r ** 0.5) + 1)): corr[i] *= (r - 1) / (r - 2)
        r += 1
C = 0.6601618158468696
for r in SMALL: C /= r * (r - 2) / (r - 1) ** 2
def logF(p, res):
    allowed = np.isin((2 * ((k * res[4]) % 4) + 1) % 8, (1, 7))
    f = 2.0 * C * corr * allowed / np.log(2.0 * k * p + 1)
    for r in SMALL:
        f = f * np.where((2 * k * res[r] + 1) % r == 0, 0.0, r / (r - 1))
    return float(np.sum(np.log1p(-f / k)))
p52 = 136279841
lo, hi = 135_000_000, 137_500_000
seg = np.ones(hi - lo, bool)
for q in range(2, int(hi ** 0.5) + 1):
    if all(q % d for d in range(2, int(q ** 0.5) + 1)):
        st = max(q * q, ((lo + q - 1) // q) * q); seg[st - lo::q] = False
ps = lo + np.nonzero(seg)[0]
ps = ps[ps % 4 == p52 % 4]
keys = {}
vals = []
for p in ps:
    key = (int(p % 4),) + tuple(int(p % r) for r in SMALL)
    if key not in keys:
        keys[key] = logF(1.36e8, {4: key[0], **{r: key[j + 1] for j, r in enumerate(SMALL)}})
    vals.append(keys[key])
vals = np.exp(np.array(vals)); vals /= vals.mean()
key52 = (p52 % 4,) + tuple(p52 % r for r in SMALL)
w52 = math.exp(keys[key52]) / np.mean(np.exp([keys[(int(p % 4),) + tuple(int(p % r) for r in SMALL)] for p in ps]))
print(f'no. 52 = {p52:,}: residues mod 4, 3, 5, 7, 11, 13 = {key52}; local factor {w52:.3f} (1 = average of the {len(ps):,} primes = {p52 % 4} mod 4 in [{lo:,}, {hi:,})); share of those primes with a smaller factor: {np.mean(vals < w52):.3f}')
