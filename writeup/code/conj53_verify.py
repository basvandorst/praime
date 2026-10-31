# Conjecture 53 (Capelle): checks every computed claim of writeup/conj-53.md (about 15 seconds).
# r < q < p consecutive primes, d = q - r, e = p - q. Identities 1-19 from the page and the new identity
#   floor(pq/r - qr/p) = 2(p - r) - [p - q > q - r]   (answer to Capelle's question)
import numpy as np
from fractions import Fraction as Fr
N = 200_000_001; sv = np.ones(N, bool); sv[:2] = False
for i in range(2, int(N**.5) + 1):
    if sv[i]: sv[i*i::i] = False
P = np.nonzero(sv)[0].astype(np.int64); del sv

# 1. the sufficient condition (p - r)^2 < r fails only for r <= 1327 below 2e8
r, p = P[:-2], P[2:]
assert int(P[np.nonzero((p - r)**2 >= r)[0]].max()) == 1327

# 2. every identity from its definition (exact rationals), all r < 10^6; failures as listed
def F(x): return x.numerator // x.denominator
fails = {}
Ps = [int(x) for x in P[P < 10**6 + 500]]
for i in range(len(Ps) - 2):
    r_, q_, p_ = Ps[i], Ps[i + 1], Ps[i + 2]
    if r_ >= 10**6: break
    a, b, c = Fr(p_ * q_, r_), Fr(q_ * r_, p_), Fr(r_ * p_, q_); fa, fb, fc = F(a), F(b), F(c)
    chk = {1: fa == p_ + q_ - r_, 2: fb == q_ + r_ - p_, 3: fc == r_ + p_ - q_ - 1, 4: F(a + b) == 2 * q_, 5: F(a + c) == 2 * p_,
           6: F(b + c) == 2 * r_, 7: F(a - c) == 2 * (q_ - r_), 8: F(b - c) == 2 * (q_ - p_), 9: F(a + b + c) == p_ + q_ + r_,
           10: fa + fb == 2 * q_, 11: fa + fc == 2 * p_ - 1, 12: fb + fc == 2 * r_ - 1, 13: fa - fb == 2 * (p_ - r_),
           14: fa - fc == 2 * (q_ - r_) + 1, 15: fb - fc == 2 * (q_ - p_) + 1, 16: fa + fb + fc == p_ + q_ + r_ - 1,
           17: fa * fb == q_**2 - (p_ - r_)**2, 18: fa * (fc + 1) == p_**2 - (q_ - r_)**2, 19: fb * (fc + 1) == r_**2 - (p_ - q_)**2,
           20: F(a - b) == 2 * (p_ - r_) - (1 if p_ - q_ > q_ - r_ else 0)}
    for k, v in chk.items():
        if not v: fails.setdefault(k, []).append(r_)
assert fails[1] == [2, 3, 5, 7, 11, 13, 19, 23, 31, 47, 83, 89, 113, 199, 1327]          # the page's list
assert fails[2] == [2, 3, 5, 11, 17, 19, 29, 43, 47, 83, 109, 199, 283] and fails[3] == [5, 19]
assert fails[20] == [2, 3, 5, 7, 13, 19, 23, 109, 113]
assert max(max(v) for v in fails.values()) == 1327                                       # nothing fails above 1327
assert all(set(fails[k]) <= set(fails[1]) | set(fails[2]) | set(fails[3]) or k in (4, 5, 6, 7, 8, 9, 20) for k in fails)

# 3. above 2e8 (to 4e18) every gap is <= 1476, so (p - r)^2 <= 2952^2 = 8714304 < 2e8 < r
assert 2952**2 == 8714304 < 2 * 10**8
print("all claims verified")
