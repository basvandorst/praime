"""Round 36 (diagnostic): Hardy-Littlewood (pairs only) prediction for the spread of prime counts across the 48 classes
mod 210 over a whole window, against round 13's measurement (r13f.json). Per class a the density is (210/48)/lam and
the pair density at distance 210k is c S(0,210k)/lam^2 with the class conditioning, so (as in review5/r13_link.py)
v_a = 1 - n r_a + 2 sum_k (n-k) r2(k) / (n r_a), spread = (48 v_a - v(W)) / 47 with v(W) the whole-window value."""
import os, json, math
import numpy as np
import theory2 as T2
import common as C

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
C2x2 = 2 * 0.66016181584686957393
r13f = json.load(open(os.path.join(OUT, 'r13f.json')))
res = {}
for win, (L, W) in (('train', C.WINDOWS['train']), ('A', C.WINDOWS['A']), ('B', C.WINDOWS['B']), ('C', (10 ** 15, 2 * 10 ** 7))):
    n = W // 210
    k = np.arange(1, n, dtype=np.int64)
    # S(0, 210k) = 2 C2 prod_{p | 210k, p > 2} (p-1)/(p-2); for a fixed class a mod 210 the local factors at 2,3,5,7
    # are replaced by the conditioned ones: both points are units mod 210, factor prod_{p|210} 1/(1-1/p)^2 against
    # (1-nu/p)/(1-1/p)^2 with nu = 1 (same class): S_a(0, 210k) = S(0,210k) * prod_{p|210} 1/(1 - 1/p).
    S = np.full(len(k), C2x2)
    for p in (3, 5, 7):
        S *= (p - 1) / (p - 2)
    for p in C.small_primes(n + 1):
        p = int(p)
        if p <= 7:
            continue
        S[(p - 1)::p] *= (p - 1) / (p - 2)               # p | k  (index k-1)
    Sa = S * np.prod([1 / (1 - 1 / p) for p in (2, 3, 5, 7)])
    lam = math.log(L + W / 2)
    ra = (210 / 48) / lam
    va = 1 - n * ra + 2 * float(np.sum((n - k) * Sa / lam ** 2)) / (n * ra)
    vW = r13f[win]['primes']
    spread = (48 * va - 0.0) / 47
    res[win] = {'v_a_HL_pairs': va, 'predicted_spread_HL': spread, 'measured_spread': r13f[win]['primes']}
    print(win, {k2: round(v, 3) for k2, v in res[win].items()}, flush=True)
json.dump(res, open(os.path.join(OUT, 'r36.json'), 'w'), indent=1)
