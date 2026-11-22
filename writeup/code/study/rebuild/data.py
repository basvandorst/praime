"""Own segmented sieve: primes in [1e9, 1.1e9 + 4000); repeat stats for p_n in [1e9, 1.1e9)."""
import json, math
import numpy as np

LO, HI = 10**9, 11 * 10**8
TOP = HI + 4000
def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]
small = primes_upto(int(math.isqrt(TOP)) + 1)[1:]   # odd primes
SEG = 10**7
chunks = []
for s0 in range(LO, TOP, SEG):          # s0 even; index j <-> s0 + 2j + 1
    s1 = min(s0 + SEG, TOP)
    n = (s1 - s0) // 2
    isp = np.ones(n, bool)
    for q in small:
        f = max(q * q, ((s0 + 1 + q - 1) // q) * q)
        if f % 2 == 0: f += q
        isp[(f - s0 - 1) // 2::q] = False
    chunks.append(s0 + 1 + 2 * np.nonzero(isp)[0].astype(np.int64))
pr = np.concatenate(chunks)
N = int(np.sum(pr < HI))
print('primes in window', N, 'first', pr[0], 'last<HI', pr[N-1])
res = {'N': N}
hist = {}
for k in range(1, 5):
    D = pr[k:N + k] - pr[:N]
    hist[k] = np.bincount(D, minlength=241)[:241].tolist()
    res[k] = {}
    for Dm in (200, 240):
        sel = D[D <= Dm]
        res[k][Dm] = {'n': int(len(sel)), 'frac_all': len(sel) / N,
                      'rep': {q: int(np.sum(sel % q == 0)) for q in (3, 5, 7, 11, 13)}}
res['hist'] = hist
json.dump(res, open('/home/bas/scratch/rebuild-C/data.json', 'w'))
for k in range(1, 5):
    r = res[k][240]; print(k, r['n'], {q: round(c / r['n'], 5) for q, c in r['rep'].items()})
