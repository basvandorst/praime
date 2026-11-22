"""Round 17 check: primes in [1e18, 1e18+2e6) by numpy small-prime sieve + sympy.isprime (deterministic < 2^64)."""
import json, math
import numpy as np, sympy
L, W = 10 ** 18, 2 * 10 ** 6
EXT = 3000
bp = np.array(list(sympy.primerange(2, 100000)))
M = np.ones(W + EXT, dtype=bool)
for p in bp:
    p = int(p); st = (-L) % p
    M[st::p] = False
cand = np.nonzero(M)[0]
P = np.array([int(c) for c in cand if sympy.isprime(L + int(c))], dtype=np.int64)
inwin = P[P < W]
P = P[:len(inwin) + 1]                  # close the last gap
print('primes in window', len(inwin), 'candidates tested', len(cand))
g = np.diff(P)
lam = math.log(L + W / 2)
TH = [int(round(k * lam / 2)) * 2 for k in (3, 4, 5)]
d = {}
for h in (2, 4, 6, 8, 10, 12, 18, 30):
    d['S.gap%d' % h] = float(np.mean(g == h))
d['N.eqgap'] = float(np.mean(g[1:] == g[:-1]))
for q in (3, 4, 5, 7, 8):
    r = (P + L % q) % q
    d['A.same%d' % q] = float(np.mean(r[1:] == r[:-1]))
for t in TH:
    d['T1.tail%d' % t] = float(np.mean(g >= t))
c = np.bincount(np.minimum(g // 2, 20), minlength=21)[1:].astype(float); c = c[c > 0]; p = c / c.sum()
d['T4.H_class'] = float(-np.sum(p * np.log2(p)) + (len(c) - 1) / (2 * c.sum() * math.log(2)))
R = json.load(open('/home/bas/praime/writeup/data/study/r17.json'))
for k, v in d.items():
    print('%-12s mine %.6g  repo %.6g  %s' % (k, v, R['stats'][k]['x'], 'OK' if abs(v - R['stats'][k]['x']) < 1e-12 else 'DIFF'))
# law check: exact HL sums for gaps <= 30 at lam vs repo order-4 law
E = json.load(open('/home/bas/scratch/study-review3/enum.json'))
for h in (2, 6, 10, 30):
    U = E['U'][str(h)]
    ex = sum((-1) ** j * U[j] / lam ** (j + 1) for j in range(len(U)))
    print('gap %d law: exact %.6f repo %.6f' % (h, ex, R['stats']['S.gap%d' % h]['T']))
# sub-range recount: [1e18, 1e18+2e5) by gmpy2-free plain Miller-Rabin bases (deterministic < 3.3e24)
def mr(n):
    if n % 2 == 0: return False
    dd, s = n - 1, 0
    while dd % 2 == 0: dd //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, dd, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True
sub = sum(1 for n in range(L, L + 200000) if mr(n))
print('sub-range [1e18,1e18+2e5): MR count', sub, ' sieve+sympy count', int(np.sum(P < 200000)))
