"""Impostors: odd composite n with 2^(n-1) = 1 mod n (base-2 Fermat pseudoprimes, OEIS A001567) and Carmichael numbers
(A002997, Korselt's criterion: squarefree, at least 3 prime factors, p - 1 | n - 1 for every p | n), counted up to X.
Vectorised square-and-multiply on uint64 (products below 2^64 for n < 2^32).  Literature counts to compare with
(Pinch): pseudoprimes to base 2: 10^4: 22, 10^5: 78, 10^6: 245, 10^7: 750, 10^8: 2057, 10^9: 5597;
Carmichael numbers: 10^4: 7, 10^5: 16, 10^6: 43, 10^7: 105, 10^8: 255, 10^9: 646.   usage: impostors.py X out.txt"""
import sys, math
import numpy as np
X = int(sys.argv[1]); out = open(sys.argv[2], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
r = int(math.isqrt(X)) + 1
bs = np.ones(r + 1, bool); bs[:2] = False
for i in range(2, int(r ** .5) + 1):
    if bs[i]: bs[i * i::i] = False
base = np.nonzero(bs)[0]
psp = []
seg = 2 * 10 ** 6
for lo in range(3, X + 1, seg):
    hi = min(X + 1, lo + seg); comp = np.zeros(hi - lo, bool)
    for p in base:
        st = max(p * p, -(-lo // p) * p)
        if st < hi: comp[st - lo::p] = True
    n = np.arange(lo, hi, dtype=np.uint64); m = comp & (n % 2 == 1)
    n = n[m]
    if len(n) == 0: continue
    e = n - 1; res = np.ones(len(n), np.uint64); b = np.full(len(n), 2, np.uint64) % n
    for bit in range(int(X).bit_length()):
        sel = ((e >> np.uint64(bit)) & np.uint64(1)) == 1
        res = np.where(sel, (res * b) % n, res); b = (b * b) % n
    psp.extend(int(v) for v in n[res == 1])
psp = np.array(psp)
def is_carm(n):
    m = n; ps = []; d = 2
    while d * d <= m:
        if m % d == 0:
            m //= d; ps.append(d)
            if m % d == 0: return False
        d += 1 if d == 2 else 2
    if m > 1: ps.append(m)
    return len(ps) >= 3 and all((n - 1) % (p - 1) == 0 for p in ps)
# every odd Carmichael number passes the base-2 test, so they are found in the pseudoprime list
carm = np.array([v for v in psp if is_carm(int(v))])
pr('X = %d: base-2 pseudoprimes %d, Carmichael numbers %d; first pseudoprimes %s; first Carmichael %s' % (X, len(psp), len(carm), psp[:6].tolist(), carm[:6].tolist()))
for e in range(3, int(math.log10(X)) + 1):
    y = 10 ** e
    pr('  up to 10^%d: pseudoprimes %d, Carmichael %d' % (e, int(np.sum(psp <= y)), int(np.sum(carm <= y))))
np.save('psp_%d.npy' % X, psp)
