"""Brun's ladder: B(x) = sum over twin pairs (p, p + 2) with p + 2 <= x... here p <= x, of 1/p + 1/(p + 2), and the
classical extrapolation B(x) + 4 C2 / ln x (C2 = 0.6601618158 the twin prime constant; Brent, Nicely), against Brun's
constant 1.902160583 (Nicely's estimate from twins below 2 * 10^16; later work to 10^18 agrees).  Prediction written
before the run: the extrapolation stays within 10^-3 of B from 10^6 on, while B(x) itself is still about 0.03 short at
10^9.   usage: brun.py X out.txt"""
import sys, math
import numpy as np
X = int(sys.argv[1]); out = open(sys.argv[2], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
C2 = 0.6601618158468696; BRUN = 1.902160583
r = int(math.isqrt(X + 3)) + 1; s = np.ones(r + 1, bool); s[:2] = False
for i in range(2, int(r ** .5) + 1):
    if s[i]: s[i * i::i] = False
base = np.nonzero(s)[0]
marks = [10 ** e for e in range(2, int(math.log10(X)) + 1)]; tot = 0.0; acc = {}; cnt = 0; cnts = {}
seg = 10 ** 7
for lo in range(0, X + 1, seg):
    hi = min(X + 3, lo + seg + 2); m = np.ones(hi - lo, bool)
    for p in base:
        st = max(p * p, -(-lo // p) * p)
        if st < hi: m[st - lo::p] = False
    if lo == 0: m[:2] = False
    top = min(seg, X + 1 - lo)
    tw = np.nonzero(m[:top] & m[2:top + 2])[0] + lo
    v = 1.0 / tw + 1.0 / (tw + 2)
    cs = np.cumsum(v)
    for mk in marks:
        if lo <= mk < lo + top:
            k = np.searchsorted(tw, mk, 'right'); acc[mk] = tot + (cs[k - 1] if k else 0.0); cnts[mk] = cnt + k
    tot += float(cs[-1]) if len(cs) else 0.0; cnt += len(tw)
pr('     x   twin pairs   B(x)          B(x) + 4 C2 / ln x   minus Brun 1.902160583')
for mk in marks:
    b = acc[mk]; e = b + 4 * C2 / math.log(mk); pr('%8.0e %11d %13.9f %18.9f %+15.2e' % (mk, cnts[mk], b, e, e - BRUN))
