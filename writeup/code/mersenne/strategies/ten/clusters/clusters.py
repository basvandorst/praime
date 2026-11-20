"""Prime star clusters: the densest admissible patterns (twins (0,2); triplets (0,2,6) and (0,4,6); quadruplets (0,2,6,8);
quintuplets (0,2,6,8,12) and (0,4,6,10,12); sextuplets (0,4,6,10,12,16)), counted by first member p <= x, against
Hardy-Littlewood: C * int_2^x dt / ln^k t, C = prod_p (1 - nu(p)/p) / (1 - 1/p)^k (nu(p) = residues the pattern covers
mod p; product over p < 10^6).  Segmented sieve up to X.  Prediction written before the run: every ratio count/prediction
within 2 % of 1 at 10^9; the two mirror patterns race fairly (the leader changes, the difference ~ sqrt of the count).
usage: clusters.py X out.txt"""
import sys, math
import numpy as np
from scipy.integrate import quad
X = int(sys.argv[1]); out = open(sys.argv[2], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
PAT = {'twins': (0, 2), 'triplet a': (0, 2, 6), 'triplet b': (0, 4, 6), 'quadruplet': (0, 2, 6, 8), 'quintuplet a': (0, 2, 6, 8, 12), 'quintuplet b': (0, 4, 6, 10, 12), 'sextuplet': (0, 4, 6, 10, 12, 16)}
M = 10 ** 6; sp = np.ones(M + 1, bool); sp[:2] = False
for i in range(2, int(M ** .5) + 1):
    if sp[i]: sp[i * i::i] = False
SP = np.nonzero(sp)[0]
def const(pat):
    k = len(pat); c = 1.0
    for p in SP:
        nu = len(set(h % p for h in pat)); c *= (1 - nu / p) / (1 - 1 / p) ** k
    return c
C = {n: const(p) for n, p in PAT.items()}
r = int(math.isqrt(X + 20)) + 1; base = SP[SP <= r]
marks = [10 ** e for e in range(3, int(math.log10(X)) + 1)]
cnt = {n: {m: 0 for m in marks} for n in PAT}; race = {'triplet': [0, 0, 0], 'quintuplet': [0, 0, 0]}   # a - b, changes, last sign
seg = 10 ** 7; first = {n: [] for n in PAT}
for lo in range(0, X + 1, seg):
    hi = min(X + 17, lo + seg + 17); m = np.ones(hi - lo, bool)
    for p in base:
        st = max(p * p, -(-lo // p) * p)
        m[st - lo::p] = False
    if lo == 0: m[:2] = False
    top = min(seg, X + 1 - lo)
    ev = {}
    for n, pat in PAT.items():
        ok = m[:top].copy()
        for h in pat[1:]: ok &= m[h:h + top]
        idx = np.nonzero(ok)[0] + lo
        if len(first[n]) < 6: first[n] += [int(v) for v in idx[:6 - len(first[n])]]
        for mm in marks: cnt[n][mm] += int(np.sum(idx <= mm))
        ev[n] = idx
    for key, (na, nb) in (('triplet', ('triplet a', 'triplet b')), ('quintuplet', ('quintuplet a', 'quintuplet b'))):
        e = np.concatenate([np.stack([ev[na], np.ones(len(ev[na]))], 1), np.stack([ev[nb], -np.ones(len(ev[nb]))], 1)])
        if len(e) == 0: continue
        e = e[np.argsort(e[:, 0], kind='stable')]; d = race[key][0] + np.cumsum(e[:, 1]); s = np.sign(d); s = s[s != 0]
        if len(s):
            prev = race[key][2]; ss = np.concatenate([[prev], s]) if prev else s
            race[key][1] += int(np.sum(ss[1:] != ss[:-1])); race[key][2] = int(s[-1])
        race[key][0] = int(d[-1])
pr('X = %d; Hardy-Littlewood constants: %s' % (X, ', '.join('%s %.6f' % (n, c) for n, c in C.items())))
for n, pat in PAT.items():
    k = len(pat)
    row = []
    for mm in marks:
        I = quad(lambda t: 1 / math.log(t) ** k, 2, mm, limit=400)[0]; row.append('%d/%.0f=%.4f' % (cnt[n][mm], C[n] * I, cnt[n][mm] / (C[n] * I)))
    pr('%-13s %s first %s' % (n, '  '.join(row[-4:]), first[n]))
for key, (d, ch, s) in race.items():
    pr('race %s a - b at X: %+d (sqrt of the count %.0f); leader changes %d' % (key, d, math.sqrt(cnt[key + ' a'][marks[-1]]), ch))
