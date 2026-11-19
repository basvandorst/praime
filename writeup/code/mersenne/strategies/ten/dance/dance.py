"""Last digits of consecutive primes (mod 10, primes > 5) as a dance between 1, 3, 7, 9: the 4x4 transition shares
and the share of repeats (25 % if independent), for the primes and for fakes (true primes up to 100, then numbers n
avoiding the multiples of W with chance prod p/(p-1) / ln n).  Prediction written before the run: fakes avoiding 2
and 5 already show most of the repeat deficit (a repeat needs a gap >= 10); avoiding 3 as well brings them closer.
Windows: all primes in (100, X].   usage: dance.py X runs out.txt"""
import sys
import numpy as np
X, RUNS = int(sys.argv[1]), int(sys.argv[2]); out = open(sys.argv[3], 'w')
def pr(*a):
    s = ' '.join(str(v) for v in a); print(s); out.write(s + '\n')
sv = np.ones(X + 1, bool); sv[:2] = False
for i in range(2, int(X ** .5) + 1):
    if sv[i]: sv[i * i::i] = False
P = np.nonzero(sv)[0]
D = {1: 0, 3: 1, 7: 2, 9: 3}
def mat(seq):
    s = seq[seq > 100]; d = s % 10; ok = (d != 5) & (d % 2 == 1); s = s[ok]; d = np.vectorize(D.get)(s % 10)
    M = np.zeros((4, 4)); np.add.at(M, (d[:-1], d[1:]), 1); return M / M.sum()
Mp = mat(P); rep = np.trace(Mp)
pr('X = %d; primes: repeats %.4f (25%% if independent); matrix rows from 1,3,7,9 to 1,3,7,9:' % (X, rep))
for i in range(4): pr('   ' + ' '.join('%.4f' % v for v in Mp[i]))
rng = np.random.default_rng(13); m = np.arange(X + 1, dtype=float)
for W in ((2, 5), (2, 3, 5), (2, 3, 5, 7), (2, 3, 5, 7, 11, 13), tuple(int(p) for p in P[P <= 100])):
    W = np.array(W); ok = np.ones(X + 1, bool); ok[:101] = False
    for p in W: ok[::p] = False
    c = float(np.prod(W / (W - 1.0))); pb = np.zeros(X + 1); pb[ok] = np.minimum(1, c / np.log(m[ok]))
    reps = []; Ms = []
    for r in range(RUNS):
        f = np.concatenate([P[P <= 100], np.nonzero(rng.random(X + 1) < pb)[0]]); M = mat(f); Ms.append(M); reps.append(np.trace(M))
    Mm = np.mean(Ms, 0)
    pr('fakes avoiding %-26s: repeats %.4f (spread %.4f over %d runs); largest cell difference to the primes %.4f' % ('2..%d' % W.max() if len(W) > 6 else tuple(W.tolist()), np.mean(reps), np.std(reps), RUNS, np.abs(Mm - Mp).max()))
