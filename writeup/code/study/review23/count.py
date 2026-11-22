"""Independent recount: P(p_{n+k} = p_n mod q | p_{n+k} - p_n <= DMAX), pairs with p_n, p_{n+k} in [L0, L1).
Own segmented sieve; pairs indexed by their FIRST prime with a look-ahead (r81 indexes by the second, with a carry).
usage: count.py L0 L1 DMAX K q1,q2,... [segsize] [outfile]"""
import sys, math, json, time
import numpy as np

L0 = int(float(sys.argv[1])); L1 = int(float(sys.argv[2])); DMAX = int(sys.argv[3]); K = int(sys.argv[4])
QS = [int(x) for x in sys.argv[5].split(',')]
SEG = int(float(sys.argv[6])) if len(sys.argv) > 6 else 20_000_000
OUTF = sys.argv[7] if len(sys.argv) > 7 else None
LOOK = 4000

r = math.isqrt(L1) + 2
s = np.ones(r + 1, bool); s[:2] = False
for i in range(2, math.isqrt(r) + 1):
    if s[i]:
        s[i * i::i] = False
BASE = np.nonzero(s)[0]


def sieve(lo, hi):
    """primes in [lo, hi), odd only trick not used: plain bool array"""
    m = np.ones(hi - lo, bool)
    for p in BASE:
        p = int(p)
        if p * p >= hi:
            break
        st = max(p * p, ((lo + p - 1) // p) * p)
        m[st - lo::p] = False
    if lo < 2:
        m[:2 - lo] = False
    return np.nonzero(m)[0].astype(np.int64) + lo


cnt = np.zeros((K + 1, len(QS)), np.int64); nin = np.zeros(K + 1, np.int64); ntot = np.zeros(K + 1, np.int64)
hist = {k: np.zeros(DMAX // 2 + 1, np.int64) for k in range(1, K + 1)}
t0 = time.time()
for lo in range(L0, L1, SEG):
    hi = min(lo + SEG, L1); hi2 = min(hi + LOOK, L1)
    P = sieve(lo, hi2)
    nfirst = int(np.searchsorted(P, hi))          # primes with p_n in [lo, hi)
    assert hi2 == L1 or len(P) - nfirst >= K, 'look-ahead too short'
    for k in range(1, K + 1):
        last = min(nfirst, len(P) - k)            # p_n index < nfirst and p_{n+k} exists in window
        d = P[k:k + last] - P[:last]
        ntot[k] += len(d)
        d = d[d <= DMAX]; nin[k] += len(d)
        hist[k] += np.bincount(d // 2, minlength=DMAX // 2 + 1)
        for j, q in enumerate(QS):
            cnt[k, j] += int(np.count_nonzero(d % q == 0))
out = {'L0': L0, 'L1': L1, 'DMAX': DMAX, 'res': {}}
for k in range(1, K + 1):
    row = {}
    for j, q in enumerate(QS):
        m = cnt[k, j] / nin[k]
        row[str(q)] = {'p': m, 'se': math.sqrt(m * (1 - m) / nin[k])}
    out['res'][str(k)] = {'n_in': int(nin[k]), 'n_tot': int(ntot[k]), 'mass': nin[k] / ntot[k], 'rep': row,
                          'hist_even_D': hist[k].tolist()}
    print(k, nin[k], ntot[k], ' '.join('%d:%.5f' % (q, row[str(q)]['p']) for q in QS), flush=True)
print('secs', round(time.time() - t0))
if OUTF:
    json.dump(out, open(OUTF, 'w'))
