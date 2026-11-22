"""Null test for the k-th neighbour law (r81): fake prime sequences in [1e9, 1.1e9).
models: cramer | sieve13 | sieveZ (Z = 1000, 10000)
usage: fakes.py MODEL [seed]"""
import sys, json, math, time
import numpy as np

L0, L1, SEG = 10**9, 11 * 10**8, 10**7
KS = (1, 2, 3, 4); QS = (3, 5, 7, 11, 13); DMAX = 240
MODEL = sys.argv[1]; SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 1
rng = np.random.default_rng(SEED)

def small_primes(z):
    s = np.ones(z + 1, bool); s[:2] = False
    for i in range(2, int(z ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]

if MODEL == 'cramer':
    SP = np.array([2]); DEN = 0.5
elif MODEL == 'sieve13':
    SP = small_primes(13); DEN = np.prod(1 - 1 / SP)
elif MODEL.startswith('sieve'):
    z = int(MODEL[5:]); SP = small_primes(z); DEN = np.prod(1 - 1 / SP.astype(float))
else:
    raise SystemExit('model?')
print(MODEL, 'sieve primes', len(SP), 'prod(1-1/p)', DEN, 'rho at mid', 1 / (math.log(1.05e9) * DEN), flush=True)

cnt = {k: {q: 0 for q in QS} for k in KS}; nin = {k: 0 for k in KS}; ntot = {k: 0 for k in KS}
carry = np.zeros(0, np.int64); nprimes = 0; t0 = time.time()
for s0 in range(L0, L1, SEG):
    m = np.ones(SEG, bool)
    for p in SP:
        p = int(p); m[(-s0) % p::p] = False
    idx = np.nonzero(m)[0].astype(np.int64)
    n = idx + s0
    keep = rng.random(len(n)) < (SEG / len(n)) / np.log(n.astype(float))   # exact density 1/ln n per segment
    P = np.concatenate([carry, n[keep]]); nprimes += int(keep.sum())
    for k in KS:
        d = P[k:] - P[:-k]
        if len(carry):
            d = d[max(0, len(carry) - k):]
        ntot[k] += len(d); d = d[d <= DMAX]; nin[k] += len(d)
        for q in QS:
            cnt[k][q] += int(np.count_nonzero(d % q == 0))
    carry = P[-max(KS):]
print('fake primes', nprimes, round(time.time() - t0), 's', flush=True)
json.dump({'model': MODEL, 'seed': SEED, 'nprimes': nprimes, 'cnt': {str(k): {str(q): v for q, v in c.items()} for k, c in cnt.items()},
           'nin': {str(k): v for k, v in nin.items()}, 'ntot': {str(k): v for k, v in ntot.items()}},
          open('/home/bas/scratch/null-D/out_%s_%d.json' % (MODEL, SEED), 'w'))
