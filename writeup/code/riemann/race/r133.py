"""Round 133 (preregister/round133.txt): class 1 last in the races mod 8 and mod 12, from the zeros of the three real
characters (gamma <= 50) and their prime-squares terms; measured from a sieve to 1e7 (counts from p = 3 on; a first run counted only from 23, a bug). usage: r133.py"""
import os, sys, json, math
import numpy as np, cypari2
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'study'))
import common as C
pari = cypari2.Pari(); pari.allocatemem(800 * 10 ** 6)
H = os.path.dirname(os.path.abspath(__file__))
P = C.small_primes(10 ** 7).astype(np.int64); P = P[P >= 23]; lnP = np.log(P.astype(float))
out = {}
for q, Ds in ((8, (-4, 8, -8)), (12, (-4, -3, 12))):
    W = {}
    for D in Ds:
        g = np.array([float(t) for t in pari.lfunzeros(pari.lfuncreate(D), 50) if float(t) > 0])
        W[D] = 1 + 2 * np.real((np.exp(1j * np.outer(lnP, g)) / (0.5 + 1j * g)[None, :]).sum(axis=1))
    units = [a for a in range(2, q) if math.gcd(a, q) == 1]
    last_pred = np.ones(len(P), bool)
    for a in units:
        s = sum(W[D] for D in Ds if int(pari.kronecker(D, a)) == -1)
        last_pred &= s > 0
    Pall = C.small_primes(10 ** 7).astype(np.int64); Pall = Pall[Pall >= 3]; keep = Pall >= 23
    cnt = {a: np.cumsum(Pall % q == a)[keep] for a in [1] + units}
    last_meas = np.ones(len(P), bool)
    for a in units: last_meas &= cnt[1] < cnt[a]
    out[str(q)] = {'pred_share': float(last_pred.mean()), 'meas_share': float(last_meas.mean())}
    print('mod', q, out[str(q)], flush=True)
rd = [all(v['pred_share'] >= 0.99 for v in out.values()), all(v['meas_share'] == 1.0 for v in out.values())]
print('readings', rd)
out['readings'] = rd
json.dump(out, open(os.path.join(H, 'r133.json'), 'w'), indent=1)
