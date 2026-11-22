"""v = Var X / E X per class, own code, against the repo's predicted values (r43.json); spot-check rho by brute force."""
import json, math, sys
import numpy as np
import theory as T

R = json.load(open('/home/bas/praime/writeup/data/study/r43.json'))
WINS = {'train': (10 ** 7, 9 * 10 ** 7), '1e8-1e9': (10 ** 8, 9 * 10 ** 8), '1e9-1e10': (10 ** 9, 9 * 10 ** 9),
        'A': (10 ** 10, 5 * 10 ** 7), 'B': (10 ** 12, 5 * 10 ** 7), 'C': (10 ** 15, 2 * 10 ** 7)}


def brute_rho(d, Q):
    out = 0.0
    for p in T.PR:
        p = int(p)
        if Q % p == 0: continue
        nu = len({0, 2 % p, d % p, (d + 2) % p})
        out += math.log((1 - nu / p) / (1 - 2 / p) ** 2)
    return math.exp(out - 4.0 / (T.PMAX * math.log(T.PMAX)))


todo = [(w, Q) for w in WINS for Q in (210, 2310, 30030) if not (Q == 210 and w == '1e9-1e10')]
if len(sys.argv) > 1:
    todo = [t for t in todo if t[0] in sys.argv[1:]]
out = {}
for w, Q in todo:
    L0, W = WINS[w]
    S = T.setup(Q, L0, W)
    F0, rho, d = T.F(S, 0)
    r = S['r']
    v = (float(np.sum(r * (1 - r))) + F0) / float(r.sum())
    pred = R[w][str(Q)]['predicted']
    if w == 'train' and Q == 2310:
        for i in (S['n'] + 1, S['n'] + 12, S['n'] + 168, S['n'] + 220, S['n'] + 1000):
            print('  rho check d=%d own %.10f brute %.10f' % (d[i], rho[i], brute_rho(int(d[i]), Q)))
    out[f'{w} {Q}'] = dict(v=v, repo=pred, rel=v / pred - 1, c=S['c'], rbar=float(r.mean()), n=S['n'])
    print(w, Q, 'own v %.5f repo %.5f rel diff %.2e  c=%.4f rbar=%.4f' % (v, pred, v / pred - 1, S['c'], r.mean()), flush=True)
json.dump(out, open('run_v_%s.json' % ('_'.join(sys.argv[1:]) or 'all'), 'w'), indent=1)
