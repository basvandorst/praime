"""Round 41 (preregister/round41.txt): round 40 in the new window [1e8, 1e9) (W/x = 9 as at train, x ten times larger).
Counts by exact segmented sieve (segments of 1e7, the gap across a segment boundary kept). usage: r41.py"""
import os, json, math, pickle
import numpy as np
import common as C
import theory_r3 as R3
import r7
import r13
import r39
import r40

Q = 210
L0, W = 10 ** 8, 9 * 10 ** 8
GAPS = r40.GAPS


def counts(U):
    idx = np.full(Q, -1); idx[np.array(U)] = np.arange(len(U))
    N = np.zeros((48, 20)); Y = np.zeros(48); Xg = {g: np.zeros(48) for g in GAPS}
    last = None
    for s in range(0, W, 10 ** 7):
        P = np.nonzero(C.primes_mask(L0 + s, 10 ** 7))[0].astype(np.int64) + L0 + s
        if last is not None:
            P = np.concatenate([[last], P])
        last = int(P[-1])
        a = idx[P[:-1] % Q]; g = np.diff(P)
        np.add.at(N, (a, np.minimum(g // 2, 20) - 1), 1)
        Y += np.bincount(a, minlength=48)
        for gg in GAPS:
            Xg[gg] += np.bincount(a[g == gg], minlength=48)
    return N, Y, Xg


if __name__ == '__main__':
    cache = pickle.load(open(r7.CACHE7, 'rb'))
    U = r7.units(Q)
    T = r7.trans(cache, Q, L0, W); Pa = T.sum(2).T; Pa = Pa / Pa.sum(1, keepdims=True)
    Pa_rows = {g: [float(Pa[i, g // 2 - 1]) if math.gcd(int(a) + g, Q) == 1 else 0.0 for i, a in enumerate(U)] for g in GAPS}
    N, Y, Xg = counts(U)
    c, d = r13.chi2(N, Pa)
    pr = r40.predict(L0, W, Pa_rows)
    Ybar = Y.mean(); v = (Y - Ybar) / math.sqrt(Ybar)
    per = {}; pu = []; mu = []; pc = []; mc = []; wts = []; cell_pred = 0.0; cell_null = 0.0
    for g in GAPS:
        rows = [i for i in range(48) if Pa_rows[g][i] > 0]
        u = np.array([(Xg[g][i] - Pa_rows[g][i] * Ybar) / math.sqrt(Pa_rows[g][i] * Ybar) for i in rows])
        d_ = {'meas_u2': float(np.mean(u * u)), 'pred_u2': float(np.mean([pr[g]['vX'][i] for i in rows])),
              'meas_uv': float(np.mean(u * v[rows])), 'pred_uv': float(np.mean([pr[g]['cXY'][i] for i in rows])), 'rows': len(rows)}
        per[str(g)] = d_
        pu.append(d_['pred_u2']); mu.append(d_['meas_u2']); pc.append(d_['pred_uv']); mc.append(d_['meas_uv']); wts.append(len(rows))
        for i in rows:                                   # round-39 cell prediction from the parts
            P = Pa_rows[g][i]
            cell_pred += pr[g]['vX'][i] - 2 * math.sqrt(P) * pr[g]['cXY'][i] + P * pr['vY']; cell_null += 1 - P
    P40 = [p for p in Pa[:, 19] if p > 0]
    implied = (cell_pred + sum(1 - p for p in P40)) / (cell_null + sum(1 - p for p in P40))
    dl = np.array(mu) - np.array(pu)
    res = {'window': [L0, W], 'per_gap': per, 'primes': int(Y.sum()),
           'u2': {'pred': float(np.average(pu, weights=wts)), 'meas': float(np.average(mu, weights=wts)),
                  'diff': float(np.average(dl, weights=wts)), 'se_from_scatter': float(np.std(dl, ddof=1) / math.sqrt(len(dl)))},
           'uv': {'pred': float(np.average(pc, weights=wts)), 'meas': float(np.average(mc, weights=wts))},
           'v2': {'pred': pr['vY'], 'meas': float(np.mean(v * v) * 48 / 47)},
           'chi2_per_dof': {'measured': c / d, 'dof': d, 'implied_prediction': implied}}
    print(json.dumps({k: res[k] for k in ('primes', 'u2', 'uv', 'v2', 'chi2_per_dof')}), flush=True)
    json.dump(res, open(os.path.join(R3.OUT, 'r41.json'), 'w'), indent=1)
