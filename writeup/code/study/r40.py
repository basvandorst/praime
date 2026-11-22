"""Round 40 (diagnostic, preregister/round40.txt): the three parts of round 39's cell prediction, measured across the
residue rows mod 210. u_a = (X_a - P_a Ybar)/sqrt(P_a Ybar), v_a = (Y_a - Ybar)/sqrt(Ybar); predicted
E u^2 = Var X/(P m), E v^2 = Var Y/m, E u v = Cov(X, Y)/(m sqrt P) (common-mode terms of order 1/48 neglected).
usage: r40.py"""
import os, json, math, pickle
import numpy as np
from scipy.signal import fftconvolve
import common as C
import theory_r3 as R3
import r7
import r39

Q = 210
GAPS = list(range(2, 39, 2))


def predict(L0, W, Pa_rows):
    n = W // Q
    r = (Q / 48) / np.log(L0 + Q * np.arange(n) + Q / 2)
    w = fftconvolve(r, r[::-1])[n - 1:]
    m = float(r.sum()); s1 = float(np.sum(r * (1 - r))); s2 = float(np.sum(r * r))
    varY = s1 + 2 * float(np.sum(w[1:] * (r39.rho(n, 1, 'YY')[1:] - 1)))
    out = {'vY': varY / m}
    for g in GAPS:
        SXX = 2 * float(np.sum(w[1:] * (r39.rho(n, g, 'XX')[1:] - 1)))
        SXY = float(np.sum(w[1:] * ((r39.rho(n, g, 'XYp')[1:] - 1) + (r39.rho(n, g, 'XYm')[1:] - 1))))
        vx, cx = [], []
        for P in Pa_rows[g]:
            if P <= 0:
                vx.append(None); cx.append(None); continue
            varX = P * m - P * P * s2 + P * P * SXX
            cov = P * s1 + P * SXY
            vx.append(varX / (P * m)); cx.append(cov / (m * math.sqrt(P)))
        out[g] = {'vX': vx, 'cXY': cx}
    return out


def measure(L0, W, U, Pa_rows):
    M = C.primes_mask(L0, W); P = np.nonzero(M)[0].astype(np.int64) + L0
    a = P[:-1] % Q; g = np.diff(P)
    idx = {int(x): i for i, x in enumerate(U)}
    ai = np.array([idx[int(x)] for x in a])
    Y = np.bincount(ai, minlength=48).astype(float)
    Ybar = Y.mean(); v = (Y - Ybar) / math.sqrt(Ybar)
    out = {'v2': float(np.mean(v * v) * 48 / 47)}
    for gg in GAPS:
        X = np.bincount(ai[g == gg], minlength=48).astype(float)
        rows = [i for i in range(48) if Pa_rows[gg][i] > 0]
        u = np.array([(X[i] - Pa_rows[gg][i] * Ybar) / math.sqrt(Pa_rows[gg][i] * Ybar) for i in rows])
        out[gg] = {'rows': rows, 'u2': (u * u).tolist(), 'uv': (u * v[rows]).tolist(), 'v2_rows': (v[rows] ** 2).tolist()}
    return out


if __name__ == '__main__':
    cache = pickle.load(open(r7.CACHE7, 'rb'))
    U = r7.units(Q)
    res = {}
    for win in ('train', 'A'):
        L0, W = R3.WIN[win]
        T = r7.trans(cache, Q, L0, W); Pa = T.sum(2).T; Pa = Pa / Pa.sum(1, keepdims=True)
        Pa_rows = {g: [float(Pa[i, g // 2 - 1]) if math.gcd(int(a) + g, Q) == 1 else 0.0 for i, a in enumerate(U)] for g in GAPS}
        pr = predict(L0, W, Pa_rows); ms = measure(L0, W, U, Pa_rows)
        per = {}
        pu = []; mu = []; pc = []; mc = []; pv = []; mv = []
        for g in GAPS:
            rows = ms[g]['rows']
            d = {'pred_u2': float(np.mean([pr[g]['vX'][i] for i in rows])), 'meas_u2': float(np.mean(ms[g]['u2'])),
                 'pred_uv': float(np.mean([pr[g]['cXY'][i] for i in rows])), 'meas_uv': float(np.mean(ms[g]['uv'])),
                 'pred_v2': pr['vY'], 'meas_v2_rows': float(np.mean(ms[g]['v2_rows'])), 'rows': len(rows)}
            per[str(g)] = d
            for lst, k in ((pu, 'pred_u2'), (mu, 'meas_u2'), (pc, 'pred_uv'), (mc, 'meas_uv')):
                lst.append(d[k])
        wts = np.array([per[str(g)]['rows'] for g in GAPS], float)
        summ = {}
        for name, p_, m_ in (('u2', pu, mu), ('uv', pc, mc)):
            p_ = np.array(p_); m_ = np.array(m_); dlt = m_ - p_
            summ[name] = {'pred': float(np.average(p_, weights=wts)), 'meas': float(np.average(m_, weights=wts)),
                          'diff': float(np.average(dlt, weights=wts)),
                          'se_from_scatter': float(np.std(dlt, ddof=1) / math.sqrt(len(dlt)))}
        summ['v2'] = {'pred': pr['vY'], 'meas': ms['v2']}
        res[win] = {'per_gap': per, 'summary': summ}
        print(win, json.dumps({k: {kk: round(vv, 3) for kk, vv in v.items()} for k, v in summ.items()}), flush=True)
    json.dump(res, open(os.path.join(R3.OUT, 'r40.json'), 'w'), indent=1)
