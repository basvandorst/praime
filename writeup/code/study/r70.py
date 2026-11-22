"""Round 70 (theory, no data): how the discrete leading-order law's information between consecutive gaps approaches
I_inf, at lam from 21 to 300 (beyond any data). P(g1, g2) proportional to S({0, g1, g1+g2}) e^(-(g1+g2)/lam), summed to
g1 + g2 <= G = 24 lam. usage: r70.py lam1,lam2,... | fit   (fit: lam (I - I_inf) against a + b ln lam and
a + b ln lam + c ln^2 lam, stored in r70.json under 'fit')"""
import os, sys, json, math
import numpy as np
import theory2 as T2

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
PATH = os.path.join(OUT, 'r70.json')
SH = 1 << 16


def mi(lam):
    G = int(24 * lam) // 2 * 2
    ev = np.arange(2, G + 1, 2, dtype=np.int64)
    num = 0.0
    tot = 0.0; Hj = 0.0; m1 = np.zeros(len(ev)); m2 = np.zeros(len(ev))
    for i0 in range(0, len(ev), 64):
        a = ev[i0:i0 + 64]
        A, B = np.meshgrid(a, ev, indexing='ij'); A = A.ravel(); B = B.ravel(); m = A + B <= G; A, B = A[m], B[m]
        if len(A) == 0:
            continue
        S = np.concatenate([T2.S_many(np.stack([np.zeros(len(A[j:j + 200000]), np.int64), A[j:j + 200000], A[j:j + 200000] + B[j:j + 200000]], 1)) for j in range(0, len(A), 200000)])
        w = S * np.exp(-(A + B) / lam); k = w > 0
        w, A, B = w[k], A[k], B[k]
        tot += w.sum(); Hj += float(np.sum(w * np.log(w)))
        np.add.at(m1, A // 2 - 1, w); np.add.at(m2, B // 2 - 1, w)
    H = lambda x: -float(np.sum((x[x > 0] / tot) * np.log2(x[x > 0] / tot)))
    Hjoint = -(Hj / tot - math.log(tot)) / math.log(2)
    return H(m1) + H(m2) - Hjoint


def fit(res):
    ks = [k for k in res if k[0].isdigit()]
    lam = np.array([float(k) for k in ks]); y = np.array([res[k]['lambda_times_diff'] for k in ks])
    o = np.argsort(lam); lam, y = lam[o], y[o]; L = np.log(lam)
    b, a = np.polyfit(L, y, 1); r1 = y - (a + b * L)
    c2, b2, a2 = np.polyfit(L, y, 2); r2 = y - (a2 + b2 * L + c2 * L * L)
    return {'a_bits': float(a), 'b_bits': float(b), 'max_residual': float(np.max(np.abs(r1))), 'residual_signs': ''.join('+' if v > 0 else '-' for v in r1),
            'quadratic': {'a': float(a2), 'b': float(b2), 'c': float(c2), 'max_residual': float(np.max(np.abs(r2)))},
            'form': 'lam (I - I_inf) against a + b ln lam (and + c ln^2 lam); a description of lam = 21-300, not a law'}


if __name__ == '__main__':
    res = json.load(open(PATH)) if os.path.exists(PATH) else {}
    if sys.argv[1] == 'fit':
        res['fit'] = fit(res); json.dump(res, open(PATH, 'w'), indent=1); print(res['fit']); sys.exit()
    Iinf = json.load(open(os.path.join(OUT, 'r56.json')))['I_inf_bits']
    for lam in map(float, sys.argv[1].split(',')):
        if '%.1f' % lam in res:
            continue
        I = mi(lam)
        res['%.1f' % lam] = {'I_discrete': I, 'diff': I - Iinf, 'lambda_times_diff': lam * (I - Iinf)}
        json.dump(res, open(PATH, 'w'), indent=1)
        print(lam, {k: round(v, 5) for k, v in res['%.1f' % lam].items()}, flush=True)
