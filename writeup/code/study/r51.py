"""Round 51 (preregister/round51.txt): the first occurrences of all 138 gaps 180..454 with the conditioned law.
usage: r51.py o2 | o3 | summary   (o2/o3 resumable, ~7 minutes per call)"""
import os, sys, json, math, time
import numpy as np
from scipy.stats import kstest
import r47

OUT = r47.OUT
PATH = os.path.join(OUT, 'r51.json')
L2 = np.arange(16.0, 27.01, 1.0)
L3NEW = np.array([18.0, 22.0, 26.0])
P = {int(g): v for g, v in json.load(open(os.path.join(OUT, 'oes_first_occurrence.json')))['P'].items()}
GS = sorted(g for g in P if 180 <= g <= 454)


def o3_delta(g, lams):
    st = r47.setup(g)
    pts = np.arange(2, g, 2, dtype=np.int64); n = len(pts)
    s1, R, M = st['s1'], st['R'], st['M']
    K = np.zeros((n, n, n), dtype=np.float32)
    for i in range(n - 2):
        jj, kk = np.triu_indices(n - i - 1, 1); jj = jj + i + 1; kk = kk + i + 1
        z = np.zeros(len(jj), np.int64)
        s3 = r47.S_gt(np.stack([z, np.full(len(jj), g), np.full(len(jj), pts[i]), pts[jj], pts[kk]], 1)) / st['S2']
        kap = s3 / (s1[i] * s1[jj] * s1[kk]) - 1 - R[i, jj] - R[i, kk] - R[jj, kk]
        K[i, jj, kk] = kap; K[i, kk, jj] = kap
    d = []
    for lam in lams:
        rho0 = 1 / (lam * r47.RHO_DEN); p = rho0 * s1; w = p / (1 - p)
        Lb = M @ np.log1p(-p)
        Cm = w[:, None] * w[None, :] * R
        Lo2 = Lb + 0.5 * np.sum((M @ Cm) * M, axis=1)
        W = (M * w[None, :]).astype(np.float32)
        t3 = np.zeros(M.shape[0])
        for i in range(n - 2):
            Wi = W[:, i + 1:]
            t3 += W[:, i] * 0.5 * np.sum((Wi @ K[i, i + 1:, i + 1:]) * Wi, axis=1)
        v = M @ Cm; v2 = M @ (Cm * Cm)
        ov = -0.5 * (0.5 * np.sum(v2 * M, axis=1) + np.sum(M * (v * v - v2), axis=1))
        d.append(float(math.log(np.sum(np.exp(Lo2 - t3 + ov)) / np.sum(np.exp(Lo2)))))
    return d


if __name__ == '__main__':
    res = json.load(open(PATH)) if os.path.exists(PATH) else {'o2': {}, 'o3': {}}
    t0 = time.time(); mode = sys.argv[1]
    if mode in ('o2', 'o3'):
        for g in GS:
            key = str(g)
            if key in res[mode] or time.time() - t0 > 400:
                continue
            if mode == 'o2':
                st = r47.setup(g); res['o2'][key] = [math.log(r47.share(st, lam, 2)) for lam in L2]
            else:
                if g <= 300:
                    u3 = json.load(open(os.path.join(OUT, 'r47_u3_y13.json')))
                    res['o3'][key] = {'lams': u3['lams'], 'delta': u3['delta_ln'][key]}
                else:
                    res['o3'][key] = {'lams': L3NEW.tolist(), 'delta': o3_delta(g, L3NEW)}
            json.dump(res, open(PATH, 'w'), indent=1)
            print(mode, g, round(time.time() - t0), flush=True)
        print(mode, 'done', sum(str(g) in res[mode] for g in GS), 'of', len(GS), flush=True)
    else:
        out = {}
        for name, use3 in (('order2', False), ('order3_overlap', True)):
            U = []
            for g in GS:
                u = np.linspace(math.log(1e7), math.log(P[g]), 4001)
                ls = np.interp(u, L2, res['o2'][str(g)])
                if use3:
                    o = res['o3'][str(g)]; ls = ls + np.interp(u, o['lams'], o['delta'])
                U.append(float(np.trapezoid(np.exp(ls) / u * np.exp(u), u)))
            U = np.array(U); garr = np.array(GS); ks = kstest(U, 'expon')
            out[name] = {'n': len(U), 'mean': float(U.mean()), 'null_se': 1 / math.sqrt(len(U)), 'z_null': float((U.mean() - 1) * math.sqrt(len(U))),
                         'ks_p': float(ks.pvalue), 'mean_180_300': float(U[garr <= 300].mean()), 'mean_302_454': float(U[garr > 300].mean()),
                         'mean_six': float(U[garr % 6 == 0].mean()), 'mean_not_six': float(U[garr % 6 != 0].mean()), 'U': U.tolist()}
            print(name, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in out[name].items() if k != 'U'}, flush=True)
        res['summary'] = out
        json.dump(res, open(PATH, 'w'), indent=1)
