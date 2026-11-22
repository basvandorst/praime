"""Round 76 (preregister/round76.txt): round 74's theory chain on the classes mod 510510 (y = 17), to test whether the
lag of the y = 13 chain is the memory from q >= 17. The reveal time is still that of p mod 30030 (q = 3..13).
usage: R47_Y=17 r76.py table | chain   (table resumable, stored in STUDY_SCRATCH)"""
import os, sys, json, math, time
import numpy as np
os.environ['R47_Y'] = '17'
import r47

assert r47.Y == 17
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
Q = r47.PP
CLS = np.array([a for a in range(Q) if math.gcd(a, Q) == 1], dtype=np.int64)
IDX = np.full(Q, -1); IDX[CLS] = np.arange(len(CLS))
LAM = math.log(1.05e9)
G = int(20 * LAM) // 2 * 2
TAB = os.path.join(os.environ.get('STUDY_SCRATCH', '/tmp'), 'r76_table_y17.npy')
KMAX = 150; NCH = 200_000; QS = (3, 5, 7, 11, 13)


def per_class(g):
    st = r47.setup(g)
    rho0 = 1 / (LAM * r47.RHO_DEN); p = rho0 * st['s1']; M = st['M']
    L = M @ np.log1p(-p) if M.shape[1] else np.zeros(M.shape[0])
    if M.shape[1] > 1:
        w = p / (1 - p)
        L = L + 0.5 * np.sum((M @ ((w[:, None] * w[None, :]) * st['R'])) * M, axis=1)
    cls = CLS[np.gcd(CLS + g, Q) == 1]
    out = np.zeros(len(CLS), np.float32); out[IDX[cls]] = rho0 * st['S2'] * np.exp(L)
    return out


if __name__ == '__main__':
    if sys.argv[1] == 'table':
        T = np.load(TAB) if os.path.exists(TAB) else np.full((len(CLS), G // 2), np.nan, np.float32)
        t0 = time.time()
        for j, g in enumerate(range(2, G + 1, 2)):
            if not np.isnan(T[0, j]) or time.time() - t0 > 420:
                continue
            T[:, j] = per_class(g); np.save(TAB, T)
        done = int(np.sum(~np.isnan(T[0])))
        print('columns done', done, 'of', G // 2, 'row sums', float(np.nanmin(np.nansum(T, 1))), float(np.nanmax(np.nansum(T, 1))), flush=True)
    else:
        T = np.load(TAB).astype(np.float64); rs = T.sum(1)
        cum = np.cumsum(T / rs[:, None], 1); cum[:, -1] = 1.0
        gaps = np.arange(2, G + 1, 2); rng = np.random.default_rng(76)
        Tall = []; Tq = {q: [] for q in QS}; rep = {q: 0 for q in QS}; nst = 0
        for b0 in range(0, NCH, 5_000):
            n = 5_000; cur = CLS[rng.integers(0, len(CLS), n)]; D = np.zeros(n, np.int64)
            seen = {q: np.zeros((n, q), bool) for q in QS}; fq = {q: np.full(n, -1) for q in QS}
            for q in QS:
                seen[q][:, 0] = True
            for k in range(1, KMAX + 1):
                row = IDX[(-cur) % Q]; u = rng.random(n)
                j = np.argmax(cum[row] >= u[:, None], axis=1); g = gaps[j]
                for q in QS:
                    rep[q] += int(np.sum(g % q == 0))
                nst += n; D += g; cur = (cur - g) % Q
                for q in QS:
                    seen[q][np.arange(n), D % q] = True
                    fq[q][(seen[q].sum(1) == q - 1) & (fq[q] < 0)] = k
            st = np.stack([fq[q] for q in QS]); t = np.max(st, 0); t[np.any(st < 0, 0)] = KMAX + 1
            Tall.append(t)
            for q in QS:
                Tq[q].append(np.where(fq[q] < 0, KMAX + 1, fq[q]))
        Tall = np.concatenate(Tall)
        d73 = json.load(open(os.path.join(OUT, 'r73.json'))); r74 = json.load(open(os.path.join(OUT, 'r74.json')))
        cd = np.array(d73['cdf_data']); cdf = np.array([np.mean(Tall <= k) for k in range(KMAX + 1)])
        res = {'y': 17, 'classes': len(CLS), 'G': G, 'row_sum_min': float(rs.min()), 'row_sum_max': float(rs.max()),
               'mean_T_chain_y17': float(np.mean(Tall)), 'sem': float(np.std(Tall) / math.sqrt(len(Tall))), 'unresolved': int(np.sum(Tall > KMAX)),
               'mean_T_chain_y13': r74['mean_T_chain'], 'mean_T_data': d73['mean_T_data'],
               'fraction_of_lag_closed': (r74['mean_T_chain'] - float(np.mean(Tall))) / (r74['mean_T_chain'] - d73['mean_T_data']),
               'max_cdf_diff_vs_data': float(np.max(np.abs(cdf - cd))),
               **{'mean_T_q%d' % q: float(np.mean(np.concatenate(Tq[q]))) for q in QS},
               **{'P_gap0_q%d' % q: rep[q] / nst for q in QS}, 'cdf': cdf.tolist()}
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in res.items() if k != 'cdf'}, flush=True)
        json.dump(res, open(os.path.join(OUT, 'r76.json'), 'w'), indent=1)
