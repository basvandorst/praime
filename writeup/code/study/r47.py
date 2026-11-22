"""Round 47 (own method): the gap law by a CONDITIONED expansion.
S(B) = S_{<=y}(B) * S_{>y}(B) exactly, and by the Chinese remainder theorem S_{<=y}(B)/lam^k = E_a[prod_b 1(a+b coprime
to P)] rho0^k, P = prod_{q<=y} q, rho0 = 1/(lam prod_{q<=y}(1 - 1/q)). Hence, exactly,
  P(gap = h | p prime) = E_{a coprime to P} [ 1(a+h coprime) rho0 S_{>y}({0,h}) Pi_a ],
  Pi_a = sum_{B subset R_a} (-1)^|B| rho0^|B| S_{>y}({0,h} u B) / S_{>y}({0,h}),
R_a = the interior points rough in class a. Pi_a is an avoidance probability with point chances rho0 s_i and only the
weak correlations of the primes q > y (at most about 1/(q-2) <= 1/15 for y = 13). It is expanded in Bahadur form:
  ln Pi_a = sum_{i in R_a} ln(1 - p_i) + sum_{i<j in R_a} p_i p_j r_ij / ((1 - p_i)(1 - p_j)) [+ order 3],
p_i = rho0 s_i, s_i = S_{>y}({0,h,i})/S_{>y}({0,h}), r_ij = S_{>y}({0,h,i,j}) S_{>y}({0,h}) / (S_{>y}({0,h,i}) S_{>y}({0,h,j})) - 1.
The classes are averaged exactly (y = 13: 5760 classes). usage: r47.py check | u (first occurrences of the gaps 180..300; resumable, env R47_Y) | o3 g1,g2,... (size of order 3) | u3 (order 3 + overlap correction, resumable)"""
import os, sys, json, math
from itertools import combinations
import numpy as np
import common as C

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
PR = C.small_primes(20_000_000).astype(np.int64)
Y = int(os.environ.get('R47_Y', 13))
SMALL = [int(q) for q in PR if q <= Y]
PP = math.prod(SMALL)
MID = [int(q) for q in PR if Y < q <= 400]
BIG = PR[PR > 400].astype(np.float64)
TAIL = {k: math.exp(float(np.sum(np.log1p(-k / BIG) - k * np.log1p(-1 / BIG)))) for k in range(1, 7)}
RHO_DEN = math.prod(1 - 1 / q for q in SMALL)


def S_gt(sets):
    """S_{>y} of each row of an int array (N, k): product over primes y < q <= 400 with the actual nu, tail beyond."""
    sets = np.asarray(sets, dtype=np.int64)
    N, k = sets.shape
    out = np.full(N, TAIL[k])
    for q in MID:
        r = np.sort(sets % q, axis=1)
        nu = 1 + np.count_nonzero(np.diff(r, axis=1), axis=1)
        out *= (1 - nu / q) / (1 - 1 / q) ** k
    return out


def setup(h):
    pts = np.arange(2, h, 2, dtype=np.int64)
    n = len(pts)
    S2 = float(S_gt(np.array([[0, h]]))[0])
    s1 = S_gt(np.stack([np.zeros(n, np.int64), np.full(n, h), pts], 1)) / S2
    ij = np.array(list(combinations(range(n), 2)), dtype=np.int64) if n > 1 else np.zeros((0, 2), np.int64)
    R = np.zeros((n, n))
    if len(ij):
        s2 = S_gt(np.stack([np.zeros(len(ij), np.int64), np.full(len(ij), h), pts[ij[:, 0]], pts[ij[:, 1]]], 1)) / S2
        rij = s2 / (s1[ij[:, 0]] * s1[ij[:, 1]]) - 1
        R[ij[:, 0], ij[:, 1]] = rij; R[ij[:, 1], ij[:, 0]] = rij
    a = np.array([x for x in range(PP) if math.gcd(x, PP) == 1], dtype=np.int64)
    a = a[np.gcd(a + h, PP) == 1]
    M = (np.gcd(a[:, None] + pts[None, :], PP) == 1).astype(np.float64)      # rough interior points per class
    return {'h': h, 'S2': S2, 's1': s1, 'R': R, 'M': M, 'nclass_all': int(sum(1 for x in range(PP) if math.gcd(x, PP) == 1))}


def share(st, lam, order=2):
    rho0 = 1 / (lam * RHO_DEN)
    p = rho0 * st['s1']
    M = st['M']
    L = M @ np.log1p(-p)
    if order >= 2:
        w = p / (1 - p)
        Ct = (w[:, None] * w[None, :]) * st['R']
        L = L + 0.5 * np.sum((M @ Ct) * M, axis=1)
    return float(np.sum(rho0 * st['S2'] * np.exp(L)) / st['nclass_all'])


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'check':
        ref = json.load(open(os.path.join(OUT, 'r28.json')))
        res = {}
        for hs in ref:
            st = setup(int(hs))
            for lam in (17.8, 23.0):
                ex = ref[hs]['lam%.1f' % lam]['exact']
                o1 = share(st, lam, 1); o2 = share(st, lam, 2)
                res['%s_%.1f' % (hs, lam)] = {'exact': ex, 'cond_o1_rel': o1 / ex - 1, 'cond_o2_rel': o2 / ex - 1,
                                              'bahadur_o4_rel': ref[hs]['lam%.1f' % lam]['order4_rel'],
                                              'bahadur_o4_overlap_rel': ref[hs]['lam%.1f' % lam]['order4_overlap_rel']}
                print(hs, lam, {k: round(v, 5) for k, v in res['%s_%.1f' % (hs, lam)].items()}, flush=True)
        json.dump(res, open(os.path.join(OUT, 'r47_check_y%d.json' % Y), 'w'), indent=1)
    elif mode == 'u':
        import time
        from scipy.stats import kstest
        LAMS = np.arange(16.0, 26.01, 1.0)
        path = os.path.join(OUT, 'r47_u_y%d.json' % Y)
        res = json.load(open(path)) if os.path.exists(path) else {'lams': LAMS.tolist(), 'ln_share': {}}
        P = {int(g): v for g, v in json.load(open(os.path.join(OUT, 'oes_first_occurrence.json')))['P'].items()}
        t0 = time.time()
        for g in range(180, 301, 2):
            if str(g) in res['ln_share']:
                continue
            if time.time() - t0 > 420:
                break
            st = setup(g)
            res['ln_share'][str(g)] = [math.log(share(st, lam, 2)) for lam in LAMS]
            json.dump(res, open(path, 'w'), indent=1)
            print(g, 'classes', st['M'].shape[0], 'secs', round(time.time() - t0), flush=True)
        done = [g for g in range(180, 301, 2) if str(g) in res['ln_share']]
        if len(done) == 61:
            U = []
            for g in done:
                u = np.linspace(math.log(1e7), math.log(P[g]), 4001)
                ls = np.interp(u, LAMS, res['ln_share'][str(g)])
                U.append(float(np.trapezoid(np.exp(ls) / u * np.exp(u), u)))
            U = np.array(U); ks = kstest(U, 'expon')
            g = np.array(done)
            res['U'] = dict(zip(map(str, done), U.tolist()))
            res['summary'] = {'n': 61, 'mean': float(U.mean()), 'null_se': 1 / math.sqrt(61), 'z_null': float((U.mean() - 1) * math.sqrt(61)),
                              'sample_se': float(U.std(ddof=1) / math.sqrt(61)), 'ks_p': float(ks.pvalue),
                              'mean_six': float(U[g % 6 == 0].mean()), 'mean_not_six': float(U[g % 6 != 0].mean())}
            json.dump(res, open(path, 'w'), indent=1)
            print('summary', res['summary'], flush=True)
        else:
            print('done', len(done), 'of 61', flush=True)
    elif mode == 'o3':
        # order 3: - sum_{i<j<k in R_a} w_i w_j w_k (t_ijk - r_ij - r_ik - r_jk), w = p/(1-p), 1 + t = s_ijk/(s_i s_j s_k)
        out = {}
        for g in map(int, sys.argv[2].split(',')):
            st = setup(g)
            pts = np.arange(2, g, 2, dtype=np.int64); n = len(pts)
            tri = np.array(list(combinations(range(n), 3)), dtype=np.int64)
            i, j, k = tri.T
            z = np.zeros(len(tri), np.int64)
            s3 = S_gt(np.stack([z, np.full(len(tri), g), pts[i], pts[j], pts[k]], 1)) / st['S2']
            s1, R = st['s1'], st['R']
            kap = s3 / (s1[i] * s1[j] * s1[k]) - 1 - R[i, j] - R[i, k] - R[j, k]
            M = st['M'].astype(np.float32)
            res = {}
            for lam in (16.0, 18.0, 20.0):
                rho0 = 1 / (lam * RHO_DEN); p = rho0 * s1; w = (p / (1 - p)).astype(np.float32)
                W = M * w[None, :]
                o3 = np.zeros(M.shape[0])
                for c0 in range(0, len(tri), 4000):
                    sl = slice(c0, c0 + 4000)
                    o3 += (W[:, i[sl]] * W[:, j[sl]] * W[:, k[sl]]) @ kap[sl].astype(np.float32)
                L1 = st['M'] @ np.log1p(-p)
                Ct = (p / (1 - p))[:, None] * (p / (1 - p))[None, :] * R
                L2 = L1 + 0.5 * np.sum((st['M'] @ Ct) * st['M'], axis=1)
                L3 = L2 - o3
                sh2 = float(np.sum(rho0 * st['S2'] * np.exp(L2))); sh3 = float(np.sum(rho0 * st['S2'] * np.exp(L3)))
                sh1 = float(np.sum(rho0 * st['S2'] * np.exp(L1)))
                res[str(lam)] = {'o2_minus_o1_ln': math.log(sh2 / sh1), 'o3_minus_o2_ln': math.log(sh3 / sh2)}
                print(g, lam, {k2: round(v, 5) for k2, v in res[str(lam)].items()}, flush=True)
            out[str(g)] = res
        json.dump(out, open(os.path.join(OUT, 'r47_o3_y%d.json' % Y), 'w'), indent=1)
    elif mode == 'u3':
        # order 3 + the polymer (overlap) correction, as a correction to the order-2 values of mode u; the triple sum
        # per class is sum_i W_i * 1/2 W_{>i} K_i W_{>i}^T with K_i[j,k] = kappa_ijk (symmetric, j,k > i)
        import time
        from scipy.stats import kstest
        L6 = np.array([16.0, 18.0, 20.0, 22.0, 24.0, 26.0])
        path = os.path.join(OUT, 'r47_u3_y%d.json' % Y)
        base = json.load(open(os.path.join(OUT, 'r47_u_y%d.json' % Y)))
        res = json.load(open(path)) if os.path.exists(path) else {'lams': L6.tolist(), 'delta_ln': {}}
        t0 = time.time()
        for g in range(180, 301, 2):
            if str(g) in res['delta_ln']:
                continue
            if time.time() - t0 > 400:
                break
            st = setup(g)
            pts = np.arange(2, g, 2, dtype=np.int64); n = len(pts)
            s1, R, M = st['s1'], st['R'], st['M']
            K = np.zeros((n, n, n))
            for i in range(n - 2):
                jj, kk = np.triu_indices(n - i - 1, 1)
                jj = jj + i + 1; kk = kk + i + 1
                z = np.zeros(len(jj), np.int64)
                s3 = S_gt(np.stack([z, np.full(len(jj), g), np.full(len(jj), pts[i]), pts[jj], pts[kk]], 1)) / st['S2']
                kap = s3 / (s1[i] * s1[jj] * s1[kk]) - 1 - R[i, jj] - R[i, kk] - R[jj, kk]
                K[i, jj, kk] = kap; K[i, kk, jj] = kap
            d = []
            for lam in L6:
                rho0 = 1 / (lam * RHO_DEN); p = rho0 * s1; w = p / (1 - p)
                L1 = M @ np.log1p(-p)
                Cm = w[:, None] * w[None, :] * R
                L2 = L1 + 0.5 * np.sum((M @ Cm) * M, axis=1)
                W = M * w[None, :]
                o3 = np.zeros(M.shape[0])
                for i in range(n - 2):
                    Wi = W[:, i + 1:]
                    o3 += W[:, i] * 0.5 * np.sum((Wi @ K[i, i + 1:, i + 1:]) * Wi, axis=1)
                v = M @ Cm; v2 = M @ (Cm * Cm)
                ov = -0.5 * (0.5 * np.sum(v2 * M, axis=1) + np.sum(M * (v * v - v2), axis=1))
                sh2 = np.sum(np.exp(L2)); sh3 = np.sum(np.exp(L2 - o3 + ov))
                d.append(float(math.log(sh3 / sh2)))
            res['delta_ln'][str(g)] = d
            json.dump(res, open(path, 'w'), indent=1)
            print(g, 'n', n, 'delta', np.round(d, 4), 'secs', round(time.time() - t0), flush=True)
        done = [g for g in range(180, 301, 2) if str(g) in res['delta_ln']]
        if len(done) == 61:
            P = {int(g): v for g, v in json.load(open(os.path.join(OUT, 'oes_first_occurrence.json')))['P'].items()}
            LAMS = np.array(base['lams'])
            U = []
            for g in done:
                u = np.linspace(math.log(1e7), math.log(P[g]), 4001)
                ls = np.interp(u, LAMS, base['ln_share'][str(g)]) + np.interp(u, L6, res['delta_ln'][str(g)])
                U.append(float(np.trapezoid(np.exp(ls) / u * np.exp(u), u)))
            U = np.array(U); ks = kstest(U, 'expon'); garr = np.array(done)
            res['U'] = dict(zip(map(str, done), U.tolist()))
            res['summary'] = {'n': 61, 'mean': float(U.mean()), 'null_se': 1 / math.sqrt(61), 'z_null': float((U.mean() - 1) * math.sqrt(61)),
                              'sample_se': float(U.std(ddof=1) / math.sqrt(61)), 'ks_p': float(ks.pvalue),
                              'mean_six': float(U[garr % 6 == 0].mean()), 'mean_not_six': float(U[garr % 6 != 0].mean())}
            json.dump(res, open(path, 'w'), indent=1)
            print('summary', res['summary'], flush=True)
        else:
            print('done', len(done), 'of 61', flush=True)

