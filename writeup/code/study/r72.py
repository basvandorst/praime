"""Round 72 (preregister/round72.txt): an online predictor of the next prime gap from the previous gaps only.
The past primes p - D_j are coprime to every q <= 13, so p mod 30030 lies in A_k = {a coprime : a - D_j coprime, j <= k};
the predictor averages the conditioned law P(g | a) (round 47, order 2, y = 13) over A_k. Its learning curve on real
primes is compared with the coupon-collector shape of round 61 restricted to q = 3, 5, 7, 11, 13.
usage: r72.py table | run"""
import os, sys, json, math, time
import numpy as np
import common as C
import r47
import r61

OUT = r47.OUT
Q = r47.PP
CLS = np.array([a for a in range(Q) if math.gcd(a, Q) == 1], dtype=np.int64)
IDX = np.full(Q, -1); IDX[CLS] = np.arange(len(CLS))
COP = np.zeros(Q, bool); COP[CLS] = True
LAM = math.log(1.05e9)
G = int(30 * LAM) // 2 * 2
TAB = os.path.join(os.environ.get('STUDY_SCRATCH', '/tmp'), 'r72_table.npy')


def per_class(g):
    """P(gap = g | p = a mod Q) for every coprime class a (order 2 of the conditioned expansion)."""
    st = r47.setup(g)
    rho0 = 1 / (LAM * r47.RHO_DEN); p = rho0 * st['s1']; M = st['M']
    L = M @ np.log1p(-p) if M.shape[1] else np.zeros(M.shape[0])
    if M.shape[1] > 1:
        w = p / (1 - p)
        L = L + 0.5 * np.sum((M @ ((w[:, None] * w[None, :]) * st['R'])) * M, axis=1)
    cls = CLS[np.gcd(CLS + g, Q) == 1]
    out = np.zeros(len(CLS)); out[IDX[cls]] = rho0 * st['S2'] * np.exp(L)
    return out


if __name__ == '__main__':
    if sys.argv[1] == 'table':
        T = np.load(TAB) if os.path.exists(TAB) else np.full((len(CLS), G // 2), np.nan)
        t0 = time.time()
        for j, g in enumerate(range(2, G + 1, 2)):
            if not np.isnan(T[0, j]) or time.time() - t0 > 450:
                continue
            T[:, j] = per_class(g)
            np.save(TAB, T)
        print('table columns done', int(np.sum(~np.isnan(T[0]))), 'of', G // 2, 'row sums', np.nanmin(np.nansum(T, 1)), np.nanmax(np.nansum(T, 1)), flush=True)
    else:
        T = np.load(TAB); T = T / T.sum(1, keepdims=True)          # renormalise each class over g <= G
        mix = T.mean(0)
        H = lambda p: float(-np.sum(p[p > 0] * np.log2(p[p > 0])))
        I_model = H(mix) - float(np.mean([H(r) for r in T]))
        P = np.nonzero(C.primes_mask(10 ** 9, 10 ** 8))[0].astype(np.int64) + 10 ** 9
        rng = np.random.default_rng(72); K = 60
        tg = rng.choice(np.arange(K + 1, len(P) - 1), 20000, replace=False)
        L = np.zeros(K + 1); Ltrue = 0.0; n = 0
        for b0 in range(0, len(tg), 1000):
            t = tg[b0:b0 + 1000]
            gnext = P[t + 1] - P[t]; ok = gnext <= G; t = t[ok]; gnext = gnext[ok]
            col = T[:, gnext // 2 - 1].T                              # (targets, classes)
            mask = np.ones((len(t), len(CLS)), bool)
            for k in range(K + 1):
                if k > 0:
                    D = P[t] - P[t - k]
                    mask &= COP[(CLS[None, :] - D[:, None]) % Q]
                pk = (col * mask).sum(1) / mask.sum(1)
                L[k] += float(np.sum(-np.log2(pk)))
            Ltrue += float(np.sum(-np.log2(col[np.arange(len(t)), IDX[P[t] % Q]])))
            n += len(t)
        L /= n; Ltrue /= n
        gain = L[0] - Ltrue
        shape_data = [(L[0] - L[k]) / gain for k in range(K + 1)]
        ks = list(range(1, K + 1)); qs = (3, 5, 7, 11, 13)
        full = sum(((q - 2) / (q - 1)) * math.log((q - 1) / (q - 2)) for q in qs)
        Ik = {k: sum(r61.I_k_q(q, [k])[k] for q in qs) for k in ks}
        shape_th = [0.0] + [Ik[k] / full for k in ks]
        dev = max(abs(a - b) for a, b in zip(shape_data[1:], shape_th[1:]))
        out = {'lambda': LAM, 'G': G, 'targets': n, 'L_k': L.tolist(), 'L_true': Ltrue, 'gain': gain, 'I_model': I_model,
               'gain_minus_I_model': gain - I_model, 'shape_data': shape_data, 'shape_coupon': shape_th, 'max_shape_dev': dev}
        for k in (0, 1, 2, 3, 5, 10, 20, 40, 60):
            print('k %2d  L %.4f  shape data %.3f  coupon %.3f' % (k, L[k], shape_data[k], shape_th[k]))
        print('L_true %.4f gain %.4f I_model %.4f diff %+.4f max shape dev %.3f' % (Ltrue, gain, I_model, gain - I_model, dev), flush=True)
        json.dump(out, open(os.path.join(OUT, 'r72.json'), 'w'), indent=1)
