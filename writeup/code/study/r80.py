"""Round 80 (preregister/round80.txt): a second-order chain from round 79's two-gap law. Backwards from p, the gap h_k is
drawn from P_c(h_{k-1}, h_k) / sum, c the reflected class of the prime two steps back (h_1 from the one-gap law);
mean reveal time of p mod 30030, and the repeat probability three steps apart against the primes in [1e9, 1.1e9).
The first-order chain (one-gap law only) is run alongside with the same code. usage: r80.py"""
import os, json, math, time
import numpy as np
import common as C
import r47

assert r47.Y == 13
OUT = r47.OUT
PP = r47.PP
CLS = np.array([a for a in range(PP) if math.gcd(a, PP) == 1], dtype=np.int64)
IDX = np.full(PP, -1); IDX[CLS] = np.arange(len(CLS))
LAM = math.log(1.05e9); RHO0 = 1 / (LAM * r47.RHO_DEN)
G = 300; NG = G // 2
QS = (3, 5, 7, 11, 13)
KMAX = 150; NCH = 200_000


def tables():
    """P1[c, g] (one gap) and J3[c, g1, g2] (two consecutive gaps), forward from a prime in class c, spans <= G."""
    P1 = np.zeros((len(CLS), NG)); J3 = np.zeros((len(CLS), NG, NG), np.float32)
    for g in range(2, G + 1, 2):
        st = r47.setup(g); pts = np.arange(2, g, 2, dtype=np.int64); n = len(pts)
        rows = IDX[CLS[np.gcd(CLS + g, PP) == 1]]
        M = st['M']; p = RHO0 * st['s1']
        L = M @ np.log1p(-p) if n else np.zeros(len(rows))
        if n > 1:
            w = p / (1 - p); V = M @ ((w[:, None] * w[None, :]) * st['R']); Qd = 0.5 * np.sum(V * M, axis=1)
        else:
            V = np.zeros((len(rows), n)); Qd = np.zeros(len(rows))
        P1[rows, g // 2 - 1] = RHO0 * st['S2'] * np.exp(L + Qd)
        for k in range(n):
            g1 = int(pts[k]); ok = M[:, k] == 1
            S3 = float(r47.S_gt(np.array([[0, g1, g]]))[0])
            oth = np.delete(np.arange(n), k)
            if len(oth):
                s1p = r47.S_gt(np.stack([np.zeros(len(oth), np.int64), np.full(len(oth), g1), np.full(len(oth), g), pts[oth]], 1)) / S3
                Lp = M[ok][:, oth] @ np.log1p(-RHO0 * s1p)
            else:
                Lp = np.zeros(int(ok.sum()))
            J3[rows[ok], g1 // 2 - 1, (g - g1) // 2 - 1] = RHO0 ** 2 * S3 * np.exp(Lp + Qd[ok] - V[ok, k])
    return P1, J3


def draw(W, rng):
    cum = np.cumsum(W, axis=1); u = rng.random(len(W)) * cum[:, -1]
    return 2 * (np.argmax(cum >= u[:, None], axis=1) + 1)


def run(P1, J3, order, rng):
    Tall = []; Tq = {q: [] for q in QS}; lag3 = {q: 0 for q in QS}; lag3n = 0; lag2 = {q: 0 for q in QS}; lag2n = 0
    for b0 in range(0, NCH, 10_000):
        n = 10_000; a = CLS[rng.integers(0, len(CLS), n)]
        D = np.zeros((KMAX + 1, n), np.int64); h = np.zeros((KMAX + 1, n), np.int64)
        seen = {q: np.zeros((n, q), bool) for q in QS}; fq = {q: np.full(n, -1) for q in QS}
        for q in QS:
            seen[q][:, 0] = True
        for k in range(1, KMAX + 1):
            if order == 1 or k == 1:
                c = IDX[(-(a - D[k - 1])) % PP]                  # reflected class of the current prime
                g = draw(P1[c], rng)
            else:
                c = IDX[(-(a - D[k - 2])) % PP]                  # reflected class of the prime two steps back
                g = draw(J3[c, h[k - 1] // 2 - 1, :].astype(np.float64), rng)
            h[k] = g; D[k] = D[k - 1] + g
            for q in QS:
                seen[q][np.arange(n), D[k] % q] = True
                fq[q][(seen[q].sum(1) == q - 1) & (fq[q] < 0)] = k
        st = np.stack([fq[q] for q in QS]); t = np.max(st, 0); t[np.any(st < 0, 0)] = KMAX + 1
        Tall.append(t)
        for q in QS:
            Tq[q].append(np.where(fq[q] < 0, KMAX + 1, fq[q]))
            lag3[q] += int(np.sum((D[3:] - D[:-3]) % q == 0)); lag2[q] += int(np.sum((D[2:] - D[:-2]) % q == 0))
        lag3n += (KMAX - 2) * n; lag2n += (KMAX - 1) * n
    Tall = np.concatenate(Tall)
    return {'mean_T': float(np.mean(Tall)), 'sem': float(np.std(Tall) / math.sqrt(len(Tall))), 'unresolved': int(np.sum(Tall > KMAX)),
            **{'mean_T_q%d' % q: float(np.mean(np.concatenate(Tq[q]))) for q in QS},
            'lag3_repeat': {str(q): lag3[q] / lag3n for q in QS}, 'lag2_repeat': {str(q): lag2[q] / lag2n for q in QS}}


if __name__ == '__main__':
    t0 = time.time()
    P1, J3 = tables()
    print('tables', round(time.time() - t0), 's', flush=True)
    res = {'lambda': LAM, 'G': G, 'chains': NCH}
    res['order2'] = run(P1, J3, 2, np.random.default_rng(80))
    res['order1'] = run(P1, J3, 1, np.random.default_rng(81))
    P = np.nonzero(C.primes_mask(10 ** 9, 10 ** 8))[0].astype(np.int64) + 10 ** 9
    res['data'] = {}
    for q in QS:
        m = float(np.mean((P[3:] % q) == (P[:-3] % q))); se = math.sqrt(m * (1 - m) / (len(P) - 3))
        res['data'][str(q)] = {'lag3': m, 'se': se, 'z_order2': (res['order2']['lag3_repeat'][str(q)] - m) / se,
                               'z_order1': (res['order1']['lag3_repeat'][str(q)] - m) / se}
    res['data_mean_T'] = 33.9922
    print(json.dumps({k: v for k, v in res.items()}, indent=0)[:3000], flush=True)
    json.dump(res, open(os.path.join(OUT, 'r80.json'), 'w'), indent=1)
