"""Round 79 (preregister/round79.txt): the repeat probability two steps apart, P(p_{n+2} = p_n mod q), from the
conditioned law for two consecutive gaps at once, against the chain of one-gap laws and the data near 1e9.
P_a(g1, g2) = rho0^2 S_{>13}({0, g1, g}) exp(L'_a + Q'_a), g = g1 + g2: L' first order over the rough interior points
other than g1 with s1 relative to {0, g1, g} (exact); Q' the order-2 term of the single-gap base {0, g} minus the
terms of g1 (approximation). usage: r79.py"""
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

if __name__ == '__main__':
    t0 = time.time()
    P1 = np.zeros((len(CLS), NG)); J = np.zeros((len(CLS), NG)); K = np.zeros((len(CLS), NG))
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
            val = RHO0 ** 2 * S3 * np.exp(Lp + Qd[ok] - V[ok, k])
            J[rows[ok], g // 2 - 1] += val; K[rows[ok], g1 // 2 - 1] += val
        if g % 50 == 0:
            print('g', g, round(time.time() - t0), 's', flush=True)
    # chain of one-gap laws, spans <= G
    Jc = np.zeros_like(J)
    for j in range(NG):
        g1 = 2 * (j + 1); nxt = IDX[(CLS + g1) % PP]; okc = nxt >= 0
        m = NG - j - 1
        if m <= 0:
            continue
        Jc[okc, j + 1:] += P1[okc, j:j + 1] * P1[nxt[okc], :m]
    spans = np.arange(2, G + 1, 2)
    res = {'lambda': LAM, 'G': G,
           'one_gap_row_sums': [float(P1.sum(1).min()), float(P1.sum(1).max())],
           'two_gap_row_sums': [float(J.sum(1).min()), float(J.sum(1).max())],
           'first_gap_marginal_vs_one_gap_max_abs': float(np.max(np.abs(K / K.sum(1, keepdims=True) - P1 / P1.sum(1, keepdims=True)))),
           'data': {'3': 0.4911, '5': 0.2367, '7': 0.1485, '11': 0.0765, '13': 0.0598},
           'data_se': {'3': 0.00023, '5': 0.00019, '7': 0.00016, '11': 0.00012, '13': 0.00011}}
    for name, A in (('two_gap_law', J), ('chain', Jc)):
        res[name] = {str(q): float(A[:, spans % q == 0].sum() / A.sum()) for q in QS}
    res['one_step_law'] = {str(q): float(P1[:, spans % q == 0].sum() / P1.sum()) for q in QS}
    res['z_two_gap'] = {k: (res['two_gap_law'][k] - res['data'][k]) / res['data_se'][k] for k in res['data']}
    res['z_chain'] = {k: (res['chain'][k] - res['data'][k]) / res['data_se'][k] for k in res['data']}
    P = np.nonzero(C.primes_mask(10 ** 9, 10 ** 8))[0].astype(np.int64) + 10 ** 9      # data: primes in [1e9, 1.1e9)
    res['data_pairs'] = int(len(P) - 2); res['data_exact'] = {}
    for q in QS:
        m = float(np.mean((P[2:] % q) == (P[:-2] % q))); se = math.sqrt(m * (1 - m) / (len(P) - 2))
        res['data_exact'][str(q)] = {'p': m, 'se': se, 'z_two_gap': (res['two_gap_law'][str(q)] - m) / se, 'z_chain': (res['chain'][str(q)] - m) / se}
    print(json.dumps(res['data_exact'], indent=0), res['two_gap_row_sums'], res['first_gap_marginal_vs_one_gap_max_abs'], flush=True)
    json.dump(res, open(os.path.join(OUT, 'r79.json'), 'w'), indent=1)
