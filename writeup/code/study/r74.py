"""Round 74 (preregister/round74.txt): the time T for the past gaps to reveal p mod 30030, predicted by a Markov chain on
the classes mod 30030 whose backward step from class a is a gap drawn from the conditioned law P(g | -a) of round 72.
By the CRT the class mod 30030 is fixed exactly when it is fixed mod every q, so T = max_q T_q.
Needs round 72's table (r72.py table, in STUDY_SCRATCH). usage: r74.py"""
import os, json
import numpy as np
from r72 import CLS, IDX, COP, Q, TAB
from r73 import predicted_cdf

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
KMAX = 150
NCH = 200_000
QS = (3, 5, 7, 11, 13)

if __name__ == '__main__':
    T = np.load(TAB); T = T / T.sum(1, keepdims=True)
    cum = np.cumsum(T, 1); cum[:, -1] = 1.0
    NB = 4096                                                     # guide table: first index with cum >= b/NB
    guide = np.stack([np.searchsorted(c, np.arange(NB) / NB) for c in cum]).astype(np.int16)
    gaps = np.arange(2, 2 * T.shape[1] + 1, 2)
    rng = np.random.default_rng(74)
    data = json.load(open(os.path.join(OUT, 'r73.json')))
    res = {'chains': NCH}
    Tall = np.full(NCH, -1); Tq = {q: np.full(NCH, -1) for q in QS}; rep = {q: 0 for q in QS}; nsteps = 0
    for b0 in range(0, NCH, 10_000):
        n = min(10_000, NCH - b0)
        a0 = CLS[rng.integers(0, len(CLS), n)]                    # p mod Q
        cur = a0.copy(); D = np.zeros(n, np.int64)
        seen = {q: np.zeros((n, q), bool) for q in QS}; fq = {q: np.full(n, -1) for q in QS}
        for q in QS:
            seen[q][:, 0] = True
        for k in range(1, KMAX + 1):
            row = IDX[(-cur) % Q]                                  # backward law = forward law of the reflected class
            u = rng.random(n)
            j = guide[row, (u * NB).astype(np.int64)].astype(np.int64)
            while True:
                low = cum[row, j] < u
                if not low.any():
                    break
                j[low] += 1
            g = gaps[j]
            for q in QS:
                rep[q] += int(np.sum(g % q == 0))
            nsteps += n
            D += g; cur = (cur - g) % Q
            for q in QS:
                seen[q][np.arange(n), D % q] = True
                fq[q][(seen[q].sum(1) == q - 1) & (fq[q] < 0)] = k
        assert np.all(np.gcd(cur, Q) == 1)                          # the chain stays in coprime classes
        Tall[b0:b0 + n] = np.max(np.stack([fq[q] for q in QS]), 0) if all(np.all(fq[q] >= 0) for q in QS) else -1
        for q in QS:
            Tq[q][b0:b0 + n] = fq[q]
    cdf = np.array([np.mean((Tall >= 0) & (Tall <= k)) for k in range(KMAX + 1)])
    cd = np.array(data['cdf_data'])
    res['max_cdf_diff_vs_data'] = float(np.max(np.abs(cdf - cd))); res['at_k'] = int(np.argmax(np.abs(cdf - cd)))
    res['mean_T_chain'] = float(np.mean(Tall[Tall >= 0])); res['unresolved'] = int(np.sum(Tall < 0))
    res['mean_T_data'] = data['mean_T_data']; res['mean_T_coupon'] = data['mean_T_predicted']
    for q in QS:
        res['mean_T_q%d_chain' % q] = float(np.mean(Tq[q])); res['mean_T_q%d_data' % q] = data['diag']['mean_T_q%d' % q]
        res['P_gap0_q%d_chain' % q] = rep[q] / nsteps
    res['cdf_chain'] = cdf.tolist()
    res['max_cdf_diff_coupon_vs_data'] = float(np.max(np.abs(np.array(predicted_cdf(KMAX)) - cd)))
    print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in res.items() if k != 'cdf_chain'})
    json.dump(res, open(os.path.join(OUT, 'r74.json'), 'w'), indent=1)
