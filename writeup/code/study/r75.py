"""Round 75 (preregister/round75.txt): round 74's chain with the empirical backward-gap table of the primes in
[1e9, 1.1e9) (gap before p, by p mod 30030) instead of the theory table: keeps the true class dependence, drops the
dependence of a gap on earlier gaps beyond the class. usage: r75.py"""
import os, json
import numpy as np
import common as C
from r72 import CLS, IDX, Q

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
KMAX = 150
NCH = 200_000
QS = (3, 5, 7, 11, 13)


def chain(cum, gaps, step_row, rng):
    """Mean reveal times; step_row(cur) gives the table row of the backward law from class cur."""
    Tall = []; Tq = {q: [] for q in QS}
    for b0 in range(0, NCH, 10_000):
        n = min(10_000, NCH - b0)
        cur = CLS[rng.integers(0, len(CLS), n)]; D = np.zeros(n, np.int64)
        seen = {q: np.zeros((n, q), bool) for q in QS}; fq = {q: np.full(n, -1) for q in QS}
        for q in QS:
            seen[q][:, 0] = True
        for k in range(1, KMAX + 1):
            row = step_row(cur); u = rng.random(n)
            j = np.argmax(cum[row] >= u[:, None], axis=1)
            g = gaps[j]; D += g; cur = (cur - g) % Q
            for q in QS:
                seen[q][np.arange(n), D % q] = True
                fq[q][(seen[q].sum(1) == q - 1) & (fq[q] < 0)] = k
        st = np.stack([fq[q] for q in QS]); t = np.max(st, 0); t[np.any(st < 0, 0)] = KMAX + 1   # unresolved by KMAX
        Tall.append(t)
        for q in QS:
            Tq[q].append(np.where(fq[q] < 0, KMAX + 1, fq[q]))
    Tall = np.concatenate(Tall)
    return {'unresolved': int(np.sum(Tall > KMAX)), 'mean_T': float(np.mean(Tall)), **{'mean_T_q%d' % q: float(np.mean(np.concatenate(Tq[q]))) for q in QS},
            'cdf': [float(np.mean(Tall <= k)) for k in range(KMAX + 1)]}


if __name__ == '__main__':
    P = np.nonzero(C.primes_mask(10 ** 9, 10 ** 8))[0].astype(np.int64) + 10 ** 9
    b = np.diff(P); cls = IDX[P[1:] % Q]
    G = int(b.max()); gaps = np.arange(2, G + 1, 2)
    E = np.zeros((len(CLS), len(gaps)))
    np.add.at(E, (cls, b // 2 - 1), 1)
    print('primes', len(P), 'min per class', int(E.sum(1).min()), flush=True)
    cum = np.cumsum(E / E.sum(1, keepdims=True), 1); cum[:, -1] = 1.0
    r = chain(cum, gaps, lambda cur: IDX[cur], np.random.default_rng(75))
    r74 = json.load(open(os.path.join(OUT, 'r74.json')))
    res = {'empirical_table_chain': {k: v for k, v in r.items()}, 'theory_chain_mean': r74['mean_T_chain'], 'data_mean': r74['mean_T_data'],
           'P_gap0_empirical_table': {str(q): float(np.sum(E[:, (gaps % q) == 0]) / E.sum()) for q in QS}}
    cd = np.array(json.load(open(os.path.join(OUT, 'r73.json')))['cdf_data'])
    res['max_cdf_diff_vs_data'] = float(np.max(np.abs(np.array(r['cdf']) - cd)))
    print({k: v for k, v in r.items() if k != 'cdf'}, 'theory', r74['mean_T_chain'], 'data', r74['mean_T_data'], 'cdf diff', res['max_cdf_diff_vs_data'], flush=True)
    json.dump(res, open(os.path.join(OUT, 'r75.json'), 'w'), indent=1)
