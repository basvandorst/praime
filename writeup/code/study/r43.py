"""Round 43 (preregister/round43.txt): the spread of twin counts across the 1485 admissible classes mod 30030 (and
mod 210, 2310 alongside), against the Hardy-Littlewood four-point sums, with round 42's code. usage: r43.py"""
import os, json, math
import numpy as np
from scipy.stats import chi2
import common as C
import theory_r3 as R3
import r42

WINS = {'train': (10 ** 7, 9 * 10 ** 7), '1e8-1e9': (10 ** 8, 9 * 10 ** 8), '1e9-1e10': (10 ** 9, 9 * 10 ** 9),
        'A': (10 ** 10, 5 * 10 ** 7), 'B': (10 ** 12, 5 * 10 ** 7), 'C': (10 ** 15, 2 * 10 ** 7)}
QS = (210, 2310, 30030)
PRIMARY = ('train', '1e8-1e9', '1e9-1e10')

if __name__ == '__main__':
    res = {}
    for win, (L0, W) in WINS.items():
        cnt = r42.count(L0, W, QS)
        res[win] = {}
        for Q in QS:
            if Q == 210 and W > 10 ** 9:
                continue                                   # 4.3e7 points per class: beyond the memory cap (not run)
            adm = [a for a in range(Q) if math.gcd(a * (a + 2), Q) == 1]
            X = cnt[Q][adm]; K = len(adm)
            meas = float(np.var(X, ddof=1) / X.mean()); pred = r42.predict(L0, W, Q)
            res[win][str(Q)] = {'classes': K, 'twins': int(X.sum()), 'measured': meas, 'predicted': pred,
                                'ratio': meas / pred, 'se_ratio': math.sqrt(2 / (K - 1)),
                                'p_lower': float(chi2.cdf((K - 1) * meas / pred, K - 1))}
            print(win, Q, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in res[win][str(Q)].items()}, flush=True)
    for name, wins in (('primary_30030', PRIMARY), ('all_30030', tuple(WINS))):
        rs = np.array([res[w]['30030']['ratio'] for w in wins]); se = np.array([res[w]['30030']['se_ratio'] for w in wins])
        wt = 1 / se ** 2
        pr = float(np.sum(wt * rs) / wt.sum()); ps = float(1 / math.sqrt(wt.sum()))
        res[name] = {'windows': list(wins), 'ratio': pr, 'se': ps, 'z': (pr - 1) / ps}
        print(name, res[name], flush=True)
    json.dump(res, open(os.path.join(R3.OUT, 'r43.json'), 'w'), indent=1)
