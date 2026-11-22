"""Round 44 (preregister/round44.txt): the spread of twin counts across the classes mod 210 and 2310 in the nine
disjoint windows [k 1e9, (k+1) 1e9), k = 1..9, against the four-point prediction of round 42. usage: r44.py"""
import os, json, math
import numpy as np
from scipy.stats import chi2
import theory_r3 as R3
import r42

if __name__ == '__main__':
    res = {'windows': {}}
    for k in range(1, 10):
        L0, W = k * 10 ** 9, 10 ** 9
        cnt = r42.count(L0, W, (210, 2310))
        d = {}
        for Q in (210, 2310):
            adm = [a for a in range(Q) if math.gcd(a * (a + 2), Q) == 1]
            X = cnt[Q][adm]; K = len(adm)
            meas = float(np.var(X, ddof=1) / X.mean()); pred = r42.predict(L0, W, Q)
            d[str(Q)] = {'twins': int(X.sum()), 'measured': meas, 'predicted': pred, 'ratio': meas / pred,
                         'p_lower': float(chi2.cdf((K - 1) * meas / pred, K - 1))}
        # level 2: spread of the 11 subclasses mod 2310 within each class mod 210, against the prediction
        # (the within part of the mod-2310 sum of squares: SS2310 - SS210/11 on 135 - 15 = 120 degrees of freedom)
        adm2 = [a for a in range(2310) if math.gcd(a * (a + 2), 2310) == 1]
        X2 = cnt[2310][adm2]; m2 = X2.mean()
        ss2 = float(np.sum((X2 - m2) ** 2))
        adm1 = [a for a in range(210) if math.gcd(a * (a + 2), 210) == 1]
        X1 = cnt[210][adm1]; ss1 = float(np.sum((X1 - X1.mean()) ** 2)) / 9      # each class mod 210 has 9 subclasses
        v2, v1 = d['2310']['predicted'], d['210']['predicted']
        # expected: E SS2310 = 134 v2 m2 ; E SS210/9 = 14 v1 m1 / 9 with m1 = 9 m2  -> 14 v1 m2
        within_meas = (ss2 - ss1) / m2; within_pred = 134 * v2 - 14 * v1
        d['level2_within'] = {'ratio': within_meas / within_pred, 'dof': 120}
        res['windows'][str(k)] = d
        print(k, {q: round(d[q]['ratio'], 3) for q in ('210', '2310')}, 'level2 %.3f' % d['level2_within']['ratio'], flush=True)
    r = np.array([res['windows'][str(k)]['2310']['ratio'] for k in range(1, 10)]); se = math.sqrt(2 / 134)
    res['pooled_2310'] = {'ratio': float(r.mean()), 'se': se / 3, 'z': float((r.mean() - 1) / (se / 3)),
                          'scatter_sd': float(r.std(ddof=1)), 'nominal_sd': se}
    r1 = np.array([res['windows'][str(k)]['210']['ratio'] for k in range(1, 10)])
    res['pooled_210'] = {'ratio': float(r1.mean()), 'se': math.sqrt(2 / 14) / 3}
    r3 = np.array([res['windows'][str(k)]['level2_within']['ratio'] for k in range(1, 10)])
    res['pooled_level2'] = {'ratio': float(r3.mean()), 'se': math.sqrt(2 / 120) / 3}
    print({k: v for k, v in res.items() if k != 'windows'}, flush=True)
    json.dump(res, open(os.path.join(R3.OUT, 'r44.json'), 'w'), indent=1)
