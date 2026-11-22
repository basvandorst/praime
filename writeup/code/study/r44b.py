"""Round 44 addendum (after the tenth review): the spread of twin counts per class over [1e9, 1e10) split into the part
within the nine windows of 1e9 and the part across them (class deviations that do or do not cancel between windows).
In units of the decade's mean count per class m: spread = SS_dec/((K-1) m) = within + cross, within = sum_k SS_k/((K-1) m).
Predicted within = sum_k v_k m_k / m; predicted total = v_dec (round 43 for mod 2310; review9's chunked value for mod 210,
the one-piece FFT being beyond the memory cap at 4.3e7 points per class). usage: r44b.py"""
import os, json, math
import numpy as np
import theory_r3 as R3
import r42

V_DEC = {'2310': json.load(open(os.path.join(R3.OUT, 'r43.json')))['1e9-1e10']['2310']['predicted'], '210': 0.4492}
if __name__ == '__main__':
    out = {}
    cnts = [r42.count(k * 10 ** 9, 10 ** 9, (210, 2310)) for k in range(1, 10)]
    for Q in (210, 2310):
        adm = [a for a in range(Q) if math.gcd(a * (a + 2), Q) == 1]; K = len(adm)
        Xk = np.array([c[Q][adm] for c in cnts])                # 9 x K
        X = Xk.sum(0); m = X.mean()
        ss_dec = float(np.sum((X - m) ** 2)); ss_k = [float(np.sum((x - x.mean()) ** 2)) for x in Xk]
        total = ss_dec / ((K - 1) * m); within = sum(ss_k) / ((K - 1) * m)
        vk = [r42.predict(k * 10 ** 9, 10 ** 9, Q) for k in range(1, 10)]
        within_pred = sum(v * x.mean() for v, x in zip(vk, Xk)) / m
        out[str(Q)] = {'total': total, 'within': within, 'cross': total - within, 'total_pred': V_DEC[str(Q)],
                       'within_pred': within_pred, 'cross_pred': V_DEC[str(Q)] - within_pred}
        print(Q, {k: round(v, 3) for k, v in out[str(Q)].items()}, flush=True)
    json.dump(out, open(os.path.join(R3.OUT, 'r44b.json'), 'w'), indent=1)
