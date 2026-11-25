"""expected errors of the registered readings from the measured noise (no zero used). usage: c9_power.py predictions_blind4.json"""
import sys, json, math
import numpy as np
P = json.load(open(sys.argv[1])); L = P['L']; S = P['sets']; E = P['expected_se']
for k in ('all', 'sym3', 'sym3_eps-1', 'sym3_eps+1', 'sym3_S70_eps-1', 'sym3_S70_eps+1', 'sym3_S50_eps-1', 'sym3_S50_eps+1', 'g2', 'rs'):
    print('%-16s n=%3d expected se %.4f  3-se band 1 +- %.3f' % (k, len(S[k]), E[k], 3 * E[k]))
sm, sp = E['sym3_eps-1'], E['sym3_eps+1']
print('R2 watch item: if b(-1) = 0.11 (round 151), it lies %.1f se from 1 -> EFFECT; if b(-1) = 1, b/se = %.1f -> FLUCTUATION; difference se %.3f' % ((1 - 0.11) / sm, 1 / sm, math.sqrt(sm ** 2 + sp ** 2)))
X = []; W = []
for k, e in L.items():
    unit = e['T2'] / (2 * math.pi); X.append([1, e['c'] - 1]); W.append((unit / e['sd']['along']) ** 2)
X = np.array(X); W = np.array(W); cov = np.linalg.inv(X.T @ (W[:, None] * X))
print('R5 amplitude: se(alpha) %.4f, se(beta) %.4f (H-c1 at %.0f se)' % (math.sqrt(cov[0, 0]), math.sqrt(cov[1, 1]), 1 / math.sqrt(cov[1, 1])))
num = sum(abs(complex(*v)) ** 2 for e in L.values() for v in e['lines'].values())
var = sum(abs(complex(*v)) ** 2 * e['sd']['line_along'] ** 2 for e in L.values() for v in e['lines'].values())
print('R6 lines: se %.4f' % (math.sqrt(var) / num))
print('R7 Hann all: se %.4f ; R8 band Sym3 eps=-1: se %.4f' % (E['hann_all'], E['band_sym3_eps-1']))
print('falsifier band for R1: 1 +- max(4 x %.4f, 0.10) = 1 +- %.3f' % (E['all'], max(4 * E['all'], 0.10)))
print('estimated zero cost %.0f s = %.2f CPU-hours' % (P['est_total_secs'], P['est_total_secs'] / 3600))
