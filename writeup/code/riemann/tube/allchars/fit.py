import json, math
import numpy as np
P = json.load(open('predict.json')); S = json.load(open('scan_100000000.json'))
names = list(P)
w1 = np.array([P[n]['w1_over_C'] for n in names]); fr = np.array([P[n]['obs_frac'] for n in names])
g1 = np.array([abs(P[n]['gamma1']) for n in names])
# zeta: gamma1 = 14.134725, C = 2 + gamma - ln 4pi
Cz = 2 + 0.5772156649015329 - math.log(4 * math.pi); wz = 2 / (0.25 + 14.134725141734693 ** 2) / Cz
print('zeta w1/C = %.3f frac = %.3f kappa = %.3f' % (wz, S['zeta']['frac'], (S['zeta']['frac'] - wz) / (1 - wz)))
k = (fr - w1) / (1 - w1)
for n, a, b, c in zip(names, w1, fr, k):
    print('%-14s %s %-7s w1/C=%.3f frac=%.3f kappa=%.3f' % (n, 'odd ' if S[n]['kind'] == 'odd' else 'even', 'real' if n in ('chi3','chi4','chi8','chi-8','chi12') or 'ord2' in n else 'complex', a, b, c))
print('kappa mean %.3f sd %.3f min %.3f max %.3f' % (k.mean(), k.std(), k.min(), k.max()))
A = np.vstack([np.ones_like(w1), w1]).T; coef, res, *_ = np.linalg.lstsq(A, fr, rcond=None)
pred = A @ coef; print('fit frac = %.3f + %.3f * w1/C, rms resid %.4f, corr %.3f' % (coef[0], coef[1], np.sqrt(np.mean((fr - pred) ** 2)), np.corrcoef(w1, fr)[0, 1]))
print('corr(frac, gamma1) %.3f, corr(frac, 1/(1/4+g1^2)) %.3f' % (np.corrcoef(g1, fr)[0, 1], np.corrcoef(1 / (0.25 + g1 ** 2), fr)[0, 1]))
q = np.array([S[n]['q'] for n in names]); odd = np.array([S[n]['kind'] == 'odd' for n in names])
print('corr(frac, q) %.3f; mean frac odd %.3f even %.3f; mean kappa odd %.3f even %.3f' % (np.corrcoef(q, fr)[0, 1], fr[odd].mean(), fr[~odd].mean(), k[odd].mean(), k[~odd].mean()))
pr = np.array([P[n]['pred_frac'] for n in names]); print('zero-sum (|g|<30) prediction vs observed: max |diff| %.4f (in units of C), mean %.4f' % (np.abs(pr - fr).max(), (pr - fr).mean()))
Z = json.load(open('zeros_30.json'))
W2 = []
for n in names:
    g = np.array(sorted(Z[n]['zeros'], key=abs)); w = 1 / (0.25 + g ** 2)
    W2.append(w[:2].sum() / S[n]['C'])
W2 = np.array(W2); k2 = (fr - W2) / (1 - W2)
A = np.vstack([np.ones_like(W2), W2]).T; coef = np.linalg.lstsq(A, fr, rcond=None)[0]
print('two lowest zeros (pair for real chi): fit frac = %.3f + %.3f*W2/C, rms %.4f, corr %.3f; kappa2 mean %.3f sd %.3f; odd %.3f even %.3f; real %.3f complex %.3f' % (
    coef[0], coef[1], np.sqrt(np.mean((fr - A @ coef) ** 2)), np.corrcoef(W2, fr)[0, 1], k2.mean(), k2.std(), k2[odd].mean(), k2[~odd].mean(),
    k2[[Z[n]['real'] for n in names]].mean(), k2[[not Z[n]['real'] for n in names]].mean()))
for n, a, b in zip(names, W2, k2): print('  %-14s W2/C=%.3f kappa2=%.3f' % (n, a, b))
