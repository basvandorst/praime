"""sizes.py: the size of each part of the recipe at the heights of the test (development triples)."""
import os, json, math, itertools
import numpy as np
G = 0.57721566490153286061
rows = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dev_pred_triples.json')))
monos = list(itertools.product((0, 1), repeat=3))
for t in [100.0, 316.0, 1000.0]:
    fr = {0: [], 1: [], 2: [], 3: []}; dl = []; ld = []
    for r in rows:
        q = [abs(D) for D in r['Ds']]; L = [math.log(x * t / (2 * math.pi)) for x in q]
        terms = []
        for i, m in enumerate(monos):
            v = r['hm'][i]
            for j in range(3):
                if m[j] == 0: v *= (L[j] + 2 * G)
            terms.append((sum(m), v))
        P = sum(v for _, v in terms)
        for o in range(4): fr[o].append(sum(v for s, v in terms if s == o) / P)
        lead = r['hm'][0] * L[0] * L[1] * L[2]; ld.append(lead / P)
    print('t = %5.0f  P split by number of differentiated shifts |m| (median [min, max] over 56 dev triples):' % t)
    for o in range(4):
        a = np.array(fr[o]); print('     |m| = %d: %+.3f [%+.3f, %+.3f]' % (o, np.median(a), a.min(), a.max()))
    a = np.array(ld); print('     leading term h0 prod L_j / P: %.3f [%.3f, %.3f]' % (np.median(a), a.min(), a.max()))
d = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dev_out.json')))
dp = np.array([x['pred']['diag'] / x['pred']['full'] for x in d['triple']])
print('triples: AFE-diagonal / recipe over [100,1000]: median %.4f [%.4f, %.4f]  => off-diagonal (mixed swaps) share of P: %.3f [%.3f, %.3f]' % (np.median(dp), dp.min(), dp.max(), np.median(1 - dp), (1 - dp).min(), (1 - dp).max()))
pp = np.array([x['diag'] / x['pred'] for x in d['pair']])
print('pairs (dev, |D| 1000-1300): diag/recipe median %.4f [%.4f, %.4f]' % (np.median(pp), pp.min(), pp.max()))
zp = np.array([x['diag'] / x['pred'] for x in d['zetaL']])
print('zeta*L: diag/recipe median %.4f [%.4f, %.4f]' % (np.median(zp), zp.min(), zp.max()))
a3 = np.array([x['A3_over_Aij'] for x in d['triple']])
print('three-body factor h0/h0_no3 = A3(0)/prod A_ij(0): median %.3f [%.3f, %.3f]' % (np.median(a3), a3.min(), a3.max()))
h0 = np.array([x['h0'] for x in d['triple']]); print('h0 over dev triples: %.4g .. %.4g' % (h0.min(), h0.max()))
# closed form for p = 2 aligned: sum C(n+2,2)^2 y^n = (1+4y+y^2)/(1-y)^5 at y = 1/2
y = 0.5; R = (1 + 4 * y + y * y) / (1 - y) ** 5; A3 = R * (1 - 0.5) ** 9; Aij = (1 - 0.25) ** 3
print('p = 2, chi_1(2) = chi_2(2) = chi_3(2): A_3,2 = %.6f, prod A_ij,2 = %.6f, ratio %.6f' % (A3, Aij, A3 / Aij))
