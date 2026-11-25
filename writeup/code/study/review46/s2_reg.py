"""Step 2: compare own covariances / predictions with the frozen files; own readings + jackknife; compare with rcross_NN.json."""
import sys, json, math
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__file__))
from mylib import *
from myreg import fit_all, jackknife

CROSS = REPO + '/writeup/code/study/lowzero/cross/'
S = np.load(OUT + '/s1.npz'); Ds = list(S['Ds']); idx = {D: i for i, D in enumerate(Ds)}
P = np.load(CROSS + 'predictions_cross_pairs.npz'); JS = json.load(open(CROSS + 'predictions_cross.json'))
D1 = P['NN_D1']; D2 = P['NN_D2']
ok = np.array([(int(a) in idx) and (int(b) in idx) for a, b in zip(D1, D2)])
D1, D2 = D1[ok], D2[ok]; I = np.array([idx[int(a)] for a in D1]); J = np.array([idx[int(b)] for b in D2])
print('pairs', len(I))
FULLf = P['NN_full'][ok].T.astype(float); RESTf = P['NN_rest'][ok].T.astype(float); PSIf = P['NN_psi'][ok].astype(float)
Cp = S['C'][:, I, J]; FULLm = S['FULL'][:, I, J]; RESTm = S['REST'][:, I, J]
PSIm = (S['X7'][I].astype(float) * S['X7'][J].astype(float))
print('psi identical:', np.array_equal(PSIm, PSIf))
print('FULL own vs frozen: max |diff| %.2e' % np.abs(FULLm - FULLf).max(), ' REST: %.2e' % np.abs(RESTm - RESTf).max())
A_T = np.array(JS['A_T']).T; B_T = np.array(JS['B_T']).T
print('A_T own vs frozen max diff %.2e, B_T %.2e' % (np.abs(S['ATAB'] - A_T).max(), np.abs(S['BTAB'] - B_T).max()))
# swap predictor (frozen)
sw = np.load(CROSS + 'swap_pred_NN.npz'); key = {(int(a), int(b)): c for a, b, c in zip(sw['D1'], sw['D2'], sw['C'])}
SW = np.zeros_like(Cp); nsw = 0
for n_, (a, b) in enumerate(zip(D1, D2)):
    c = key.get((int(a), int(b)), key.get((int(b), int(a))))
    if c is not None: SW[:, n_] = c; nsw += 1
print('pairs with a swap prediction:', nsw, ' sum of squares %.4f' % (SW ** 2).sum())
np.savez(OUT + '/s2_pairs.npz', I=I, J=J, D1=D1, D2=D2, Cp=Cp, FULL=FULLm, REST=RESTm, PSI=PSIm, SW=SW, FULLf=FULLf)

full, se, coefs, cse, cr = jackknife(Cp, PSIm, RESTm, FULLm, SW, S['ATAB'], S['BTAB'], I, J)
R = json.load(open(CROSS + 'rcross_NN.json'))
print('%-11s %10s %8s | %10s %8s' % ('reading', 'own', 'se', 'rcross', 'se'))
for k in ('gamma_2', 'gamma_3', 'gamma_5', 'gamma_7', 'gammaB_2', 'gammaB_3', 'gamma_rest', 'gamma_full', 'small_w', 'gamma_sw'):
    print('%-11s %+10.5f %8.5f | %+10.5f %8.5f' % (k, full[k], se[k], R['readings'][k], R['se'][k]))
print('intercepts own:', ' '.join('%+.6f' % v for v in coefs[:, 0]))
print('intercepts rc :', ' '.join('%+.6f' % v for v in R['readings']['icpt']))
print('per-width swap coefficient:', ' '.join('%.3f' % v for v in coefs[:, 10]), '\n   se:', ' '.join('%.3f' % v for v in cse[:, 10]))
# intercepts against the expected (finite-range law) intercepts
E = np.array(JS['expected_readings']['NN']['H_law']['per_w_icpt'])
d = coefs[:, 0] - E
covI = np.cov(cr[:, :, 0].T, bias=True) * (len(cr) - 1)
print('intercept - expected (law, finite range):', ' '.join('%+.1e' % v for v in d))
print('  z per width:', ' '.join('%+.2f' % v for v in d / cse[:, 0]))
w1 = np.ones(16) / 16; m = d.mean(); sm = math.sqrt(w1 @ covI @ w1)
print('  mean over widths %+.2e +- %.2e (jackknife, widths correlated): z = %.2f' % (m, sm, m / sm))
chi2 = d @ np.linalg.solve(covI, d); print('  Hotelling-type chi2 on 16: %.1f' % chi2)
m0 = coefs[:, 0].mean(); print('  raw mean intercept %+.2e +- %.2e: z = %.2f' % (m0, sm, m0 / sm))
corr = covI / np.sqrt(np.outer(np.diag(covI), np.diag(covI))); print('  mean correlation between widths of the intercept: %.2f' % corr[np.triu_indices(16, 1)].mean())
json.dump({k: float(v) for k, v in full.items()}, open(OUT + '/s2_readings.json', 'w'))
