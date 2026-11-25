"""(B, retrospective) regress the measured pair covariances on swapped-term predictor variants (registered design,
SWAP column replaced), jackknife over L-functions. Variants from s4_asym.npz; frozen = swap_pred_NN.npz."""
import sys, math
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__file__))
from mylib import *
from myreg import fit_all, jackknife

S2 = np.load(OUT + '/s2_pairs.npz'); S1 = np.load(OUT + '/s1.npz'); A4 = np.load(OUT + '/s4_asym.npz')
Cp, PSI, REST, FULL, SW, I, J = S2['Cp'], S2['PSI'], S2['REST'], S2['FULL'], S2['SW'], S2['I'], S2['J']
AT, BT = S1['ATAB'], S1['BTAB']; near = A4['near']


def expand(v):
    out = np.zeros_like(SW); out[:, near] = v.T; return out


V = {'frozen (Euler 5e4, full integral)': SW}
for k in ('a2_200k', 'a3_200k', 'a3_50k', 'h1', 'hconst'):
    V[k] = expand(A4[k])
V['resonant part = frozen - a3_50k'] = SW - V['a3_50k']
for name, X in V.items():
    full, se, coefs, cse, cr = jackknife(Cp, PSI, REST, FULL, X, AT, BT, I, J)
    num = (coefs[:, 10] * (X ** 2).sum(1)).sum() / (X ** 2).sum()
    print('%-36s gamma_sw %.3f +- %.3f | sum-of-squares weighted %.3f | rms per w %.4f..%.4f' % (name, full['gamma_sw'], se['gamma_sw'], num, np.sqrt((X[:, near] ** 2).mean(1)).min(), np.sqrt((X[:, near] ** 2).mean(1)).max()))
    print('     per w: ' + ' '.join('%.2f' % c for c in coefs[:, 10]))
    print('     se   : ' + ' '.join('%.2f' % c for c in cse[:, 10]))
    print('     gamma_2 %.4f gamma_full %.4f small_w %.4f icpt mean %.2e' % (full['gamma_2'], full['gamma_full'], full['small_w'], full['icpt_mean']))
# joint fits
print('\njoint fits (both columns in the per-width OLS; mean over widths of each coefficient, jackknife se):')
for a, b in (('a3_200k', 'resonant part = frozen - a3_50k'), ('h1', 'hconst'), ('h1', 'a3_200k'), ('frozen (Euler 5e4, full integral)', 'a3_200k')):
    Xa = V[a]; Xb = V[b] - (V['h1'] if (a == 'h1' and b == 'hconst') else 0)
    full, se, coefs, cse, cr = jackknife(Cp, PSI, REST, FULL, Xa, AT, BT, I, J, extra=[Xb])
    cb = cr[:, :, 11].mean(1); nj = len(cr); seb = math.sqrt((nj - 1) / nj * ((cb - cb.mean()) ** 2).sum())
    print('  %-36s %.3f +- %.3f   %-36s %.3f +- %.3f' % (a, full['gamma_sw'], se['gamma_sw'], b + (' - h1' if (a == 'h1' and b == 'hconst') else ''), coefs[:, 11].mean(), seb))
    print('     per w (first): ' + ' '.join('%.2f' % c for c in coefs[:, 10]))
    print('     per w (second): ' + ' '.join('%.2f' % c for c in coefs[:, 11]))
