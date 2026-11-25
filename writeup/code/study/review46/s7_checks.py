"""(A2, A3) checks: intercept by parity class; contribution of the non-zero means of dn; the expected (finite-range)
intercepts; split-half check of the jackknife se; correlation of readings across widths."""
import sys, math, json
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__file__))
from mylib import *
from myreg import fit_all, jackknife

S1 = np.load(OUT + '/s1.npz'); S2 = np.load(OUT + '/s2_pairs.npz'); Ds = np.array(S1['Ds'])
Cp, PSI, REST, FULL, SW, I, J = S2['Cp'], S2['PSI'], S2['REST'], S2['FULL'], S2['SW'], S2['I'], S2['J']
AT, BT = S1['ATAB'], S1['BTAB']; MU = S1['MU']
JS = json.load(open(REPO + '/writeup/code/study/lowzero/cross/predictions_cross.json'))
E = np.array(JS['expected_readings']['NN']['H_law']['per_w_icpt'])
# 1. product of means of dn
pm = np.array([(MU[k][I] * MU[k][J]).mean() for k in range(16)])
print('mean over pairs of mu_i mu_j per w:', ' '.join('%+.1e' % v for v in pm))
print('rms mu_i per w:', ' '.join('%.4f' % v for v in np.sqrt((MU ** 2).mean(1))), ' mean mu:', ' '.join('%+.4f' % v for v in MU.mean(1)))
# 2. intercept by parity class: replace the constant by three class indicators; subtract expected intercept
s1 = np.sign(Ds[I]); s2 = np.sign(Ds[J]); cls = {'++': (s1 > 0) & (s2 > 0), '+-': s1 != s2, '--': (s1 < 0) & (s2 < 0)}
Ccorr = Cp - E[:, None]          # remove the law's own finite-range intercept
extra = [np.repeat(cls[c][None, :].astype(float), 16, 0) for c in ('+-', '--')]
full, se, coefs, cse, cr = jackknife(Ccorr, PSI, REST, FULL, SW, AT, BT, I, J, extra=extra)
nj = len(cr)
def jse(v):
    return math.sqrt((nj - 1) / nj * ((v - v.mean()) ** 2).sum())
b0 = cr[:, :, 0].mean(1); b1 = cr[:, :, 11].mean(1); b2 = cr[:, :, 12].mean(1)
print('intercept minus law, by parity (mean over widths, jackknife se): ++ %+.2e +- %.1e, +- %+.2e +- %.1e, -- %+.2e +- %.1e' % (
    coefs[:, 0].mean(), jse(b0), (coefs[:, 0] + coefs[:, 11]).mean(), jse(b0 + b1), (coefs[:, 0] + coefs[:, 12]).mean(), jse(b0 + b2)))
print('   counts ++ %d, +- %d, -- %d' % tuple(cls[c].sum() for c in ('++', '+-', '--')))
# 3. intercept with the b_p coefficients fixed at the law (subtract B_T b_p and A_T psi_p and REST and FULL-INF structure):
res = Cp - FULL
print('mean over pairs of C - FULL (law, finite range) per w:', ' '.join('%+.1e' % v for v in res.mean(1)))
full0, se0, c0, cse0, cr0 = jackknife(res, PSI, REST, FULL, SW, AT, BT, I, J)
print('   mean over w of that mean: %+.2e +- %.1e (jackknife)' % (res.mean(), jse(np.array([ (res[:, ~((I == f) | (J == f))]).mean() for f in range(len(Ds))]))))
# 4. split-half check of the jackknife se
full, se, coefs, cse, cr = jackknife(Cp, PSI, REST, FULL, SW, AT, BT, I, J)
rng = np.random.default_rng(3); keys = ('gamma_2', 'gamma_3', 'gamma_5', 'gamma_7', 'gamma_rest', 'gamma_full', 'small_w', 'gamma_sw', 'icpt_mean')
D = {k: [] for k in keys}
for r in range(40):
    perm = rng.permutation(len(Ds)); A_ = set(perm[:92].tolist())
    inA = np.array([(i in A_) and (j in A_) for i, j in zip(I, J)]); inB = np.array([(i not in A_) and (j not in A_) for i, j in zip(I, J)])
    rA, _ = fit_all(Cp[:, inA], PSI[inA], REST[:, inA], FULL[:, inA], SW[:, inA], AT, BT)
    rB, _ = fit_all(Cp[:, inB], PSI[inB], REST[:, inB], FULL[:, inB], SW[:, inB], AT, BT)
    for k in keys: D[k].append(rA[k] - rB[k])
print('split-half: sd(rA - rB)/(2 se_jack) [~1 if the jackknife se is right; var ~ 1/n_functions]:')
print('   ' + ', '.join('%s %.2f' % (k, np.std(D[k]) / (2 * se[k])) for k in keys))
# 5. intercepts of the explicit-formula synthetic (n <= 1e4, no zeros) against the data, per width
Sy = np.load(OUT + '/s6_synth_10000.npz')
rs, cs = fit_all(Sy['Csyn'][:, I, J], PSI, REST, FULL, SW, AT, BT)
rr, cc = fit_all(Cp, PSI, REST, FULL, SW, AT, BT)
print('intercept per w: data  ' + ' '.join('%+.1e' % v for v in cc[:, 0]))
print('                 synth ' + ' '.join('%+.1e' % v for v in cs[:, 0]))
print('                 law (frozen, n<=300 finite range) ' + ' '.join('%+.1e' % v for v in E))
rm, cm = fit_all((S1['C'] - Sy['Csyn'])[:, I, J], PSI, REST, FULL, SW, AT, BT)
print('                 data - synth ' + ' '.join('%+.1e' % v for v in cm[:, 0]) + '   mean %+.1e' % cm[:, 0].mean())
