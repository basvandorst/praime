"""Power of the registered kappa3 readings, from the NOISE of the used zero sets (residual covariance after the indicator
regression, noise_clean_resid.npz; no measured mean enters) and the frozen predictions (pred3_*.json).
Noise model for an L-function: the residual covariance across l of a reference sample, scaled by (reference mean zeros /
own zeros). Expected values = the predictions. For every reading: its expected standard error and the expected z
against each alternative (GUE, Q0, two-point only).
usage: power3.py -> power3_out.txt, power3.json"""
import os, json, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(HERE, 'noise_clean_resid.npz'))
LS = [0.5 * i for i in range(1, 17)]; PR = (2, 3, 5, 7, 11, 13); P3 = (2, 3, 5)
J = json.load(open(os.path.join(HERE, 'kappa3_tables.json'))); IXT = [int(round(l / 0.25)) - 1 for l in LS]
GUE = np.array(J['gue_k3'])[IXT]; GUEF = np.array(J['gue_kf3'])[IXT]


def zeros_band(D, lo, hi):
    q = abs(D); f = lambda t: t / (2 * math.pi) * math.log(q * t / (2 * math.pi * math.e)); return f(hi) - f(lo)


def cov(ref, key):
    R = Z['%s_%s' % (ref, key)]; return R.T @ R / float(Z[ref + '_dof']), float(Z[ref + '_NZ'].mean())


def analyse(names, refs, label):
    rows = {k: [] for k in ('D', 'full', 'Q0', 'twopt', 'kf3', 'diagonly', 'Wk', 'Wf')}
    for nm, rf in zip(names, refs):
        d = json.load(open(os.path.join(HERE, 'pred3_%s.json' % nm))); Ck, nz = cov(rf, 'K3'); Cf, _ = cov(rf, 'KF'); lo, hi = d['height']
        for i, D in enumerate(d['D']):
            s = nz / zeros_band(D, lo, hi)
            rows['D'].append(D); rows['full'].append(d['full'][i]); rows['Q0'].append(d['Q0'][i]); rows['twopt'].append(d['twopt'][i])
            rows['kf3'].append(d['kf3'][i]); rows['diagonly'].append(d['diagonly'][i]); rows['Wk'].append(Ck * s); rows['Wf'].append(Cf * s)
    F, Q0, TW, KF, DG = (np.array(rows[k]) for k in ('full', 'Q0', 'twopt', 'kf3', 'diagonly')); n = len(F)
    Ck = np.mean(rows['Wk'], axis=0); Cf = np.mean(rows['Wf'], axis=0); sd = np.sqrt(np.diag(Ck)); sdf = np.sqrt(np.diag(Cf))
    X = np.array([[1] + [int(abs(D) % p == 0) for p in PR] for D in rows['D']], float); A = np.linalg.pinv(X.T @ X); H = A @ X.T
    groups = [int(X[:, j + 1].sum()) for j in range(len(PR))]
    out = {'label': label, 'n': n, 'groups': groups}
    print('\n==== %s: n = %d, groups (2,3,5,7,11,13) %s' % (label, n, groups))
    print('  model noise sd of kappa3 per L-function, l = 0.5..8: ' + ' '.join('%.4f' % x for x in sd))
    # T1 level
    Cl = Ck / n; Cli = np.linalg.inv(Cl); se1 = np.sqrt(np.diag(Cl)); m = F.mean(0)
    print('T1 level: mean predicted kappa3 (full)      ' + ' '.join('%+.4f' % x for x in m))
    print('          se of the mean                    ' + ' '.join('%.5f' % x for x in se1))
    out['T1_se'] = se1.tolist(); out['T1_full'] = m.tolist()
    for nm, alt in (('GUE', np.tile(GUE, (n, 1))), ('Q0', Q0), ('two-point only', TW)):
        dv = m - alt.mean(0); z = dv / se1; chi = float(dv @ Cli @ dv)
        print('          full - %-15s z per l: %s | GLS chi2 %.0f (z %.0f)' % (nm, ' '.join('%+.0f' % x for x in z), chi, math.sqrt(chi)))
        out['T1_z_' + nm] = z.tolist()
    # T2 arithmetic amplitude of the level
    a = m - GUE; seA = 1 / math.sqrt(a @ Cli @ a); print('T2 gamma_A (level, arithmetic amplitude): se %.4f -> GUE (0) at %.0f se' % (seA, 1 / seA))
    out['T2_se'] = seA
    # T3 per-prime coefficients
    cf, cq, ct = (H @ F)[1:], (H @ Q0)[1:], (H @ TW)[1:]
    se_c = np.sqrt(np.diag(A))[1:, None] * sd[None, :]; out['T3'] = {}
    for j, p in enumerate(PR[:4]):
        Cc = A[j + 1, j + 1] * Ck; Cci = np.linalg.inv(Cc)
        print('T3 p=%-2d coefficient (full)  %s' % (p, ' '.join('%+.4f' % x for x in cf[j])))
        print('        Delta(full - Q0)    %s' % ' '.join('%+.4f' % x for x in cf[j] - cq[j]))
        print('        se                  %s' % ' '.join('%.4f' % x for x in se_c[j]))
        zz = {}
        for nm, alt in (('GUE', 0 * cf[j]), ('Q0', cq[j]), ('two-point only', ct[j])):
            dv = cf[j] - alt; zz[nm] = (float(np.abs(dv / se_c[j]).max()), float(math.sqrt(dv @ Cci @ dv)))
            print('        vs %-15s max |z| per cell %.1f, combined over 16 lengths %.1f' % (nm, zz[nm][0], zz[nm][1]))
        out['T3'][p] = {'coef': cf[j].tolist(), 'se': se_c[j].tolist(), 'z': zz}
    # T4 gamma_Q: (a) from the 48 cells p = 2, 3, 5 (diag weights); (b) pooled per L-function
    dq = (cf - cq)[:3]; w = 1 / se_c[:3] ** 2; seQa = 1 / math.sqrt((dq ** 2 * w).sum())
    x = F - Q0; x = x - x.mean(0); wl = 1 / np.diag(Ck); seQb = 1 / math.sqrt((x ** 2 * wl).sum())
    print('T4 gamma_Q (triple-orbit term removed by p | k): cells p=2,3,5: se %.4f (Q0 at %.0f se); pooled per L-function: se %.4f (%.0f se)' % (seQa, 1 / seQa, seQb, 1 / seQb))
    out['T4_se_cells'] = seQa; out['T4_se_pooled'] = seQb
    dt = (cf - ct)[:3]; se2 = 1 / math.sqrt((dt ** 2 * w).sum()); print('    gamma_2pt from the same cells: se %.5f (two-point-only at %.0f se)' % (se2, 1 / se2))
    cd = (H @ DG)[1:]; do = (cf - cd)[:3]; seO = 1 / math.sqrt((do ** 2 * w).sum()); out['T4_se_off'] = seO
    print('    gamma_off (bad-prime off-diagonal remainder 3dSigma2 + dk[3]_P): cells p=2,3,5: se %.4f (diagonal-only at %.0f se)' % (seO, 1 / seO))
    for j, p in enumerate(P3):
        print('       p=%d remainder coefficient (full - diagonly) %s' % (p, ' '.join('%+.4f' % x for x in do[j])))
    # T5 factorial cumulant
    Clf = Cf / n; mf = KF.mean(0); sef = np.sqrt(np.diag(Clf)); KF1 = KF - (F - TW)   # k[3] of k = 1 at the same L (two-point-only law)
    print('T5 k[3] level: mean predicted (full) minus GUE ' + ' '.join('%+.3f' % x for x in mf - GUEF))
    print('          se ' + ' '.join('%.4f' % x for x in sef))
    print('          full vs GUE: z per l ' + ' '.join('%+.0f' % x for x in (mf - GUEF) / sef))
    cff, cf1 = (H @ KF)[1:], (H @ KF1)[1:]; se_cf = np.sqrt(np.diag(A))[1:, None] * sdf[None, :]
    for j, p in enumerate(PR[:3]):
        print('   per-prime k[3] p=%d: full %s | vs two-point-only max |z| %.0f' % (p, ' '.join('%+.3f' % x for x in cff[j][::3]), np.abs((cff[j] - cf1[j]) / se_cf[j]).max()))
    return out


res = {}
res['546'] = analyse(['R120', 'B', 'Dn'], ['R120', 'B', 'Dn'], 'pooled 546 real L-functions, 100-600 (zeros on disk; kappa3 means never computed)')
res['C'] = analyse(['C'], ['C'], 'C: R120 discriminants at 600-1000 (zeros on disk)')
res['Fn'] = analyse(['Fn'], ['Dn'], 'Fn (fresh): 1000 < |D| <= 1300, 100-600; noise model from Dn')
res['Bp'] = analyse(['Bp'], ['C'], 'Bp (fresh): B discriminants at 600-1000; noise model from C')
json.dump(res, open(os.path.join(HERE, 'power3.json'), 'w'))
