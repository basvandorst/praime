"""Numbers for round_draft.txt: the expected value of every reading under every hypothesis (data replaced by that
hypothesis' frozen predictions), with the expected se from the noise model of power3.py. No data enter.
Hypotheses: full, GUE, Q0, twopt (two-point only), diagonly (bad primes act only through the triple-orbit term).
usage: registry.py -> registry_out.txt"""
import os, json, math
import numpy as np
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(HERE, 'noise_clean_resid.npz'))
LS = [0.5 * i for i in range(1, 17)]; PR = (2, 3, 5, 7, 11, 13)
J = json.load(open(os.path.join(HERE, 'kappa3_tables.json'))); IXT = [int(round(l / 0.25)) - 1 for l in LS]
GUE = np.array(J['gue_k3'])[IXT]


def zeros_band(D, lo, hi):
    q = abs(D); f = lambda t: t / (2 * math.pi) * math.log(q * t / (2 * math.pi * math.e)); return f(hi) - f(lo)


for name, ref in (('Fn', 'Dn'), ('Bp', 'C')):
    d = json.load(open(os.path.join(HERE, 'pred3_%s.json' % name))); lo, hi = d['height']; n = len(d['D'])
    R = Z[ref + '_K3']; C0 = R.T @ R / float(Z[ref + '_dof']); nz = float(Z[ref + '_NZ'].mean())
    Ck = np.mean([C0 * nz / zeros_band(D, lo, hi) for D in d['D']], axis=0)
    X = np.array([[1] + [int(abs(D) % p == 0) for p in PR] for D in d['D']], float); A = np.linalg.pinv(X.T @ X); Hm = A @ X.T
    H = {k: np.array(d[k]) for k in ('full', 'Q0', 'twopt', 'diagonly')}; H['GUE'] = np.tile(GUE, (n, 1))
    coef = {k: (Hm @ v)[1:] for k, v in H.items()}
    sd = np.sqrt(np.diag(Ck)); se_c = np.sqrt(np.diag(A))[1:, None] * sd[None, :]; w = 1 / se_c[:3] ** 2
    Cl = Ck / n; Cli = np.linalg.inv(Cl); a = H['full'].mean(0) - GUE
    print('\n=== test %s (n = %d), noise model from %s ===' % (name, n, ref))
    print('T1 level: full mean kappa3   ' + ' '.join('%+.4f' % x for x in H['full'].mean(0)))
    print('          GUE                ' + ' '.join('%+.4f' % x for x in GUE))
    print('          Q0 - full          ' + ' '.join('%+.4f' % x for x in (H['Q0'] - H['full']).mean(0)))
    print('          se                 ' + ' '.join('%.5f' % x for x in np.sqrt(np.diag(Cl))))
    fcrit = (n - 1) * 16 / (n - 16) * stats.f.ppf(0.999, 16, n - 16)
    print('          T^2 0.1%% point %.1f (chi2_16 0.1%% point %.1f)' % (fcrit, stats.chi2.ppf(0.999, 16)))
    row = lambda k, f: ' | '.join('%s %s' % (h, f(h)) for h in ('full', 'GUE', 'Q0', 'twopt', 'diagonly'))
    def gA(h): b = H[h].mean(0) - GUE; return '%.3f' % ((a @ Cli @ b) / (a @ Cli @ a))
    def T1(h): m = (H[h] - H['full']).mean(0); return '%.0f' % (n * 0 + m @ Cli @ m + 16)
    def gam(h, alt):
        dd = (coef['full'] - coef[alt])[:3]; y = (coef[h] - coef[alt])[:3]; return '%.3f' % ((w * y * dd).sum() / (w * dd * dd).sum())
    def se_g(alt):
        dd = (coef['full'] - coef[alt])[:3]; return 1 / math.sqrt((w * dd * dd).sum())
    def chi2p(h, j):
        Cc = A[j + 1, j + 1] * Ck; dv = coef[h][j] - coef['full'][j]; return dv @ np.linalg.inv(Cc) @ dv + 16
    print('T1 expected T^2 (incl. 16 dof):     ' + row('', T1))
    print('T2 gamma_A (se %.4f):              ' % (1 / math.sqrt(a @ Cli @ a)) + row('', gA))
    print('T4 gamma_Q (se %.4f):              ' % se_g('Q0') + row('', lambda h: gam(h, 'Q0')))
    print('T4 gamma_2pt (se %.5f):            ' % se_g('twopt') + row('', lambda h: gam(h, 'twopt')))
    print('T4 gamma_off (se %.4f):            ' % se_g('diagonly') + row('', lambda h: gam(h, 'diagonly')))
    for j, p in enumerate((2, 3, 5)):
        print('T3 chi2_%d against full (expected; 0.1%% point %.1f): ' % (p, stats.chi2.ppf(0.999, 16)) + ' | '.join('%s %.0f' % (h, chi2p(h, j)) for h in ('full', 'GUE', 'Q0', 'twopt', 'diagonly')))
    for j, p in enumerate((2, 3, 5)):
        print('T3 p=%d predicted coefficient (full) ' % p + ' '.join('%+.4f' % x for x in coef['full'][j]))
        print('                   Q0               ' + ' '.join('%+.4f' % x for x in coef['Q0'][j]))
        print('                   diagonly         ' + ' '.join('%+.4f' % x for x in coef['diagonly'][j]))
        print('                   se               ' + ' '.join('%.4f' % x for x in se_c[j]))
    i6 = [LS.index(l) for l in (3, 3.5, 4, 4.5, 5.5)]
    print('T6 p=2 coefficient at l = 3, 3.5, 4, 4.5, 5.5: ' + ' '.join('%+.4f' % coef['full'][0][i] for i in i6) + '  (z ' + ' '.join('%.0f' % (coef['full'][0][i] / se_c[0][i]) for i in i6) + ')')
