# finite-X model (a): Bhargava local factor at l <= 7 replaced by the observed type frequencies of the same fields (independence over l)
import sys, os, time, json; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pred4 as P
from collections import Counter
H = os.path.dirname(os.path.abspath(__file__))
ue = np.array([0.0] + [(2 * j + 1) / 24 for j in range(24)]); ye = ue**3
L1 = json.load(open(os.path.join(H, 'look1.json')))
base = P.quartic_model(T=1000.0)
for r2 in (0, 1, 2):
    z = np.load(os.path.join(H, 'data', 'disc_r%d.npz' % r2)); ty = z['types']; n = len(ty)
    ov = {}
    for k, l in enumerate((2, 3, 5, 7)):
        fr = {nm: c / n for nm, c in Counter(ty[:, k]).items()}
        ov[l] = (lambda X, fr=fr, l=l: P.E_local(l, X, fr))
    m_emp = P.quartic_model(T=1000.0, local_override=ov)
    mus = P.SIG_MU[r2]; pe = P.bin_means(m_emp, mus, ye); pb = P.bin_means(base, mus, ye)
    m = np.array(L1[str(r2)]['m']); s = np.array(L1[str(r2)]['s'])
    print('r2=%d: A(1) limit %.4f, empirical-l<=7 %.4f; y^-1/3 coeff limit %.3f emp %.3f' % (r2, P.A1_exact(), P.A1_exact() * np.prod([ov[l](np.array([1.0]))[0] / P.E_local(l, np.array([1.0]))[0] for l in ov]).real if False else float('nan'),
          P.tail_coeffs(base, mus)['0.3333'], P.tail_coeffs(m_emp, mus)['0.3333']))
    for name, pr in (('limit', pb), ('emp l<=7', pe)):
        for lo, hi in ((1, 4), (4, 13), (13, 24), (4, 24)):
            sl = slice(lo, hi); w = 1 / s[sl]**2; a = np.sum(pr[sl] * m[sl] * w) / np.sum(pr[sl]**2 * w); ae = 1 / np.sqrt(np.sum(pr[sl]**2 * w))
            print('   %-9s bins %2d-%2d (u %.3f-%.3f): amplitude %.3f +- %.3f, chi2 %.1f (shape %.1f) / %d' % (name, lo, hi - 1, ue[lo], ue[hi], a, ae, np.sum((m[sl] - pr[sl])**2 * w), np.sum((m[sl] - a * pr[sl])**2 * w), hi - lo))
    print('   pred emp:', np.round(pe[1:], 2).tolist())
