"""Independent recomputation (no import of elib or ecross_test):
 (a) c_E = lim F(1+x)/x by the direct partial product prod_{good p<=P}(1 - lambda(p^2)/p) prod_{p|N}(1 + 1/p), P = 1e6
     (a_p from PARI ellap), against the frozen c (B/L(1,Sym^2 E));
 (b) J for 6 development and 6 fresh pairs: local factors h_p = Y_p/(F_1 F_2) from the Satake parameters by the
     constant term on the unit circle (FFT, 512 points), p <= 2000; first-order lines 2000 < p <= 2e5 (a_p to 2e5);
     trapezoid in u with step 0.02 to 8 sigma; against the frozen J;
 (c) X_ij for 8 development pairs by the direct double sum over zero pairs with the exact window weight
     w(m) = [erf((W2 - m)/s) - erf((W1 - m)/s)]/2 and theta from mpmath loggamma, against run_dev_data;
 (d) the one-point sums of 5 development curves with mpmath theta, against run_dev_data."""
import os, math, json, sys
import numpy as np, cypari2, mpmath as mp
from scipy.special import erf
HERE = os.path.dirname(os.path.abspath(__file__)) + '/'
pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'deg2') + '/'
sets = {s: json.load(open(D + '%s_set.json' % s))['curves'] for s in ('dev', 'fresh')}
out = open(HERE + 'validate_out.txt', 'w')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
def curve(c, M):
    pari('E = ellinit(%s); PP = primes([2, %d]);' % (str(c['ainv']), M))
    ap = np.array([int(x) for x in pari('Vecsmall(apply(p -> ellap(E, p), PP))')], float)
    pr = np.array([int(x) for x in pari('Vecsmall(PP)')], np.int64)
    return pr, ap
# (a)
for sname in ('dev', 'fresh'):
    F = np.load(HERE + 'pred_%s.npz' % sname); ids = [str(s) for s in F['ids']]
    rel = []
    for c in sets[sname][:6]:
        pr, ap = curve(c, 10 ** 6); good = ~np.isin(pr, c['bad'])
        lam2 = ap[good] ** 2 / pr[good]
        direct = math.exp(np.sum(np.log1p(-(lam2 - 1) / pr[good])) + sum(math.log1p(1 / q) for q in c['bad']))
        rel.append(direct / F['c'][ids.index(c['id'])])
    P('(a) %s: direct partial product to 1e6 / frozen c_E for 6 curves: %s' % (sname, ' '.join('%.4f' % r for r in rel)))
# (b)
NFFT = 512; phi = np.exp(2j * math.pi * np.arange(NFFT) / NFFT)
def local(p, a, bad):
    if bad:
        al = np.array([a / math.sqrt(p), 0.0])
    else:
        th = math.acos(max(-1, min(1, a / (2 * math.sqrt(p))))); al = np.array([np.exp(1j * th), np.exp(-1j * th)])
    return al
def Lp(al, X):      # 1/((1 - al0 X)(1 - al1 X))
    return 1 / ((1 - al[0] * X) * (1 - al[1] * X))
def hp(p, a1, b1, a2, b2, u):
    al1, al2 = local(p, a1, b1), local(p, a2, b2); out = np.empty(len(u), complex)
    q = p ** -0.5
    F1 = np.mean(Lp(al1, q * phi) / Lp(al1, q / phi))      # F_p(1/p) = sum lambda(p^j) mu(p^j) p^-j = CT
    F2 = np.mean(Lp(al2, q * phi) / Lp(al2, q / phi))
    for k, uu in enumerate(u):
        y = p ** (1j * uu)
        A = Lp(al1, q * y * phi); B = 1 / Lp(al1, q / y / phi); C = Lp(al2, q / phi); Dd = 1 / Lp(al2, q * phi)
        out[k] = np.mean(A * B * C * Dd)
    return (out / (F1 * F2)).real if False else out / (F1 * F2)
PMAX = 2 * 10 ** 5
for sname, pairs in (('dev', [(0, 1), (2, 5), (3, 17), (6, 7), (10, 19), (12, 13)]), ('fresh', [(0, 1), (5, 77), (20, 30), (44, 100), (60, 137), (90, 91)])):
    F = np.load(HERE + 'pred_%s.npz' % sname); ids = [str(s) for s in F['ids']]
    byid = {c['id']: c for c in sets[sname]}
    I, J = F['I'], F['J']
    for (i, j) in pairs:
        ci, cj = byid[ids[i]], byid[ids[j]]
        pr, a1 = curve(ci, PMAX); _, a2 = curve(cj, PMAX)
        u = np.arange(0, 8.0001, 0.02); lnH = np.zeros(len(u), complex)
        for k, p in enumerate(pr):
            p = int(p); b1 = p in ci['bad']; b2 = p in cj['bad']
            if p <= 2000 or b1 or b2:
                lnH += np.log(hp(p, a1[k], b1, a2[k], b2, u))
            else:
                lnH += 2 * (a1[k] / math.sqrt(p)) * (a2[k] / math.sqrt(p)) * np.cos(u * math.log(p)) / p
        H = np.exp(lnH); m = int(np.nonzero((I == i) & (J == j))[0][0])
        for s in (0.5, 1.0):
            K = np.exp(-u ** 2 / (2 * s * s)); wts = np.full(len(u), 0.02); wts[0] = 0.01
            Jv = 2 * np.sum(H.real * K * wts)
            P('(b) %s pair (%s, %s) sigma %.1f: J independent %.5f, frozen %.5f, diff %.1e; pred diff %.1e' % (
                sname, ids[i], ids[j], s, Jv, F['J_%.1f' % s][m], Jv - F['J_%.1f' % s][m],
                (Jv - F['J_%.1f' % s][m]) * F['c'][i] * F['c'][j] / (4 * math.pi ** 2)))
# (c), (d)
F = np.load(HERE + 'pred_dev.npz'); ids = [str(s) for s in F['ids']]; I, J = F['I'], F['J']
Z = {k: np.sort(np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'r148', 'dev_zeros.npz'))[k]) for k in ids}
def theta(N, t):
    return np.array([float(tt / 2 * mp.log(N) - tt * mp.log(2 * mp.pi) + mp.im(mp.loggamma(1 + 1j * tt))) for tt in t])
W1, W2 = 20.0, 396.0; T = W2 - W1
for s_ in (0.5, 1.0):
    X = np.load(HERE + 'run_dev_data/X_%.1f.npy' % s_); s = s_ / math.sqrt(2); md = 0
    for (i, j) in [(0, 1), (2, 5), (3, 17), (6, 7), (10, 19), (12, 13), (4, 9), (8, 15)]:
        gi = Z[ids[i]]; gj = Z[ids[j]]; gi = gi[(gi > W1 - 7 * s) & (gi < W2 + 7 * s)]; gj = gj[(gj > W1 - 7 * s) & (gj < W2 + 7 * s)]
        pi_ = np.exp(2j * theta(int(F['N'][i]), gi)); pj = np.exp(2j * theta(int(F['N'][j]), gj))
        U = gj[None, :] - gi[:, None]; Mm = (gj[None, :] + gi[:, None]) / 2
        wgt = 0.5 * (erf((W2 - Mm) / s) - erf((W1 - Mm) / s))
        x = np.sum(pi_[:, None] * np.conj(pj)[None, :] * np.exp(-U ** 2 / (2 * s_ * s_)) * wgt) / T
        md = max(md, abs(x - X[i, j]))
    P('(c) sigma %.1f: X by direct double sum (exact window, mpmath theta) vs ecross_test grid: max |diff| over 8 dev pairs %.1e' % (s_, md))
O = np.load(HERE + 'run_dev_data/onepoint.npz'); md = 0
for k in range(5):
    g = Z[ids[k]]; b = g[(g > 20) & (g <= 400)]; G0 = np.sum(np.exp(2j * theta(int(F['N'][k]), b)))
    md = max(md, abs(G0 - O['Gl'][k, 0, 0]))
P('(d) one-point sums (r = 1, sharp) with mpmath theta vs ecross_test: max |diff| over 5 dev curves %.1e' % md)
out.close()
