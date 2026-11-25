"""Pilot on NON-reserved complex zeros (lzeros/q7_k1.txt = Conrey chi_7(3,.), order 6, odd, |gamma| <= 1100;
lzeros/q13_k1.txt = chi_13(2,.), order 12, odd, |gamma| <= 800; identified by |L(1/2 + i gamma)| < 1e-12):
conventions of the one-point law (eps, both half-lines), the resonance predictor, and the phase law of X12 on the
6 pairs among chi7, conj chi7, chi13, conj chi13 (conj chi = the mirrored negative zeros of chi)."""
import sys, math, os
import numpy as np
from scipy.special import loggamma
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hlib import *

ZD = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'data', 'riemann', 'lzeros')) + '/'
pari = pari_instance()


def theta(t, q, a):
    t = np.asarray(t, float)
    return t / 2 * math.log(q / math.pi) + np.imag(loggamma(0.25 + a / 2 + 0.5j * t))


def resonance(chi_p, pr, q, a, w, nmax):
    """R_w = -(i^a/sqrt q) sum_{n = p^k <= nmax} chi(p)^k ln p e(-n/q) w(2 pi n / q)"""
    tot = 0j
    for k in range(1, 40):
        m = pr.astype(float) ** k <= nmax
        if not m.any(): break
        n = pr[m].astype(np.int64) ** k
        tot += np.sum(chi_p[m] ** k * np.log(pr[m]) * np.exp(-2j * math.pi * (n % q) / q) * w(2 * math.pi * n / q))
    return -(1j ** a / math.sqrt(q)) * tot


TT1, TT2 = 100.0, 790.0
def hann(t):
    x = (np.asarray(t, float) - TT1) / (TT2 - TT1); return np.where((x > 0) & (x <= 1), np.sin(math.pi * x) ** 2, 0.0)
def sharp(t):
    t = np.asarray(t, float); return ((t > TT1) & (t <= TT2)).astype(float)


funcs = []      # (label, q, a, eps, chi table at primes, zeros > 0)
P = 30000
for N, m, f in ((7, 3, 'q7_k1.txt'), (13, 2, 'q13_k1.txt')):
    z = np.loadtxt(ZD + f)
    q, a, eps, eps_tau = conrey_info(pari, N, m)
    pr, c = conrey_table(pari, N, m, int(q * TT2 / (2 * math.pi)) + 10)
    print('chi_%d(%d,.): q %d a %d eps %s (tau: %s) kappa %.4f; zeros + %d, - %d' % (N, m, q, a, np.round(eps, 6), np.round(eps_tau, 6), phi_ratio(q), (z > 0).sum(), (z < 0).sum()))
    funcs.append(('chi%d' % N, q, a, eps, pr, c, np.sort(z[z > 0])))
    funcs.append(('conj chi%d' % N, q, a, np.conj(eps), pr, np.conj(c), np.sort(-z[z < 0])))
for lab, q, a, eps, pr, c, g in funcs:
    t = np.linspace(1, g.max() - 1, 200); s = np.searchsorted(g, t) - theta(t, q, a) / math.pi
    u = g / (2 * math.pi) * np.log(q * g / (2 * math.pi * math.e))
    kap = phi_ratio(q); ph = np.exp(2j * theta(g, q, a))
    out = '%-10s N - theta/pi: mean %+.2f, range [%+.2f, %+.2f], largest unfolded gap %.2f |' % (lab, s.mean(), s.min(), s.max(), np.diff(u[g > 20]).max())
    for nm, w, wint in (('sharp', sharp, TT2 - TT1), ('hann', hann, (TT2 - TT1) / 2)):
        G = np.sum(ph * w(g)); M = -eps * kap * wint / (2 * math.pi); R = resonance(c, pr, q, a, w, q * TT2 / (2 * math.pi) + 1)
        out += ' %s: G %.2f%+.2fi, law -eps k T/2pi %.2f%+.2fi, G/law %.3f%+.3fi, R %.2f%+.2fi, |G-R| %.2f |' % (nm, G.real, G.imag, M.real, M.imag, (G / M).real, (G / M).imag, R.real, R.imag, abs(G - R))
    print(out)

# X12 for the 6 pairs, window [TT1 + 6, TT2 - 6] (all partners present)
W1, W2 = TT1 + 6, TT2 - 6
for sig in (0.5, 1.0):
    s = sig / math.sqrt(2); dt = 0.02; tg = np.arange(W1, W2 + dt / 2, dt); T = W2 - W1
    Gm = []
    for lab, q, a, eps, pr, c, g in funcs:
        b = g[(g > W1 - 6 * s) & (g < W2 + 6 * s)]; ph = np.exp(2j * theta(b, q, a))
        row = np.zeros(len(tg), complex); cidx = np.round((b - W1) / dt).astype(int); half = int(6 * s / dt) + 1
        for off in range(-half, half + 1):
            j = cidx + off; ok = (j >= 0) & (j < len(tg))
            np.add.at(row, j[ok], ph[ok] * np.exp(-(tg[j[ok]] - b[ok]) ** 2 / (2 * s * s)))
        Gm.append(row)
    Gm = np.array(Gm); Xm = (Gm @ Gm.conj().T) * dt / (math.sqrt(math.pi) * s * T)
    C = np.array([f[5][:len(primes_upto(P))] if len(f[5]) >= len(primes_upto(P)) else None for f in funcs], dtype=object)
    pr30 = primes_upto(P)
    tabs = []
    for lab, q, a, eps, pr, c, g in funcs:
        _, cc = conrey_table(pari, 7 if '7' in lab else 13, 3 if '7' in lab else 2, P)
        tabs.append(np.conj(cc) if lab.startswith('conj') else cc)
    tabs = np.array(tabs)
    Jm = J_general_all(tabs, pr30, (sig,))[sig]
    print('sigma %.1f, window [%.0f, %.0f]:' % (sig, W1, W2))
    for i in range(4):
        for j in range(i + 1, 4):
            fi, fj = funcs[i], funcs[j]
            pred = predict(phi_ratio(fi[1]), phi_ratio(fj[1]), fi[3], fj[3], Jm[i, j])
            pnr = phi_ratio(fi[1]) * phi_ratio(fj[1]) * Jm[i, j] / (4 * math.pi ** 2)
            pcj = predict(phi_ratio(fi[1]), phi_ratio(fj[1]), fi[3], fj[3], np.conj(Jm[i, j]))
            print('  X(%s, %s) = %.4f%+.4fi   law %.4f%+.4fi   [no root numbers %.4f%+.4fi; conj J %.4f%+.4fi]  J = %.3f%+.3fi' % (
                fi[0], fj[0], Xm[i, j].real, Xm[i, j].imag, pred.real, pred.imag, pnr.real, pnr.imag, pcj.real, pcj.imag, Jm[i, j].real, Jm[i, j].imag))
