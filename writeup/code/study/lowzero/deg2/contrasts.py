"""Per-prime contrasts of Sigma2(l) at fixed L (single prime changed, everything else held): the diagonal term split by
traversal (k = 1, k = 2, k >= 3; closed form (2/pi^2) dW_p(k) sin^2(pi k l ln p/L)/(k^2 p^k)) and the off-diagonal change
2 int_0^l (l - r) (2/L^2) Re[e^{-2 pi i r} K(i eps) (T_p'(i eps)/T_p(i eps) - 1)] dr, with K of a host curve.
usage: contrasts.py"""
import sys, os, math, pickle
import numpy as np, cypari2
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deg2lib as D
pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)
H = os.path.dirname(os.path.abspath(__file__)); Z = pickle.load(open(os.path.join(H, 'zgrid.pkl'), 'rb'))
LSHOW = (0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0); IX = [int(round(l / D.W)) - 1 for l in LSHOW]
out = []
def rep(s): print(s, flush=True); out.append(s)


def W(p, ap, k):
    if ap is None: return p ** (-k)                       # multiplicative
    th = math.acos(ap / math.sqrt(p) / 2); return (2 * math.cos(k * th)) ** 2


def lnT(p, ap, x):
    """ln T_p(x) = ln O_p(x) - 2 ln(1 - 1/p)"""
    if ap is None:
        y = np.exp(x * math.log(p)); O = 1 + (1 - 1 / y) ** 2 * y ** 2 * np.exp(-(1 + x) * math.log(p)) / p / (1 - y ** 2 * np.exp(-(1 + x) * math.log(p)) / p)
        return np.log(O) - 2 * math.log1p(-1 / p)
    d, la, _, _, _ = D.local_good(p, ap / math.sqrt(p), x)
    lam = ap / math.sqrt(p); th = math.acos(lam / 2); a2 = np.exp(2j * th); b2 = np.conj(a2); lp = math.log(p)
    zp = np.exp(-(1 + x) * lp); zm = np.exp(-(1 - x) * lp); z0 = 1 / p
    lnden = -(np.log(1 - zp) + np.log(1 - zm)) - (np.log(1 - a2 * zp) + np.log(1 - zp) + np.log(1 - b2 * zp)) \
            - (np.log(1 - a2 * zm) + np.log(1 - zm) + np.log(1 - b2 * zm)) + 2 * (np.log(1 - a2 * z0) + np.log(1 - z0) + np.log(1 - b2 * z0)).real
    return la + lnden


def diag(p, apA, apB, L, ks):
    l = D.LGRID
    return sum((2 / math.pi ** 2) * (W(p, apA, k) - W(p, apB, k)) * np.sin(math.pi * k * l * math.log(p) / L) ** 2 / (k * k * p ** k)
               for k in ks if p ** (-k) > 1e-18)


def off(c, p, apA, apB, L):
    r = D.R_NODES; e = 2 * math.pi * r / L; x = 1j * e
    Kfull = c('KSA', e) + Z('Z2m1', e) * np.exp(c('lnSA', e)) + 1 / np.where(e > 0, e, 1e-300) ** 2
    hostap = None if p in c.bad else c.ap[p]
    rA = lnT(p, apA, x) - lnT(p, hostap, x); rB = lnT(p, apB, x) - lnT(p, hostap, x)
    dK = Kfull * (np.exp(rA) - np.exp(rB))
    f = (2 / L ** 2) * (np.exp(-2j * math.pi * r) * dK).real
    return D.cum_sigma(f)


for host in ([1, -1, 1, -10, 14], [1, 1, 1, -15, 16]):
    c = D.Curve(pari, host)
    rep('host curve %s, N = %d (its own K(i eps) carries the other primes)' % (host, c.N))
    for L in (12.0, 14.0):
        rep(' L = %.1f; l = ' % L + ' '.join('%7.1f' % l for l in LSHOW))
        for p, pairs in ((2, [(0, 2), (0, 1), (1, 2), (None, 1)]), (3, [(0, 3), (0, 1), (1, 2), (2, 3), (None, 1)])):
            for A_, B_ in pairs:
                nm = 'p=%d |a|=%s vs %s' % (p, 'bad' if A_ is None else A_, B_)
                d1 = diag(p, A_, B_, L, [1]); d2 = diag(p, A_, B_, L, [2]); d3 = diag(p, A_, B_, L, range(3, 61))
                o = off(c, p, A_, B_, L); tot = d1 + d2 + d3 + o
                rep('  %-18s k=1  ' % nm + ' '.join('%+7.4f' % d1[i] for i in IX))
                rep('  %-18s k=2  ' % '' + ' '.join('%+7.4f' % d2[i] for i in IX))
                rep('  %-18s k>=3 ' % '' + ' '.join('%+7.4f' % d3[i] for i in IX))
                rep('  %-18s off  ' % '' + ' '.join('%+7.4f' % o[i] for i in IX))
                rep('  %-18s total' % '' + ' '.join('%+7.4f' % tot[i] for i in IX))
open(os.path.join(H, 'contrasts_out.txt'), 'w').write('\n'.join(out) + '\n')
