"""Independent recomputation of key numbers of predict.py (no code shared with bklib.py):
- Psi_off in the paper's literal form (1 + (p^{ie/2} - p^{-ie/2})^2/(p - p^{-ie}))^{-1}
- |zeta(1+ie)|^2 and (zeta'/zeta)' by mpmath at every node (the latter by mpmath numerical differentiation of log zeta)
- Phi_off as a direct product over primes < 1e5, d^2 ln Phi_diag by central finite differences
- integrals by scipy.quad on each half-period
Prints the values next to those of bklib (imported only for the comparison column). usage: validate.py"""
import os, sys, math, json
import numpy as np, mpmath as mp
from scipy import integrate
mp.mp.dps = 25
HERE = os.path.dirname(os.path.abspath(__file__))

def primes(n):
    s = [True] * (n + 1); s[0] = s[1] = False
    for i in range(2, int(n ** .5) + 1):
        if s[i]:
            for j in range(i * i, n + 1, i): s[j] = False
    return np.array([i for i in range(n + 1) if s[i]], float)
P = primes(100000); LP = np.log(P)

def zeta_abs2(e): return float(abs(mp.zeta(mp.mpc(1, e))) ** 2)
def phi_off(e): return np.prod(1 - (np.exp(1j * e * LP) - 1) ** 2 / (P - 1) ** 2)
def psi_off_paper(p, e):
    a = p ** (0.5j * e) - p ** (-0.5j * e)
    return 1 / (1 + a * a / (p - p ** (-1j * e)))
def ln_phi_diag(e):
    s = 0.0
    for m in range(1, 60):
        sel = P ** m < 1e17
        if not sel.any(): break
        s += np.sum((1 - m) / (m * m * P[sel] ** m) * np.cos(m * e * LP[sel]))
    return 2 * s
def d2(f, e, h=1e-3): return (f(e + h) - 2 * f(e) + f(e - h)) / h ** 2
def ln_psi(B, e): return -2 * sum(sum(math.cos(m * e * math.log(p)) / (m * m * p ** m) for m in range(1, 64) if p ** m < 1e18) for p in B)

def dR_off(r, L, p, rest=()):
    """change of R(r) - 1 (unfolded) from adding p to the bad set `rest`, off-diagonal part only"""
    e = 2 * math.pi * r / L
    X0 = zeta_abs2(e) * phi_off(e) * np.prod([psi_off_paper(q, e) for q in rest]) if rest else zeta_abs2(e) * phi_off(e)
    dX = X0 * (psi_off_paper(p, e) - 1)
    return (4 * math.pi ** 2 / L ** 2) * (1 / (4 * math.pi ** 2)) * 2 * (np.exp(1j * e * L) * dX).real

def sig_int(f, l):
    """2 int_0^l (l - r) f(r) dr, quad per half-period"""
    edges = list(np.arange(0, l, 0.5)) + [l]
    return sum(integrate.quad(lambda r: 2 * (l - r) * f(r), a, b, limit=100, epsabs=1e-12)[0] for a, b in zip(edges[:-1], edges[1:]) if b > a)

def R_minus_1_zeta(r, L):
    """full R(r) - 1 for zeta (k = 1) from the paper's formulas, direct"""
    e = 2 * math.pi * r / L
    dlz = mp.diff(lambda x: mp.log(mp.zeta(mp.mpc(1, x))), e, 2)          # d^2/de^2 ln zeta(1+ie)
    R2diag = -(1 / (4 * math.pi ** 2)) * (2 * float(mp.re(dlz)) + d2(ln_phi_diag, e))
    R2off = (1 / (4 * math.pi ** 2)) * 2 * (np.exp(1j * e * L) * zeta_abs2(e) * phi_off(e)).real
    return (4 * math.pi ** 2 / L ** 2) * (R2diag + R2off)

sys.path.insert(0, HERE)
import bklib
G = bklib.Grid(eps_max=7.6, h=0.0025, pmax=10 ** 5, verbose=False)
ix = lambda l: int(round(l / bklib.W)) - 1
out = []
print('off-diagonal single prime, L = 9.38: independent | bklib')
for p, l in ((2, 0.5), (2, 1.0), (2, 7.0), (2, 7.5), (3, 4.5), (3, 5.0), (13, 0.5), (13, 3.5)):
    v = sig_int(lambda r: dR_off(r, 9.38, p), l); w = bklib.dsig_off(G, 9.38, p)[ix(l)]
    print('  p=%2d l=%.1f  %+.6f | %+.6f  diff %.1e' % (p, l, v, w, v - w)); out.append(abs(v - w))
print('off-diagonal with other bad primes (rest), L = 9.0: independent | bklib')
for p, rest, l in ((2, (3, 5), 7.0), (2, (3, 5), 7.5), (2, (29,), 7.0), (3, (2, 11), 4.5), (5, (2, 3), 6.0)):
    v = sig_int(lambda r: dR_off(r, 9.0, p, rest), l); w = bklib.dsig_off(G, 9.0, p, rest=rest)[ix(l)]
    print('  p=%2d rest=%s l=%.1f  %+.6f | %+.6f  diff %.1e' % (p, rest, l, v, w, v - w)); out.append(abs(v - w))
print('multi-prime consistency in bklib: Sigma2_full(B) - Sigma2_full(B - {p}) - diag(p) vs dsig_off(p, rest)')
for B, p in (((2, 3, 5), 2), ((2, 3, 5), 3), ((2, 7, 13), 13)):
    rest = tuple(q for q in B if q != p)
    a = bklib.sigma2_full(G, 9.0, B) - bklib.sigma2_full(G, 9.0, rest) - bklib.dsig_diag_grid(p, 9.0); b = bklib.dsig_off(G, 9.0, p, rest=rest)
    print('  B=%s p=%d: max |diff| over l in (0, 8]: %.1e' % (B, p, np.abs(a - b).max())); out.append(np.abs(a - b).max())
print('diagonal per prime: closed form vs numerical integral of dR_diag = -(2/L^2) sum_m (ln p)^2 p^-m cos(2 pi m r ln p/L)')
for p, l in ((2, 7.0), (13, 3.5), (3, 0.5)):
    f = lambda r: -(2 / 9.38 ** 2) * sum(math.log(p) ** 2 * p ** -m * math.cos(2 * math.pi * m * r * math.log(p) / 9.38) for m in range(1, 64) if p ** m < 1e18)
    v = sig_int(f, l); w = bklib.dS_diag(p, 9.38, l)
    print('  p=%2d l=%.1f  %+.6f | %+.6f  diff %.1e' % (p, l, v, w, v - w)); out.append(abs(v - w))
print('zeta, full Sigma2 at fixed L = 8.1728: independent (GUE by quad + correction by quad) | bklib')
L = 8.172796842142144
for l in (0.5, 1.0, 4.5, 8.0):
    gue = l - sig_int(lambda r: (math.sin(math.pi * r) / (math.pi * r)) ** 2 if r > 0 else 1.0, l)
    corr = sig_int(lambda r: R_minus_1_zeta(r, L) + ((math.sin(math.pi * r) / (math.pi * r)) ** 2 if r > 0 else 1.0), l)
    v = gue + corr; w = bklib.sigma2_full(G, L, ())[ix(l)]
    print('  l=%.1f  %.6f | %.6f  diff %.1e' % (l, v, w, v - w)); out.append(abs(v - w))
print('max |independent - bklib| over all checks: %.1e' % max(out))
