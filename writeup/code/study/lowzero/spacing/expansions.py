"""Task 1(b): small-distance expansion of the repo's two-point function (bklib, Bogomolny-Keating 2013) for a primitive
character with bad-prime set S, coefficient by coefficient.
  R(r) - 1 = -(sin pi r/pi r)^2 + C(r),  C(r) = (2/L^2)[ g1(eps) + Re(e^{2 pi i r} h(eps)) ],  eps = 2 pi r/L,
  g1 = Re(zeta'/zeta)'(1+i eps) + 1/eps^2 - (1/2)(d2 ln Phi_diag + d2 ln Psi_S),  h = (eps^2|zeta|^2 Phi_off Psi_off,S - 1)/eps^2.
Claim (derived in derivation.md 2.1): g1(0) = -Lambda_S, h(0) = Lambda_S, Im h'(0) = Q_S, so
  C(r) = -(4 Lambda_S/L^2) sin^2(pi r) - (4 pi Q_S/L^3) r sin(2 pi r) + O(L^-4),
  Lambda_S = gamma0^2 + 2 gamma1 + sum_p ln^2 p/(p-1)^2 + sum_{p in S} ln^2 p/(p-1),
  Q_S      = sum_p ln^3 p/(p-1)^2 - sum_{p in S} ln^3 p/(p-1)^2.
"""
import sys, math, json
import numpy as np, mpmath as mp
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/lib')
import loadlibs
bk = loadlibs.load_bklib(); tp = loadlibs.load_tplib(); G = tp.G
TWO_PI = 2 * math.pi

# exact constants
mp.mp.dps = 30
c0 = float(mp.euler ** 2 + 2 * mp.stieltjes(1))
P = bk.primes_upto(2 * 10 ** 7).astype(float); lp = np.log(P)
PM = P[-1]
def tail(a):   # sum_{p > PM} ln^a p/(p-1)^2 ~ int_PM^inf ln^(a-1) x / x^2 dx  (PNT density 1/ln x)
    f = lambda x: mp.log(x) ** (a - 1) / x ** 2
    return float(mp.quad(f, [PM, mp.inf]))
A2 = float(np.sum(lp ** 2 / (P - 1) ** 2)) + tail(2)
A3 = float(np.sum(lp ** 3 / (P - 1) ** 2)) + tail(3)
LAM = c0 + A2; Q = A3
print('gamma0^2 + 2 gamma1 = %.8f' % c0)
print('sum_p ln^2 p/(p-1)^2 = %.6f (tail %.1e)   Lambda_zeta = %.6f   (BBLM: 1.57314)' % (A2, tail(2), LAM))
print('sum_p ln^3 p/(p-1)^2 = Q_zeta = %.6f (tail %.1e)   C = Q/Lambda = %.4f (BBLM: 1.4720)' % (A3, tail(3), Q / LAM))

def lam_S(S): return LAM + sum(math.log(p) ** 2 / (p - 1) for p in S)
def q_S(S): return Q - sum(math.log(p) ** 3 / (p - 1) ** 2 for p in S)

def parts(S, e):
    e = np.asarray(e, float)
    Z2m1 = G('Z2m1', e); Phi = G('Phi', e); Phim1 = G('Phim1', e)
    Psi = np.ones_like(e, dtype=complex); Psim1 = np.zeros_like(e, dtype=complex)
    for p in S:
        fp = bk.f_p(p, e); F = 1 + e ** 2 * fp; Psim1 = Psim1 * F + fp; Psi = Psi * F
    h = Z2m1 * Phi * Psi + Phim1 * Psi + Psim1
    g1 = G('regd', e) - 0.5 * (G('d2lnPhiDiag', e) + bk.d2lnPsi(tuple(S), e))
    return g1, h

out = {'c0': c0, 'A2': A2, 'Lambda_zeta': LAM, 'Q_zeta': Q, 'rows': []}
print('\n[1] coefficients from the functions at small eps (polynomial fit on eps in [0.01, 0.2])')
print('  S          Lambda_S(fit h)  -g1(0)      formula    | Q_S = Im h\'(0)   formula   | g1 eps^2 coef  Re h eps^2 coef  Im h eps^3 coef')
e = np.linspace(0.01, 0.2, 40)
for S in [(), (2,), (3,), (5,), (7,), (11,), (13,), (2, 3), (3, 5), (2, 7, 13)]:
    g1, h = parts(S, e)
    # even fit for g1 and Re h, odd fit for Im h
    Ae = np.vstack([np.ones_like(e), e ** 2, e ** 4, e ** 6]).T; Ao = np.vstack([e, e ** 3, e ** 5, e ** 7]).T
    cg = np.linalg.lstsq(Ae, g1, rcond=None)[0]; cr = np.linalg.lstsq(Ae, h.real, rcond=None)[0]; ci = np.linalg.lstsq(Ao, h.imag, rcond=None)[0]
    row = dict(S=list(S), lam_fit=cr[0], lam_g1=-cg[0], lam_formula=lam_S(S), q_fit=ci[0], q_formula=q_S(S), g2=cg[1], h2=cr[1], h3=ci[1])
    out['rows'].append(row)
    print('  %-10s %10.6f  %10.6f  %10.6f   | %10.6f   %10.6f  | %+9.4f  %+9.4f  %+9.4f' % (str(S), cr[0], -cg[0], lam_S(S), ci[0], q_S(S), cg[1], cr[1], ci[1]))

print('\n[2] C(r) against the two-term expansion at large L (max over r in [0, 3] of the residual, scaled)')
r = np.linspace(0.0, 3.0, 301)
for S in [(), (2,), (5,)]:
    for L in (25., 50., 100., 200.):
        Cr = bk.C_of_r(G, L, tuple(S), r=r)
        two = -(4 * lam_S(S) / L ** 2) * np.sin(np.pi * r) ** 2 - (4 * math.pi * q_S(S) / L ** 3) * r * np.sin(2 * math.pi * r)
        one = -(4 * lam_S(S) / L ** 2) * np.sin(np.pi * r) ** 2
        print('  S=%-6s L=%5.0f  max|C - 1st|*L^3 = %.4f   max|C - 2 terms|*L^4 = %.4f' % (str(S), L, np.abs(Cr - one).max() * L ** 3, np.abs(Cr - two).max() * L ** 4))

print('\n[3] at the data heights: how far the expansion is from the exact C(r) (L = 9.4)')
L = 9.4
for S in [(), (2,), (3,), (5,)]:
    Cr = bk.C_of_r(G, L, tuple(S), r=r)
    two = -(4 * lam_S(S) / L ** 2) * np.sin(np.pi * r) ** 2 - (4 * math.pi * q_S(S) / L ** 3) * r * np.sin(2 * math.pi * r)
    one = -(4 * lam_S(S) / L ** 2) * np.sin(np.pi * r) ** 2
    i5 = 50; i10 = 100; i15 = 150
    print('  S=%-5s r=0.5: exact %+.5f  1-term %+.5f  2-term %+.5f | r=1: exact %+.5f | r=1.5: exact %+.5f 2-term %+.5f | max|exact|=%.4f max|exact-2term|=%.4f'
          % (str(S), Cr[i5], one[i5], two[i5], Cr[i10], Cr[i15], two[i15], np.abs(Cr).max(), np.abs(Cr - two).max()))
    # R2 <= 1 ?  (needed for a real determinantal kernel)
    rr = np.linspace(0.0, 8.0, 3201); Cf = bk.C_of_r(G, L, tuple(S), r=rr); R = 1 - np.sinc(rr) ** 2 + Cf
    print('        max R2 on [0,8] = %.5f at r = %.3f;  min(sinc^2 - C) = %+.2e' % (R.max(), rr[R.argmax()], (np.sinc(rr) ** 2 - Cf).min()))
json.dump(out, open('expansions.json', 'w'), indent=1)
