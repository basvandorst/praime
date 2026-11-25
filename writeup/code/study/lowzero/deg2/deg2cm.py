"""Extension of deg2lib.Curve to curves with additive reduction (needed for CM curves, which are never semistable).
Additive p: W_p(k) = 0, O_p = 1 (the degree-1 missing orbit); the Sym^2 local factor is read from PARI (lfuneuler on
lfunsympow, arithmetic normalisation x = p^-s, turned analytic by x^j -> (x/p)^j). Multiplicative p as in deg2lib.
Pilot only (CM curves were not registered).
usage: deg2cm.py candidates OUT.json        (CM isogeny classes with 300 <= N <= 5000 from small equations)"""
import sys, os, math, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deg2lib as D
from scipy.interpolate import CubicSpline

CMJ = {0, 1728, -3375, 8000, -32768, 54000, 287496, -884736, -12288000, 16581375, -884736000, -147197952000,
       -262537412640768000}


def sym2_local_poly(pari, S, p):
    """coefficients c_0..c_d of the analytic Sym^2 local polynomial P(X) (L_p = 1/P(p^-s))"""
    f = pari.lfuneuler(S, p); den = pari.denominator(f); num = pari.numerator(f)
    assert int(pari.poldegree(num)) == 0
    n0 = float(pari.polcoef(num, 0)); c = [float(pari.polcoef(den, j)) / n0 for j in range(int(pari.poldegree(den)) + 1)]
    return [cj / p ** j for j, cj in enumerate(c)]


def local_add(p, poly, x):
    """additive p: Dcorr_p = 0 - zeta-part - Sym2-part; lnratio_p = ln (1-1/p)^-2 - ln[zeta_p zeta_p LS_p LS_p/LS_p(1)^2]"""
    lp = math.log(p); z = np.exp(-(1 + x) * lp); zm = np.exp(-(1 - x) * lp)
    P = lambda X: sum(c * X ** j for j, c in enumerate(poly))
    dP = lambda X: sum(j * c * X ** (j - 1) for j, c in enumerate(poly) if j)
    d2P = lambda X: sum(j * (j - 1) * c * X ** (j - 2) for j, c in enumerate(poly) if j > 1)
    # (L'/L)'(s) for L = 1/P(p^-s): log L = -log P(X), X = p^-s, d/ds = -ln p X d/dX
    # d/ds log L = ln p * X P'(X)/P(X);  d^2/ds^2 log L = -(ln p)^2 [X P'/P + X^2 (P''/P - (P'/P)^2)]
    X = z
    ls2 = -(lp ** 2) * (X * dP(X) / P(X) + X ** 2 * (d2P(X) / P(X) - (dP(X) / P(X)) ** 2))
    zeta2 = (lp ** 2) * z / (1 - z) ** 2                       # (zeta'/zeta)'_p
    Dcorr = -zeta2 - ls2
    lnT = -2 * math.log1p(-1 / p)
    lnden = -(np.log(1 - z) + np.log(1 - zm)) - (np.log(P(z)) + np.log(P(zm))) + 2 * math.log(P(1 / p))
    return Dcorr, lnT - lnden


class CurveAny(D.Curve):
    def __init__(self, pari, ainv, pmax=10 ** 5, h=D.HGRID, emax=D.EMAX):
        E = pari.ellinit(ainv); self.ainv = list(ainv)
        self.N = int(pari.ellglobalred(E)[0]); fa = pari.factor(self.N)
        self.bad = [int(p) for p in fa[0]]; self.add = [int(p) for p, e in zip(fa[0], fa[1]) if int(e) >= 2]
        an = pari.ellan(E, pmax); P = D.primes_upto(pmax)
        self.ap = {int(p): int(an[int(p) - 1]) for p in P}
        eps = D.eps_grid(h, emax); self.eps = eps; x = 1j * eps; x[0] = 1e-6j
        S = pari.lfunsympow(E, 2); Li = pari.lfuninit(S, [2, 0.05, emax + 0.5], 2)
        lnLS = np.empty(len(eps), complex); dd = np.empty(len(eps), complex)
        for i, e in enumerate(eps):
            s = pari('2+%.12f*I' % e)
            v0 = complex(pari.lfun(Li, s)); v1 = complex(pari.lfun(Li, s, 1)); v2 = complex(pari.lfun(Li, s, 2))
            lnLS[i] = np.log(v0); dd[i] = v2 / v0 - (v1 / v0) ** 2
        self.LS1 = float(np.exp(lnLS[0].real))
        Dc = np.zeros(len(eps), complex); lnA = np.zeros(len(eps), complex)
        good = np.array([int(p) for p in P if int(p) not in self.bad]); lam = np.array([self.ap[int(p)] / math.sqrt(p) for p in good])
        cuts = [0] + [int(np.searchsorted(good, b)) for b in (50, 1000)] + list(range(int(np.searchsorted(good, 1000)) + 400, len(good), 400)) + [len(good)]
        for k0, k1 in zip(cuts[:-1], cuts[1:]):
            if k1 <= k0: continue
            d, la, _, _, _ = D.local_good(good[k0:k1, None], lam[k0:k1, None], x[None, :]); Dc += d.sum(0); lnA += la.sum(0)
        self.sympoly = {}
        for p in self.bad:
            if p in self.add:
                poly = sym2_local_poly(pari, S, p); self.sympoly[p] = poly; d, la = local_add(p, poly, x)
            else:
                assert self.ap[p] in (1, -1); d, la, _, _, _ = D.local_mult(p, x)
            Dc += d; lnA += la
        lnSA = 2 * (lnLS.real - lnLS[0].real) + lnA; lnSA[0] = 0.0
        KSA = D.cexpm1(lnSA) / np.where(eps > 0, eps, 1.0) ** 2
        self.raw = dict(dd=dd, Dc=Dc, lnSA=lnSA, KSA=KSA)
        for k in ('dd', 'Dc', 'KSA'):
            v = self.raw[k]; c = np.polyfit(eps[1:5], v[1:5].real, 3); ci = np.polyfit(eps[1:5], v[1:5].imag, 3)
            v[0] = np.polyval(c, 0) + 1j * np.polyval(ci, 0)
        self.sp = {k: (CubicSpline(eps, v.real), CubicSpline(eps, v.imag)) for k, v in self.raw.items()}


if __name__ == '__main__' and sys.argv[1] == 'candidates':
    import cypari2
    pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6); seen = {}
    PK = [int(p) for p in pari.primes(46)]
    eqs = [[0, 0, 0, d, 0] for d in range(-200, 201) if d] + [[0, 0, 0, 0, d] for d in range(-400, 401) if d]
    eqs += [[a1, a2, a3, a4, a6] for a1 in (0, 1) for a2 in (-1, 0, 1) for a3 in (0, 1) for a4 in range(-30, 31) for a6 in range(-30, 31)]
    for ainv in eqs:
        try: E = pari.ellinit(ainv)
        except cypari2.PariError: continue
        if len(E) == 0 or E[11] == 0: continue
        j = E[12]
        if int(pari.denominator(j)) != 1 or int(j) not in CMJ: continue
        gr = pari.ellglobalred(E); N = int(gr[0])
        if not 300 <= N <= 5000: continue
        key = (N,) + tuple(int(pari.ellap(E, p)) for p in PK)
        if key in seen: continue
        Em = pari.ellchangecurve(E, gr[1]); seen[key] = ([int(Em[i]) for i in range(5)], int(j))
    rows = [{'ainv': a, 'N': k[0], 'j': j, 'id': ('N%d_%s' % (k[0], '_'.join(map(str, a)))).replace('-', 'm'),
             'zeros_ap_small': sum(1 for p, ap in zip(PK[:25], k[1:26]) if ap == 0 and k[0] % p)} for k, (a, j) in sorted(seen.items())]
    json.dump({'band': [20.0, 400.0], 'curves': rows}, open(sys.argv[2], 'w'), indent=0)
    print(len(rows), 'CM isogeny classes with 300 <= N <= 5000:', [(r['N'], r['j']) for r in rows][:60])
