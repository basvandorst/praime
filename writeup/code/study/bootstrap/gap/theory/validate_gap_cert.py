# Ball-arithmetic verification of the round-175 dual certificates (class R: entire L-functions of real archimedean type).
# Own code; nothing from the repository is imported.  The certificate of record is the pair of exact rational files
# gap_cert_rational_T{9.23,9.225}.json (vectors s_i, r_j); on a first run they are built from cert_real_X32_T*.npz
# (eigen-factors of P_H and of Z P_G Z', Z rebuilt by the recipe of Gap.__init__, converted exactly, r_j projected exactly
# onto sum_k (-1)^k r_k = 0).  The npz matrices are only compared, never used in the certified numbers.
#
# Test function (support [-2a, 2a], a = log(32)/2, basis e^{i w_k u} on [-a, a], w_k = k pi/a, |k| <= 23):
#   ghat(t) = sum_i |shat_i(t)|^2 + (t^2 - T^2) sum_j |rhat_j(t)|^2,   G(u) = Re g(u) (even, real),
#   Ghat(t) = (ghat(t) + ghat(-t))/2 >= 0 for |t| >= T, by construction.
# Closed form on [0, 2a]:  G(u) = (2a - u) sum_{j=0}^{23} alpha_j cos(w_j u) + sum_{k=1}^{23} b_k sin(w_k u).
# Certified (arb): c0 = G(0) < 0; sup of F(kappa) = A_kappa(G) over kappa >= 0 (and over -1/2 <= kappa <= 0);
# max_theta P_p(theta), P_p(theta) = -sum_{p^m < 32} 2 log p p^(-m/2) G(m log p) cos(m theta), p <= 31, also with the
# Ramanujan radius 1 replaced by p or sqrt(p);  B = sup F + sum_p max P_p < 0.
# Float cross-checks: direct autocorrelation, Fourier route (digamma), cert.py numbers, and Weil's explicit formula
# against the zeros of L(s, Delta) and L(s, chi_-3) computed by PARI (counts checked against Riemann-von Mangoldt).
# Run: (ulimit -v 1500000; timeout 900 nice .venv/bin/python3 -I validate_gap_cert.py) > validate_gap_cert_out.txt
import sys, os, math, time, json, hashlib
from fractions import Fraction
from decimal import Decimal, Context, ROUND_CEILING, ROUND_FLOOR
import numpy as np
from flint import arb, acb, arb_poly, acb_poly, fmpq, ctx

ctx.prec = 128
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = '/home/bas/praime/writeup/data/study/gap'
NCHK = [0, 0]
def check(name, ok, info=''):
    NCHK[0] += 1
    if not ok: NCHK[1] += 1
    print('  [%s] %s %s' % ('ok' if ok else 'FAIL', name, info), flush=True)
def fup(x, d=6):
    """decimal string >= x (outward), d significant digits; x arb"""
    v = math.nextafter(float(x.upper()), math.inf)
    return str(Context(prec=d, rounding=ROUND_CEILING).plus(Decimal(v)))
def fdn(x, d=6):
    v = math.nextafter(float(x.lower()), -math.inf)
    return str(Context(prec=d, rounding=ROUND_FLOOR).plus(Decimal(v)))
def iv(x, d=8): return '[%s, %s]' % (fdn(x, d), fup(x, d))
def q_of(x): return fmpq(*Fraction(float(x)).as_integer_ratio())

PI = arb.pi(); I = acb(0, 1); EG = arb.const_euler()
A = arb(32).log()/2                     # a = log(32)/2, so the support of G is [-log 32, log 32]
PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
ZCACHE = {}

# ------------------------------------------------------------------------------------------------------------ A
class Cert:
    def __init__(self, Tstr):
        z = np.load('%s/cert_real_X32_T%s.npz' % (DATA, Tstr))
        self.Tstr = Tstr; self.file = {k: float(z[k]) for k in ('c0', 'supF', 'B')}
        PH, PG, N = z['PH'], z['PG'], int(z['N']); self.N = N; n = self.n = 2*N + 1
        assert float(z['X']) == 32.0 and abs(float(z['T']) - float(Tstr)) < 1e-12
        af = math.log(32.0)/2; ks = np.arange(-N, N + 1); self.ks = ks; self.af = af
        e = np.where(ks % 2 == 0, 1.0, -1.0)
        q, _ = np.linalg.qr(np.column_stack([e, np.eye(n)[:, :n - 1]])); Z = q[:, 1:]
        W = Z @ PG @ Z.T; W = (W + W.T)/2
        self.PHf, self.Wf = PH, W
        def factors(M, thr=1e-15):
            lam, V = np.linalg.eigh(M)
            return [np.sqrt(l)*V[:, i] for i, l in enumerate(lam) if l > thr], lam
        ef = [fmpq(1) if k % 2 == 0 else fmpq(-1) for k in ks]; self.ef = ef
        self.json = os.path.join(HERE, 'gap_cert_rational_T%s.json' % Tstr)
        if os.path.exists(self.json):      # the certificate of record (exact rationals); the npz is only compared
            d = json.load(open(self.json)); assert d['T'] == Tstr
            self.S = [[fmpq(*Fraction(x).as_integer_ratio()) for x in v] for v in d['s']]
            self.R = [[fmpq(*Fraction(x).as_integer_ratio()) for x in v] for v in d['r']]
            self.source = 'read from %s' % os.path.basename(self.json); self.dropped = None
        else:                              # first run: eigen-factors of P_H and Z P_G Z', converted exactly
            Sf, lamH = factors(PH); Rf, lamW = factors(W)
            self.dropped = (float(lamH[lamH <= 1e-15].sum()), float(lamW[lamW <= 1e-15].sum()))
            self.S = [[q_of(x) for x in v] for v in Sf]
            self.R = []
            for v in Rf:
                r = [q_of(x) for x in v]; t = sum((u*w for u, w in zip(ef, r)), fmpq(0))/n
                self.R.append([ri - t*ei for ri, ei in zip(r, ef)])
            blob = json.dumps({'T': Tstr, 's': [[str(x) for x in v] for v in self.S], 'r': [[str(x) for x in v] for v in self.R]})
            open(self.json, 'w').write(blob); self.source = 'built from the npz and written to %s' % os.path.basename(self.json)
        self.T = fmpq(*Fraction(Tstr).as_integer_ratio())
        self.exact()
    def exact(self):
        """alpha_j, b_k as exact rational parts:  alpha = aH + (pi/a)^2 aK - T^2 aR,
           b = (2a/pi)(bH - T^2 bR) + (2 pi/a) bK."""
        N, ks = self.N, self.ks
        def alpha(vecs):
            return [sum((v[N + j]**2 + (v[N - j]**2 if j else 0) for v in vecs), fmpq(0)) for j in range(N + 1)]
        def beta(vecs, wk):
            inv = {d: fmpq(1 if d % 2 == 0 else -1, d) for d in range(-2*N, 2*N + 1) if d}
            bt = []
            for ik, k in enumerate(ks):
                s = fmpq(0)
                for v in vecs:
                    inner = fmpq(0)
                    for ij, j in enumerate(ks):
                        if j != k: inner += v[ij]*inv[int(j - k)]*(int(j) if wk else 1)
                    s += v[ik]*inner*(int(k) if wk else 1)
                bt.append(s)
            return [bt[N + k] - bt[N - k] for k in range(1, N + 1)]
        self.aH, self.aR = alpha(self.S), alpha(self.R); self.aK = [self.aR[j]*j*j for j in range(N + 1)]
        self.bH, self.bR, self.bK = beta(self.S, False), beta(self.R, False), beta(self.R, True)
        Ta = arb(self.T); pa = PI/A
        self.alpha = [arb(self.aH[j]) + pa**2*arb(self.aK[j]) - Ta**2*arb(self.aR[j]) for j in range(N + 1)]
        self.b = [(2*A/PI)*(arb(self.bH[k]) - Ta**2*arb(self.bR[k])) + (2*PI/A)*arb(self.bK[k]) for k in range(N)]
        al = [acb(x) for x in self.alpha]; bl = [acb(0)] + [acb(x) for x in self.b]
        self.Pa, self.Pb = acb_poly(al), acb_poly(bl)
        qd = lambda c: acb_poly([sum(c[k + 1:], acb(0)) for k in range(len(c) - 1)])   # (P(x) - P(1))/(x - 1)
        self.Qa, self.Qb = qd(al), qd(bl)
        self.c0 = 2*A*sum(self.alpha, arb(0))
    # G and (G(u) - G(0))/u for complex balls u (entire functions; the second has no division by u)
    def G(self, u):
        u = acb(u); z = (I*PI*u/A).exp(); zi = 1/z
        return (2*A - u)*(self.Pa(z) + self.Pa(zi))/2 + (self.Pb(z) - self.Pb(zi))/(2*I)
    def Gd(self, u):
        z = (I*PI*u/A).exp(); zi = 1/z; sc = (PI*u/(2*A)).sinc()
        E1 = I*PI/A*(I*PI*u/(2*A)).exp()*sc; E2 = -I*PI/A*(-I*PI*u/(2*A)).exp()*sc     # (z-1)/u, (1/z-1)/u
        return 2*A*(self.Qa(z)*E1 + self.Qa(zi)*E2)/2 - (self.Pa(z) + self.Pa(zi))/2 + (self.Qb(z)*E1 - self.Qb(zi)*E2)/(2*I)
    def h(self, u):       # h(u) = 2u e^{-u/2} G(u)/(1 - e^{-2u}) = e^{u/2} G(u) u/sinh(u);  F'(kappa) = int e^{-kappa u} h
        u = acb(u); return (u/2).exp()*self.G(u)/(I*u).sinc()
    def F(self, kap, tol=2.0**-62):
        """A_kappa(G) = -(gamma + log pi) c0 - c0 log(1 - e^{-4a}) + int_0^{2a} [2e^{-2u}c0 - 2e^{-cu}G(u)]/(1-e^{-2u}) du,
           c = 1/2 + kappa, integrand rewritten without removable singularities (S(x) = sinh(x)/x = sinc(ix))."""
        c = arb(1)/2 + kap; c0 = self.c0
        def f(u, analytic):
            S = (I*u).sinc(); Sc = (I*(c - 2)*u/2).sinc()
            return u.exp()/S*(c0*(c - 2)*(-(c + 2)*u/2).exp()*Sc - (-c*u).exp()*self.Gd(u))
        J = acb.integral(f, 0, 2*A, rel_tol=tol, abs_tol=tol)
        assert J.imag.contains(0)
        return -(EG + PI.log())*c0 - c0*(1 - arb(2)**-10).log() + J.real
    def Pcoef(self, p):
        """w_m, m = 1..M: P_p(theta) = sum_m w_m cos(m theta), w_m = -2 log p p^(-m/2) G(m log p), p^m < 32"""
        out = []; m = 1
        while p**m < 32:
            out.append(-2*arb(p).log()*arb(p)**(-arb(m)/2)*self.G(m*arb(p).log()).real); m += 1
        return out

def cheb_to_mono(ws):
    """sum_m ws[m-1] T_m(x) -> monomial coefficients (arb)"""
    Tm = [[1], [0, 1]]
    while len(Tm) <= len(ws): t1, t0 = Tm[-1], Tm[-2]; Tm.append([0] + [2*c for c in t1]); Tm[-1] = [x - (t0[i] if i < len(t0) else 0) for i, x in enumerate(Tm[-1])]
    co = [arb(0)]*(len(ws) + 1)
    for m, w in enumerate(ws, 1):
        for i, c in enumerate(Tm[m]): co[i] += w*c
    return co

def poly_max(ws, R=arb(1), tol=1e-13):
    """certified [lower, upper] of max over |z| <= R of Re sum_m ws_m z^m = max_theta sum ws_m R^m cos(m theta)
       (maximum principle: Re of a polynomial is harmonic).  Branch and bound on x = cos(theta) in [-1, 1],
       mean-value form P(x) in P(c) + P'(X)(x - c)."""
    P = arb_poly(cheb_to_mono([w*R**m for m, w in enumerate(ws, 1)])); dP = P.derivative()
    lo = max((P(arb(float(x))).lower() for x in np.linspace(-1, 1, 257)), key=float)
    best_lo = lo; ub = None; stack = [(-1.0, 1.0)]; nint = 0
    while stack:
        l, r = stack.pop(); c = (l + r)/2; rad = (r - l)/2; nint += 1          # dyadic: exact
        Pc = P(arb(c)); U = (Pc + dP(arb(c, rad))*arb(0, rad)).upper()
        if float(Pc.lower()) > float(best_lo): best_lo = Pc.lower()
        if float(U) <= float(best_lo) + tol or rad < 2.0**-40:
            ub = U if ub is None or float(U) > float(ub) else ub
        else: stack += [(l, c), (c, r)]
    return best_lo, ub, nint

# ------------------------------------------------------------------------------------------------------------ run
def run(Tstr):
    t_start = time.time()
    print('\n=== certificate T = %s ===' % Tstr, flush=True)
    C = Cert(Tstr); N = C.N
    print('A. exact sum-of-squares form')
    print('  exact rational vectors %s' % C.source)
    print('  rank: %d vectors s_i (P_H), %d vectors r_j (Z P_G Z\')%s' % (len(C.S), len(C.R), '' if C.dropped is None else '; dropped eigenvalue mass %.1e, %.1e' % C.dropped))
    dens = [max(fmpq(x).q for v in C.S + C.R for x in v)]
    print('  largest denominator: 2^%d (r_j: dyadic float times rationals with denominator 47 from the projection)' % (int(dens[0]).bit_length() - 1))
    PHr = np.array([[float(sum((v[i]*v[j] for v in C.S), fmpq(0))) for j in range(C.n)] for i in range(C.n)])
    Wr = np.array([[float(sum((v[i]*v[j] for v in C.R), fmpq(0))) for j in range(C.n)] for i in range(C.n)])
    check('rational P_H = sum s_i s_i\' reproduces the file matrix', np.abs(PHr - C.PHf).max() < 1e-14, '(max diff %.1e)' % np.abs(PHr - C.PHf).max())
    check('rational W = sum r_j r_j\' reproduces Z P_G Z\'', np.abs(Wr - C.Wf).max() < 1e-14, '(max diff %.1e)' % np.abs(Wr - C.Wf).max())
    check('every r_j satisfies sum_k (-1)^k r_jk = 0 exactly (psi_j(+-a) = 0)', all(sum((u*w for u, w in zip(C.ef, r)), fmpq(0)) == 0 for r in C.R))
    check('T is the exact rational %s' % Tstr, C.T == fmpq(*Fraction(Tstr).as_integer_ratio()))
    print('  sha256 of %s: %s' % (os.path.basename(C.json), hashlib.sha256(open(C.json, 'rb').read()).hexdigest()[:16]))

    print('B. closed form of G against a direct autocorrelation (float, independent route)')
    af = C.af; om = C.ks*np.pi/af
    Sf = np.array([[float(x) for x in v] for v in C.S]); Rf = np.array([[float(x) for x in v] for v in C.R]); Tf = float(Tstr)
    def G_direct(u):
        x, w = np.polynomial.legendre.leggauss(400); lo, hi = -af, af - u
        xs = (x + 1)*(hi - lo)/2 + lo; ws = w*(hi - lo)/2
        E1 = np.exp(1j*np.outer(xs + u, om)); E0 = np.exp(1j*np.outer(xs, om))
        tot = 0
        for v in Sf: tot += np.sum(ws*(E1 @ v)*np.conj(E0 @ v))
        for v in Rf: tot += np.sum(ws*((E1 @ (1j*om*v))*np.conj(E0 @ (1j*om*v)) - Tf**2*(E1 @ v)*np.conj(E0 @ v)))
        return tot.real
    for u in (0.0, math.log(2), math.log(3), math.log(5), 2.5):
        g1 = float(C.G(arb(u)).real.mid()); g2 = G_direct(u)
        check('G(%.4f): closed form %.12e, direct %.12e' % (u, g1, g2), abs(g1 - g2) < 1e-12)
    for u in (0.3, 1.1, 2.9):
        check('(G(u)-G(0))/u, sinc form = difference quotient at u = %g' % u, ((C.Gd(acb(u)) - (C.G(arb(u)) - C.c0)/u).real).abs_upper() < 1e-30)
    check('G(2a) = 0 (support edge; exact by the closed form: (2a-2a)A + sum b_k sin(2 pi k))', C.G(2*A).real.abs_upper() < 1e-30)

    print('C. certified c0, prime values, prime maxima (Ramanujan |alpha| <= 1)')
    c0 = C.c0
    print('  c0 = G(0) = %s' % iv(c0, 12))
    check('c0 < 0', float(c0.upper()) < 0)
    check('c0 agrees with cert.py (%.6e)' % C.file['c0'], abs(float(c0.mid()) - C.file['c0']) < 1e-12)
    pmax = {}; pmax_lo = {}; wcoef = {}
    for p in PRIMES:
        ws = C.Pcoef(p); wcoef[p] = ws
        lo, up, nint = poly_max(ws)
        pmax[p], pmax_lo[p] = up, lo
        print('  p = %2d: m <= %d, G(log p) = %s;  max_theta P_p in [%s, %s] (%d boxes)' % (p, len(ws), iv(C.G(arb(p).log()).real, 6), fdn(lo, 8), fup(up, 8), nint))
    psum = sum(pmax.values(), arb(0)); psum_lo = sum(pmax_lo.values(), arb(0))
    print('  sum of prime maxima <= %s' % fup(psum, 8))

    print('D. certified sup of F(kappa) = A_kappa(G)')
    # D1: interval bounds of h on a partition of [0, 2a]
    K = 4096; hb = []; hu = []
    for i in range(K):
        mid = (i + 0.5)*2*af/K; rad = af/K*(1 + 1e-9) + 1e-15
        hv = C.h(acb(arb(mid, rad))).real
        hb.append(float(hv.abs_upper())); hu.append(math.nextafter(float(hv.upper()), math.inf) + 1e-15)
    hb = np.array(hb)*(1 + 1e-12); hu = np.array(hu); w = 2*af/K
    ulo = np.arange(K)*w; uhi = (np.arange(K) + 1)*w
    def M2(kl):   # sup_{kappa >= kl} |F''(kappa)| <= int_0^{2a} u e^{-kl u} |h(u)| du  (upper Riemann sum, float with margin)
        ex = np.exp(-kl*(ulo if kl >= 0 else uhi))
        return float(np.sum(hb*uhi*ex*w))*(1 + 1e-9) + 1e-15
    print('  partition of [0, 2a] into %d pieces: sup|h| <= %.4f; M2(0) = %.5f, M2(-1/2) = %.5f' % (K, hb.max(), M2(0.0), M2(-0.5)))
    # D2: tail.  h <= -eta on [0, u1], |h| <= M on [u1, 2a]  =>  F'(kappa) < 0 for kappa > log(1 + M/eta)/u1
    best = None
    for i1 in range(8, K//4):
        eta = -hu[:i1].max()
        if eta <= 0: break
        M = hb[i1:].max(); u1 = i1*w; K1 = math.log1p(M/eta)/u1
        if best is None or K1 < best[0]: best = (K1, i1, u1, eta, M)
    K1, i1, u1, eta, M = best
    K1g = math.ceil(K1*1.001 + 1e-9)                    # grid end (integer, >= K1)
    print('  tail: h <= -%.4e on [0, %.4f], |h| <= %.4f on [%.4f, 2a]  =>  F decreasing on [%.3f, oo); grid to %d' % (eta, u1, M, u1, K1, K1g))
    check('tail argument: F\'(kappa) < 0 for kappa >= K1 = %.3f' % K1, eta > 0 and K1 < K1g)
    # D3: grid scan with linear interpolation bound F <= max(F_l, F_r) + M2 h^2/8, refined near the maximum
    def scan(k_lo, k_hi, step, tol):
        Mtab = sorted({-0.5, 0.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0})
        M2l = lambda kl: M2(max([m for m in Mtab if m <= kl + 1e-12]))
        nodes = {}
        def Fat(k):
            if k not in nodes: nodes[k] = C.F(arb(k))
            return nodes[k]
        ks_ = [k_lo + i*step for i in range(int(round((k_hi - k_lo)/step)) + 1)]
        for k in ks_: Fat(k)
        L = max((nodes[k].lower() for k in ks_), key=float); argL = max(ks_, key=lambda k: float(nodes[k].lower()))
        stack = list(zip(ks_[:-1], ks_[1:])); Ubest = None; nref = 0
        while stack:
            kl, kr = stack.pop(); hh = kr - kl
            U = (Fat(kl).upper().max(Fat(kr).upper()) + arb(M2l(kl)*hh*hh/8)).upper()
            if float(U) <= float(L) + tol or hh < 2.0**-20:
                Ubest = U if Ubest is None or float(U) > float(Ubest) else Ubest
            else:
                km = (kl + kr)/2; nref += 1; Fm = Fat(km)
                if float(Fm.lower()) > float(L): L, argL = Fm.lower(), km
                stack += [(kl, km), (km, kr)]
        return L, argL, Ubest, len(nodes), nref
    t0 = time.time()
    L0, arg0, U0, nev, nref = scan(0.0, float(K1g), 1/32, 1e-10)
    print('  kappa in [0, %d]: %d integrals (%d refinements), %.0fs' % (K1g, nev, nref, time.time() - t0))
    print('  sup_{kappa >= 0} F in [%s, %s], attained near kappa = %.6f' % (fdn(L0, 10), fup(U0, 10), arg0))
    check('sup F agrees with cert.py\'s grid value %.6e within 1e-6' % C.file['supF'], abs(float(U0) - C.file['supF']) < 1e-6)
    t0 = time.time()
    Le, arge, Ue, nev_e, _ = scan(-0.5, 0.0, 1/32, 1e-8)
    print('  extension kappa in [-1/2, 0]: sup F <= %s (at %.4f; F(-1/2) = %s), %d integrals, %.0fs' % (fup(Ue, 8), arge, iv(C.F(arb(-0.5)), 8), nev_e, time.time() - t0))

    print('E. the certificate')
    Bup = (U0 + psum).upper(); Blo = L0 + psum_lo
    print('  B = sup F + sum_p max P_p  in  [%s, %s]' % (fdn(Blo, 8), fup(Bup, 8)))
    check('B < 0 and c0 < 0: THEOREM CERTIFIED for T0 = %s' % Tstr, float(Bup) < 0 and float(c0.upper()) < 0)
    check('B agrees with cert.py (%.6e) within 1e-6' % C.file['B'], abs(float(Bup) - C.file['B']) < 1e-6)
    Bext = (max(U0, Ue, key=float) + psum).upper()
    print('  with kappa allowed in (-1/2, oo): B <= %s' % fup(Bext, 8))
    print('  extension to real kappa > -1/2: %s' % ('certified (the sup over [-1/2, 0] is below the sup over [0, oo))' if float(Bext) < 0 and float(Ue) <= float(U0) else 'NOT certified'))
    # conductor slack: (log N / d) c0 + B < 0 still holds for log N / d > B/|c0|
    print('  conductor slack: the bound holds for log N/d > %s (N >= 1 is more than enough)' % fup(Bup/(-c0), 6))

    print('F. which Ramanujan hypotheses are used')
    one = arb(1)
    def Btot(Rdict):
        return (U0 + sum((poly_max(wcoef[p], Rdict.get(p, one))[1] for p in PRIMES), arb(0))).upper()
    # |alpha_j(p)| <= p holds for every L of the class (absolute convergence of the Euler product, Lemma 0)
    hyp = {}
    for S in ([2, 3, 5], [2, 5], [2, 3], [2], [3, 5], []):
        for other, oname in ((lambda p: arb(p), 'p'), (lambda p: arb(p).sqrt(), 'sqrt(p)')):
            Rd = {p: (one if p in S else other(p)) for p in PRIMES}
            bt = Btot(Rd); hyp[(tuple(S), oname)] = float(bt) < 0
            print('  Ramanujan only at %-10s |alpha_j(p)| <= %-7s at the other p <= 31:  B <= %s  -> %s' % (S, oname, fup(bt, 6), 'certified' if float(bt) < 0 else 'fails'))
    smallest = [S for S in ([2], [2, 3], [2, 5], [2, 3, 5]) if hyp[(tuple(S), 'p')]]
    check('Ramanujan at %s and the trivial bound |alpha| <= p at the other p <= 31 suffice' % (smallest[0] if smallest else None), bool(smallest))
    check('no Ramanujan at all (|alpha_j(p)| <= p for every p) is NOT enough with this certificate', not hyp[((), 'p')])
    lo_, hi_ = 1.0, 2.0
    for _ in range(30):
        mid = (lo_ + hi_)/2
        if float(Btot({p: arb(mid) for p in PRIMES})) < 0: lo_ = mid
        else: hi_ = mid
    print('  uniform relaxation: |alpha_j(p)| <= R for all p <= 31 is enough for R <= %.6f (certified at R = %.6f)' % (lo_, lo_))
    check('uniform R = %.6f certified' % lo_, float(Btot({p: arb(lo_) for p in PRIMES})) < 0)
    for p in (2, 3, 5):
        lo_, hi_ = 1.0, float(p)
        for _ in range(30):
            mid = (lo_ + hi_)/2; Rd = {q: (arb(mid) if q == p else one) for q in PRIMES}
            if float(Btot(Rd)) < 0: lo_ = mid
            else: hi_ = mid
        print('  only p = %d relaxed (|alpha| <= 1 at the others): |alpha_j(%d)| <= %.6f is enough' % (p, p, lo_))
    print('  time for this certificate: %.0fs' % (time.time() - t_start))
    return C, dict(c0=c0, U0=U0, L0=L0, arg0=arg0, pmax=pmax, Bup=Bup, Blo=Blo, Ue=Ue, wcoef=wcoef, hyp=smallest[0] if smallest else None)

# ------------------------------------------------------------------------------------------------------------ G, H
def float_routes(C, res):
    print('G. second route (float): the Fourier side  A_kappa(G) = (1/pi) int_0^oo Ghat_sym(t) [-log pi + Re psi((1/2+kappa+it)/2)] dt')
    from scipy.special import psi
    from scipy.integrate import quad
    af = C.af; om = C.ks*np.pi/af; Tf = float(C.Tstr)
    Sf = np.array([[float(x) for x in v] for v in C.S]); Rf = np.array([[float(x) for x in v] for v in C.R])
    sgn = np.where(C.ks % 2 == 0, 1.0, -1.0)
    def Ghat(t):          # ghat(t), not symmetrised
        X = np.sinc(af*(t[:, None] + om[None, :])/np.pi)*2*af
        return ((X @ Sf.T)**2).sum(1) + (t**2 - Tf**2)*((X @ Rf.T)**2).sum(1)
    Gs = lambda t: (Ghat(t) + Ghat(-t))/2
    def Qs(t):            # Ghat_sym(t) = 4 sin^2(a t) Qs(t)
        out = 0
        for s_ in (1, -1):
            V = sgn[None, :]/(s_*t[:, None] + om[None, :])
            out = out + ((V @ Sf.T)**2).sum(1) + (t**2 - Tf**2)*((V @ Rf.T)**2).sum(1)
        return out/2
    H = 1000.0; x, wq = np.polynomial.legendre.leggauss(10); edges = np.arange(0, H + 1e-9, 0.25)
    tt = (edges[:-1, None] + (x[None, :] + 1)*0.125).ravel(); wt = np.tile(wq*0.125, len(edges) - 1)
    Gv = np.concatenate([Gs(tt[i:i + 20000]) for i in range(0, len(tt), 20000)])
    check('Ghat_sym = 4 sin^2(at) Qs at t = 2500.3 (rational form used for tails)', abs(Gs(np.array([2500.3]))[0] - 4*np.sin(af*2500.3)**2*Qs(np.array([2500.3]))[0]) < 1e-15)
    Om = lambda k, t: -np.log(np.pi) + psi((0.5 + k + 1j*t)/2).real
    def qcos(f, lo, w):   # int_lo^oo f(t) cos(w t) dt: [lo, 50 lo] by QAWO (or plain), [50 lo, oo) by QAWF (or plain)
        if w == 0: return quad(f, lo, 50*lo, limit=2000, epsabs=1e-16)[0] + quad(f, 50*lo, np.inf, epsabs=1e-16, limit=200)[0]
        return quad(f, lo, 50*lo, weight='cos', wvar=w, limit=2000, epsabs=1e-16)[0] + quad(f, 50*lo, np.inf, weight='cos', wvar=w, limit=200)[0]
    def tailint(fun, u=0.0):   # (1/pi) int_H^oo 4 sin^2(at) Qs fun cos(tu) = (1/pi) int 2 Qs fun [cos(tu) - cos((2a+u)t)/2 - cos((2a-u)t)/2]
        f = lambda t: 2*Qs(np.array([t]))[0]*fun(t)
        return (qcos(f, H, u) - qcos(f, H, 2*af + u)/2 - qcos(f, H, abs(2*af - u))/2)/np.pi
    c0f = np.sum(wt*Gv)/np.pi + tailint(lambda t: 1.0)
    check('c0 by Fourier: %.13e vs certified %.13e' % (c0f, float(res['c0'].mid())), abs(c0f - float(res['c0'].mid())) < 1e-12)
    for u in (math.log(2), math.log(3)):
        gf = np.sum(wt*Gv*np.cos(tt*u))/np.pi + tailint(lambda t: 1.0, u)
        gc = float(C.G(arb(u)).real.mid())
        check('G(log %d) by Fourier inversion %.12e vs closed form %.12e' % (round(math.exp(u)), gf, gc), abs(gf - gc) < 1e-12)
    for k in (0.0, 1.0, 5.5, res['arg0'], 6.5, 20.0):
        Ff = np.sum(wt*Gv*Om(k, tt))/np.pi + tailint(lambda t: Om(k, t))
        Fc = float(C.F(arb(k)).mid())
        check('F(%.4f) by Fourier %.12e vs u-integral %.12e' % (k, Ff, Fc), abs(Ff - Fc) < 1e-11)
    # Ghat sign structure (by construction; sampled as a check of the Fourier formula)
    t1 = np.linspace(float(C.Tstr), 400, 400001); g1 = Gs(t1)
    check('Ghat_sym >= 0 on [T, 400] (sampled 4e5 points; holds by construction)', g1.min() > -1e-15, '(min %.2e)' % g1.min())
    t2 = np.linspace(0, float(C.Tstr), 20001); g2 = Gs(t2)
    print('  Ghat_sym on [0, T]: min %.4e at t = %.4f (negative: the test function is not positive-definite)' % (g2.min(), t2[g2.argmin()]))
    return Gs, Qs

def zeros_check(C, res, Gs, Qs):
    print('H. explicit formula against zeros computed by PARI (checks the normalisation of every term)')
    try:
        import cypari2
    except ImportError:
        print('  cypari2 missing: skipped'); return
    from scipy.special import psi, loggamma
    from scipy.integrate import quad
    af = C.af; Hz = 800.0
    tau = {2: -24, 3: 252, 5: 4830, 7: -16744, 11: 534612, 13: -577738, 17: -6905934, 19: 10661420, 23: 18643272, 29: 128406630, 31: -52843168}
    cases = [('Delta', 'lfunetaquo([1,24])', 1, [5.5, 6.5], lambda p: [math.acos(tau[p]/(2*p**5.5))]*1, 2),
             ('chi_-3', 'lfuncreate(-3)', 3, [1.0], None, 1)]
    def qinf(f, w=None):
        if w is None: return quad(f, Hz, 50*Hz, limit=2000, epsabs=1e-16)[0] + quad(f, 50*Hz, np.inf, epsabs=1e-16, limit=200)[0]
        return quad(f, Hz, 50*Hz, weight='cos', wvar=w, limit=2000, epsabs=1e-16)[0] + quad(f, 50*Hz, np.inf, weight='cos', wvar=w, limit=200)[0]
    for name, cmd, Ncond, kaps, ang, d in cases:
        if name not in ZCACHE:
            pari = cypari2.Pari(); pari.allocatemem(300*10**6); pari.set_real_precision(25)
            ZCACHE[name] = ([float(x) for x in pari('lfunzeros(%s, %g)' % (cmd, Hz))], [float(x) for x in pari('lfunzeros(%s, %g, 32)' % (cmd, Hz))])
        zz, zz32 = ZCACHE[name]
        # Riemann-von Mangoldt: N(T) = (1/pi)[(T/2) log N + sum_j Im log Gamma_R(1/2 + kappa_j + iT)] + 1/pi arg L + O(1/T)
        Nsm = lambda T: ((T/2)*math.log(Ncond) + sum(-(T/2)*math.log(math.pi) + loggamma((0.5 + k + 1j*T)/2).imag - loggamma((0.5 + k)/2).imag for k in kaps))/math.pi
        check('%s: %d zeros in (0, %g]; divz=32 gives %d (same values to 1e-12: %s); Riemann-von Mangoldt main term %.2f'
              % (name, len(zz), Hz, len(zz32), len(zz) == len(zz32) and max(abs(x - y) for x, y in zip(zz, zz32)) < 1e-12, Nsm(Hz)),
              len(zz) == len(zz32) and abs(len(zz) - Nsm(Hz)) < 3 and min(np.diff(zz)) > 0.05)
        zside = 2*np.sum(Gs(np.array(zz)))
        # tail: smooth density + the n = 32 resonance (log 32 = 2a matches sin^2(a t)); other n oscillate
        if name == 'Delta':
            th2 = math.acos(tau[2]/(2*2**5.5)); c32 = math.log(2)*2*math.cos(5*th2)
        else:
            c32 = math.log(2)*(-1)**5
        dens = lambda t: (math.log(Ncond) + sum(-math.log(math.pi) + psi((0.5 + k + 1j*t)/2).real for k in kaps))/(2*math.pi)
        f_s = lambda t: 2*Qs(np.array([t]))[0]*dens(t)
        sm = 2*(qinf(f_s) - qinf(f_s, 2*af))
        resn = 2*(c32/math.sqrt(32))/math.pi*qinf(lambda t: Qs(np.array([t]))[0])
        zside_t = zside + sm + resn
        # prime/arch side from the certified pieces
        Fk = sum((C.F(arb(k)) for k in kaps), arb(0))
        if name == 'Delta':
            pr = sum((sum((w*math.cos(m*ang(p)[0]) for m, w in enumerate(res['wcoef'][p], 1)), arb(0)) for p in PRIMES), arb(0))*d
        else:
            chi = lambda p: 0 if p == 3 else (1 if p % 3 == 1 else -1)
            pr = sum((sum((w*chi(p)**m for m, w in enumerate(res['wcoef'][p], 1)), arb(0)) for p in PRIMES), arb(0))
        Wv = float((res['c0']*arb(Ncond).log() + Fk + pr).mid())
        print('  %s: sum over zeros to %g: %.10e; tail (smooth %.3e, resonance n=32 %.3e) -> %.10e;  W(G) from c0, F, primes: %.10e;  diff %.2e'
              % (name, Hz, zside, sm, resn, zside_t, Wv, zside_t - Wv))
        check('%s: explicit formula holds within 1e-9 (|W| = %.2e; the tail correction is %.1e)' % (name, abs(Wv), sm + resn), abs(zside_t - Wv) < 1e-9)
        check('%s: without the n = 32 resonance term the agreement is worse' % name, abs(zside + sm - Wv) > abs(zside_t - Wv))
        if name == 'Delta':
            print('  Delta: first zero %.6f, 2 Ghat_sym(first zero) = %.6e; W(G)/2 = %.6e <= B: %s' % (zz[0], 2*Gs(np.array([zz[0]]))[0], Wv/2, Wv/2 <= float(res['Bup'])))
            check('Delta is consistent with the certificate (W/2 <= B; its zero 9.2224 lies inside (-T, T))', Wv/2 <= float(res['Bup']) and zz[0] < float(C.Tstr))

if __name__ == '__main__':
    t00 = time.time()
    print('Validation of the round-175 dual certificates (class R).  python-flint arb, %d bits.' % ctx.prec)
    out = {}
    for Tstr in ('9.23', '9.225'):
        C, res = run(Tstr); out[Tstr] = (C, res)
        Gs, Qs = float_routes(C, res)
        zeros_check(C, res, Gs, Qs)
    print('\nSUMMARY')
    for Tstr, (C, res) in out.items():
        print('  T0 = %s: c0 <= %s, sup F <= %s (kappa ~ %.4f), sum_p max P_p <= %s, B <= %s; Ramanujan needed only at %s'
              % (Tstr, fup(res['c0'], 6), fup(res['U0'], 8), res['arg0'], fup(sum(res['pmax'].values(), arb(0)), 8), fup(res['Bup'], 8), res['hyp']))
    print('total time %.0fs' % (time.time() - t00))
    print('%d checks, %d failed' % (NCHK[0], NCHK[1]))
