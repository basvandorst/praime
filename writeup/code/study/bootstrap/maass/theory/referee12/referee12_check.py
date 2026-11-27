# referee12_check.py -- referee's own recomputation for exact12.md (self-contained; reads only cert12_d{0,8,18}.json
# for the polynomials B, and uses PARI only as a cross-check of tau(n) and of Xi).
# Route: Xi_Delta(t) = 4 sum_n tau(n) Re[(2 pi n)^{-(6+it)} Gamma(6+it, 2 pi n)] (q-expansion Mellin integral split at y=1),
# tau(n) from my own product q prod(1-q^n)^24; ghat = Xi^2 B on Gauss-Legendre nodes (unit panels on [0, 56]);
# g(x) = (1/pi) int_0^oo ghat cos(xt) dt; A^(k) = (2/pi) int_0^oo ghat [-log 2pi + Re psi(k/2+it)] dt.
# Delta_B(k) = A^(12) - A^(k) by the t-route and by the u-route 2 int_0^oo g(u)(e^{-ku/2}-e^{-6u})/(1-e^{-u}) du.
# Run: (ulimit -v 1200000; timeout 600 nice .venv/bin/python -I referee12_check.py)
import mpmath as mp, json, sys, time
from fractions import Fraction
mp.mp.dps = 60
HERE = '/home/bas/praime/writeup/code/study/bootstrap/maass/theory/'
T0 = time.time(); NCHK = [0, 0]
def check(name, cond, info=''):
    NCHK[0] += 1
    if not cond: NCHK[1] += 1
    print('%s  %s  %s' % ('ok  ' if cond else 'FAIL', name, info), flush=True)

# ------------------------------------------------------------------ tau, Xi
def tau_list(N):
    poly = [0]*N; poly[0] = 1
    for n in range(1, N):
        for _ in range(24):
            for k in range(N-1, n-1, -1): poly[k] -= poly[k-n]
    return [0] + poly
TAU = tau_list(402)
NTERMS, T = 26, 56
def xi(t, tau=TAU):
    w = mp.mpc(6, t); s = mp.mpf(0)
    for n in range(1, NTERMS+1):
        a = 2*mp.pi*n; s += tau[n]*mp.re(a**(-w)*mp.gammainc(w, a))
    return 4*s
def gl_nodes(n):
    xs, ws = [], []
    for i in range(1, n+1):
        x = mp.cos(mp.pi*(i-mp.mpf(1)/4)/(n+mp.mpf(1)/2))
        for it in range(102):
            p0, p1 = mp.mpf(1), x
            for k in range(2, n+1): p0, p1 = p1, ((2*k-1)*x*p1-(k-1)*p0)/k
            dp = n*(x*p1-p0)/(x*x-1)
            if it == 101 or abs(p1/dp) < mp.mpf(10)**(-mp.mp.dps+5): break
            x -= p1/dp
        xs.append(x); ws.append(2/((1-x*x)*dp*dp))
    return xs, ws
def nodes(NG):
    xs, ws = gl_nodes(NG); t, w, v = [], [], []
    for k in range(T):
        for x, wt in zip(xs, ws):
            tt = k + (x+1)/2; t.append(tt); w.append(wt/2); v.append(xi(tt))
    return t, w, v

# ------------------------------------------------------------------ g, A, Delta, P_p, Defect
def loadB(deg):
    d = json.load(open(HERE + 'cert12_d%d.json' % deg))
    return [mp.mpf(Fraction(s).numerator)/Fraction(s).denominator for s in d['b_exact']]
PR = [p for p in range(2, 401) if all(p % q for q in range(2, int(p**0.5)+1))]
class G:
    def __init__(self, nd, b):
        self.t, self.w, xv = nd
        self.gh = [x*x*mp.fsum(bk*t**(2*k) for k, bk in enumerate(b)) for x, t in zip(xv, self.t)]
    def g(self, x):
        x = mp.mpf(x); return mp.fsum(w*h*mp.cos(x*t) for w, h, t in zip(self.w, self.gh, self.t))/mp.pi
    def tint(self, f): return 2*mp.fsum(w*h*f(t) for w, h, t in zip(self.w, self.gh, self.t))/mp.pi
    def A(self, k): return self.tint(lambda t: -mp.log(2*mp.pi) + mp.re(mp.digamma(mp.mpf(k)/2 + 1j*t)))
    def Delta(self, k): return self.tint(lambda t: mp.re(mp.digamma(6+1j*t)) - mp.re(mp.digamma(mp.mpf(k)/2 + 1j*t)))
    def Aprime(self, s): return self.tint(lambda t: mp.re(mp.psi(1, mp.mpf(s)/2 + 1j*t)))/2
    def M(self, j): return self.tint(lambda t: t**(2*j))/2
def cheb2(x, m):
    c0, c1 = mp.mpf(2), 2*x
    if m == 0: return c0
    for _ in range(m-1): c0, c1 = c1, 2*x*c1 - c0
    return c1
def setup_P(Gb, NMAX=400):
    P = {}
    for p in PR:
        m, terms = 1, []
        while p**m <= NMAX:
            terms.append((m, mp.log(p)*mp.mpf(p)**(-mp.mpf(m)/2)*Gb.g(m*mp.log(p)))); m += 1
        P[p] = terms
    return P
def Pval(terms, x): return mp.fsum(c*cheb2(x, m) for m, c in terms)
def Pmin(terms, N=2000):
    xs = [-1 + mp.mpf(2)*i/N for i in range(N+1)]; vs = [Pval(terms, x) for x in xs]
    cands = [(vs[0], xs[0]), (vs[-1], xs[-1])]
    for i in range(1, N):
        if vs[i] <= vs[i-1] and vs[i] <= vs[i+1]:
            a, b = xs[i-1], xs[i+1]
            for _ in range(160):
                c1 = b - (b-a)/mp.phi; c2 = a + (b-a)/mp.phi
                if Pval(terms, c1) < Pval(terms, c2): b = c2
                else: a = c1
            xm = (a+b)/2; cands.append((Pval(terms, xm), xm))
    return min(cands)
def xp(p, tau=TAU): return mp.mpf(tau[p])/(2*mp.mpf(p)**(mp.mpf(11)/2))
def defect(P):
    sh = {}
    for p in PR:
        v, xm = Pmin(P[p]); sh[p] = (Pval(P[p], xp(p)) - v, xm)
    return mp.fsum(v[0] for v in sh.values()), sh
def LB_Delta(Ms, k):  # the note's cosine-truncation bound, best odd truncation
    a = mp.mpf(k)/2; acc = 2*Ms[0]*(mp.digamma(6) - mp.digamma(a)); best = None
    for j in range(1, len(Ms)):
        acc += 2*(-1)**j*Ms[j]*(mp.zeta(2*j+1, a) - mp.zeta(2*j+1, 6))
        if j % 2 == 1: best = acc if best is None else max(best, acc)
    return best
def LB_Aprime(Ms, s):
    a = mp.mpf(s)/2; acc = 0; best = None
    for j in range(len(Ms)):
        acc += (-1)**j*(2*j+1)*mp.zeta(2*j+2, a)*Ms[j]
        if j % 2 == 1: best = acc if best is None else max(best, acc)
    return best

# ================================================================== 1. normalisation
print('== 1. tau, Mellin normalisation, Xi against PARI', flush=True)
import cypari2
pari = cypari2.Pari(); pari.set_real_precision(60)
check('tau(n) own product = PARI ramanujantau, n <= 400', [int(pari('ramanujantau(%d)' % n)) for n in range(1, 401)] == TAU[1:401])
s = mp.mpf(10)
ds = 2*(2*mp.pi)**(-(s+5.5))*mp.gamma(s+5.5)*mp.fsum(TAU[n]*mp.mpf(n)**(-(s+5.5)) for n in range(1, 401))
ms = 2*mp.fsum(TAU[n]*((2*mp.pi*n)**(-(s+5.5))*mp.gammainc(s+5.5, 2*mp.pi*n) + (2*mp.pi*n)**(-(6.5-s))*mp.gammainc(6.5-s, 2*mp.pi*n)) for n in range(1, NTERMS+1))
check('Lambda(10): Gamma_C(10+11/2) L(10) = 2 int Delta(iy) y^{10+11/2} dy/y (split Mellin, uses Delta(i/y)=y^12 Delta(iy))', abs(ds-ms) < mp.mpf(10)**-25*abs(ds), mp.nstr(ds, 20))
for t in ['0', '5.5', '9.2', '9.2224', '20', '40']:   # strings: a binary float 9.2 differs from PARI's decimal 9.2
    mine = xi(mp.mpf(t)); pv = mp.mpf(str(pari.real(pari('lfunlambda(lfunetaquo(Mat([1,24])), 6+%s*I)' % t))).replace(' E', 'E'))
    check('Xi(%s) own Mellin = PARI lfunlambda(eta^24, 6+it), rel 1e-45' % t, abs(mine-pv) <= mp.mpf(10)**-45*max(abs(pv), mp.mpf(10)**-20), mp.nstr(mine, 20))
ND = {}
for NG in (24, 20):
    t1 = time.time(); ND[NG] = nodes(NG); print('   nodes NG=%d computed [%.0f s]' % (NG, time.time()-t1), flush=True)
check('Xi(56) negligible (|Xi(T)|^2 < 1e-60)', abs(ND[24][2][-1])**2 < mp.mpf(10)**-60, mp.nstr(ND[24][2][-1], 4))

# ================================================================== 2. B = 1, deg 8, deg 18
print('== 2. Lemma 1, Defect, kappa*, moment bound (own Fourier inversion)', flush=True)
NOTE = {0: ('0.3583', '0.3584', '0.1600', '4.086', 2), 8: ('5.402e-4', '5.403e-4', '0.1230', '8.785e-3', 11),
        18: ('2.248e-9', '2.249e-9', '0.1217', '3.694e-8', 5)}
RES = {}
for deg in (0, 8, 18):
    b = loadB(deg); Gb = G(ND[24], b); g0 = Gb.g(0)
    P = setup_P(Gb); A12 = Gb.A(12)
    res = A12 - 2*mp.fsum(Pval(P[p], xp(p)) for p in PR)
    check('deg %d: Lemma 1 residual/g(0) < 1e-40' % deg, abs(res/g0) < mp.mpf(10)**-40, 'A12/g0 = %s, residual/g0 = %s' % (mp.nstr(A12/g0, 12), mp.nstr(res/g0, 3)))
    bad = A12 - 2*mp.fsum(Pval(P[p], xp(p) + (mp.mpf('1e-3') if p == 2 else 0)) for p in PR)
    check('deg %d: negative control x_2 + 1e-3 breaks Lemma 1 (|residual|/g0 > 1e-9)' % deg, abs(bad/g0) > mp.mpf(10)**-9, mp.nstr(bad/g0, 3))
    tau2 = list(TAU); tau2[3] += 1
    bad = A12 - 2*mp.fsum(Pval(P[p], xp(p, tau2)) for p in PR)
    check('deg %d: negative control tau(3)+1 breaks Lemma 1' % deg, abs(bad/g0) > mp.mpf(10)**-12, mp.nstr(bad/g0, 3))
    D, sh = defect(P)
    lo, hi, apn, epn, car = NOTE[deg]
    check('deg %d: Defect/g0 in the note interval [%s, %s]' % (deg, lo, hi), mp.mpf(lo) <= D/g0 <= mp.mpf(hi), mp.nstr(D/g0, 10))
    top = max(sh.items(), key=lambda kv: kv[1][0])[0]
    check('deg %d: largest share p = %d' % (deg, car), top == car)
    check('deg %d: every Defect_p >= -1e-45 g0 (x_p in [-1,1])' % deg, min(v[0] for v in sh.values()) > -mp.mpf(10)**-45*g0)
    k = mp.findroot(lambda k: Gb.Delta(k) - 2*D, mp.mpf(12) - 2*D/Gb.Aprime(12))
    us, uw = gl_nodes(24); tot = 0
    for j in range(24):
        a0, a1 = mp.mpf(j)/4, mp.mpf(j+1)/4
        for x, w in zip(us, uw):
            u = a0 + (a1-a0)*(x+1)/2
            tot += w*(a1-a0)/2*Gb.g(u)*(mp.exp(-k*u/2) - mp.exp(-6*u))/(1 - mp.exp(-u))
    check('deg %d: Delta(kappa*) u-route = t-route (2 Defect), 1e-30' % deg, abs(2*tot/(2*D) - 1) < mp.mpf(10)**-30, mp.nstr(2*tot/(2*D)-1, 3))
    Ms = [Gb.M(j) for j in range(8)]
    lbA = LB_Aprime(Ms, 12); exA = Gb.Aprime(12)
    check("deg %d: cosine bound A'(12)/g0 reproduces the note's %s, and <= exact" % (deg, apn), mp.mpf(apn) < lbA/g0 < mp.mpf(apn) + mp.mpf('1e-4') and lbA <= exA, 'LB %s, exact %s' % (mp.nstr(lbA/g0, 8), mp.nstr(exA/g0, 8)))
    ke = 12 - mp.mpf(epn); lbD = LB_Delta(Ms, ke); exD = Gb.Delta(ke)
    check('deg %d: at kappa = 12 - %s: 2 Defect < cosine bound <= exact Delta' % (deg, epn), 2*D < lbD <= exD, 'LB/2D = %s, exact/2D = %s' % (mp.nstr(lbD/(2*D), 8), mp.nstr(exD/(2*D), 8)))
    lo_, hi_ = mp.mpf(0), mp.mpf(12)
    for _ in range(70):
        mid = (lo_+hi_)/2
        if LB_Delta(Ms, 12-mid) > 2*D: hi_ = mid
        else: lo_ = mid
    check('deg %d: eps from the cosine bound <= note %s; exact crossing eps = 12 - kappa* smaller' % (deg, epn), (12-k) < hi_ <= mp.mpf(epn), 'eps_LB = %s, eps_exact = %s' % (mp.nstr(hi_, 8), mp.nstr(12-k, 8)))
    check('deg %d: Delta(12 - eps_exact/2) < 2 Defect (kappa* is the crossing)' % deg, Gb.Delta(12 - (12-k)/2) < 2*D)
    amin = min(Gb.Aprime(s) for s in [mp.mpf('0.05'), 0.5, 1, 2, 3, mp.mpf('3.99'), 4])
    check("deg %d: A'_B(s) > 0 at s in {0.05, 0.5, 1, 2, 3, 3.99, 4} (outside Lemma 5's range; not needed)" % deg, amin > 0, mp.nstr(amin/g0, 4))
    RES[deg] = dict(g0=g0, D=D, kstar=k, epsLB=hi_, sh=sh)
    if deg == 8:
        v, xm = Pmin(P[2])
        check('deg 8: argmin P_2 = -0.268, Defect_2/g0 = 1.2e-7', abs(xm + mp.mpf('0.268')) < mp.mpf('5e-4') and abs((Pval(P[2], xp(2)) - v)/g0 - mp.mpf('1.18e-7')) < mp.mpf('1e-9'), '%s, %s' % (mp.nstr(xm, 6), mp.nstr((Pval(P[2], xp(2)) - v)/g0, 4)))
    if deg == 0:
        # closed form 8 e^{6w} sum (tau*tau)(N) K0(4 pi sqrt(N) e^{w/2}) (own Bessel series) vs Fourier inversion
        tt = [0]*401
        for a in range(1, 401):
            for c in range(1, 400//a + 1): tt[a*c] += TAU[a]*TAU[c]
        for n in (1, 2, 3):
            cf = 8*mp.mpf(n)**6*mp.fsum(tt[N]*mp.besselk(0, 4*mp.pi*mp.sqrt(N*n)) for N in range(1, 401))
            check('closed form phi*phi(log %d) = Fourier inversion of Xi^2, 1e-40 rel' % n, abs(cf - Gb.g(mp.log(n))) < mp.mpf(10)**-40*g0, mp.nstr(cf, 20))
        cfbad = 8*mp.fsum(TAU[N]*TAU[N]*mp.besselk(0, 4*mp.pi*mp.sqrt(N)) for N in range(1, 401))  # additive-diagonal misreading
        check('negative control: sum tau(N)^2 K0 (wrong convolution) differs', abs(cfbad - g0) > mp.mpf(10)**-8*g0, mp.nstr(cfbad, 10))
        us_ = [mp.mpf(i)/50 for i in range(301)]
        check('B=1: g > 0 on [0, 6] (grid 0.02)', all(Gb.g(u) > 0 for u in us_))
        check('B=1: 12 - eps_LB > 4 (union of the two ranges covers (0, 12 - eps))', 12 - hi_ > 4, mp.nstr(12-hi_, 8))

# other B, values the note did not use: Lemma 1 with B = (t^2 - 85)^2 (double zero near the first zero 9.22) and B = t^2
for name, b in [('(t^2-85)^2', [mp.mpf(7225), mp.mpf(-170), mp.mpf(1)]), ('t^2', [0, mp.mpf(1)])]:
    Gb = G(ND[24], b); g0 = Gb.g(0); P = setup_P(Gb)
    res = Gb.A(12) - 2*mp.fsum(Pval(P[p], xp(p)) for p in PR)
    check('Lemma 1 for B = %s: residual/|g(0)| < 1e-35' % name, abs(res/g0) < mp.mpf(10)**-35, mp.nstr(res/g0, 3))

# quadrature convergence: NG = 20 against NG = 24 at degree 18
b = loadB(18); G20 = G(ND[20], b); g0 = G20.g(0); D20, _ = defect(setup_P(G20))
check('NG 20 vs 24, deg 18: g(0) and Defect agree to 1e-40 rel', abs(g0/RES[18]['g0'] - 1) < mp.mpf(10)**-40 and abs(D20/RES[18]['D'] - 1) < mp.mpf(10)**-25, '%s, %s' % (mp.nstr(g0/RES[18]['g0']-1, 3), mp.nstr(D20/RES[18]['D']-1, 3)))

# ================================================================== 3. Lemma 5 and its extension to x >= 1/2
print('== 3. Re psi\'(x+iy)', flush=True)
ys = [mp.mpf(i)/20 for i in range(0, 4001)] + [mp.mpf(10)**e for e in range(3, 7)]
for x in ['0.55', '0.75', '1', '1.5', '2', '3']:
    m = min(mp.re(mp.psi(1, mp.mpf(x)+1j*y)) for y in ys)
    check("Re psi'(%s+iy) > 0 on y in [0,200] and 1e3..1e6" % x, m > 0, mp.nstr(m, 3))
for y in [0, 0.3, 1, 2.5]:
    check("reflection: Re psi'(1/2+iy) = pi^2/(2 cosh^2 pi y) at y=%s" % y, abs(mp.re(mp.psi(1, 0.5+1j*mp.mpf(y))) - mp.pi**2/(2*mp.cosh(mp.pi*y)**2)) < mp.mpf(10)**-50)
check("Re psi'(0.3+0.6i) < 0 and Re psi'(0.3+50i) < 0", mp.re(mp.psi(1, mp.mpf('0.3')+0.6j)) < 0 and mp.re(mp.psi(1, mp.mpf('0.3')+50j)) < 0)

# ================================================================== 4. an unconditionally feasible family (uniform d theta/pi)
print('== 4. uniform measure: <2 cos m theta> = 0, Q = A^(k)(g) >= 0 once psi(k/2) >= log 2 pi', flush=True)
ku = 2*mp.findroot(lambda x: mp.digamma(x) - mp.log(2*mp.pi), 6.8)
check('kappa_u = 2 psi^{-1}(log 2 pi) = 13.5531 (psi(kappa_u/2) - log 2pi = 0)', abs(ku - mp.mpf('13.5531')) < mp.mpf('1e-4') and abs(mp.digamma(ku/2) - mp.log(2*mp.pi)) < mp.mpf(10)**-50, mp.nstr(ku, 10))
check('uniform measure: (1/pi) int_0^pi 2 cos(m theta) d theta = 0, m = 1..6', all(abs(mp.quad(lambda th: 2*mp.cos(m*th), [0, mp.pi])/mp.pi) < mp.mpf(10)**-50 for m in range(1, 7)))
worst = min(mp.re(mp.digamma(mp.mpf(x)+1j*mp.mpf(y))) - mp.digamma(mp.mpf(x)) for x in ['0.3', '1', '3', '6.8'] for y in [mp.mpf(k)/4 for k in range(1, 400)])
check('Re psi(x+iy) >= psi(x) (sample; termwise (n+x)/((n+x)^2+y^2) <= 1/(n+x))', worst >= 0, mp.nstr(worst, 3))
check('uniform family not feasible at kappa = 12: Omega_12(0) = psi(6) - log 2pi < 0', mp.digamma(6) - mp.log(2*mp.pi) < 0, mp.nstr(mp.digamma(6) - mp.log(2*mp.pi), 6))

print('\n%d checks, %d failed   [%.0f s]' % (NCHK[0], NCHK[1], time.time()-T0))
