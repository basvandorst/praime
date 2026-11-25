"""check_degd.py: numerical checks of Theorems G1, G2, G3 and of the singular-series identities (degd_proof.md).
Usage: check_degd.py [sections], sections from A B C D E (default: all). Output to stdout.
Every assertion prints PASS/FAIL with the measured value and the tolerance."""
import sys, os, time, json, math, cmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import mpmath as mp
from degd_lib import *

UNI = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, UNI)
import ulib  # only build() and params(): the PARI objects of the zero lists

RESULTS = []
def check(name, ok, msg):
    RESULTS.append((name, bool(ok)))
    print(('PASS ' if ok else 'FAIL ') + name + ': ' + msg, flush=True)

_pari = None
def pari():
    global _pari
    if _pari is None:
        import cypari2
        _pari = cypari2.Pari()
        _pari.allocatemem(600 * 10 ** 6)
        _pari.set_real_precision(38)
    return _pari

def spec_of(key):
    for k, spec, T2, fam in ulib.CATALOGUE:
        if k == key: return spec
    raise KeyError(key)

def consts(key):
    return json.load(open(os.path.join(UNI, 'constants_%s.json' % key)))

def zeros(key):
    return np.load(os.path.join(UNI, 'zeros', 'zeros_%s.npz' % key))['g']

# ------------------------------------------------------------------ coefficient data
def lam_X5(N):
    """chi(n) for the PARI object X5 (Conrey 2 mod 5), from lfunan, periodic extension."""
    P = pari(); L = ulib.build(P, spec_of('X5'))
    an = P.lfunan(L, 5)
    per = [complex(ulib.to_c(an[i])) for i in range(5)]  # chi(1..5)
    lam = np.zeros(N + 1, complex)
    n = np.arange(1, N + 1)
    lam[1:] = np.array(per)[(n - 1) % 5]
    return lam

def lam_D12(N):
    P = pari()
    v = P('mfcoefs(mfDelta(),%d)' % N)
    a = np.array([float(v[i]) for i in range(N + 1)])
    lam = np.zeros(N + 1, complex)
    lam[1:] = a[1:] / np.arange(1, N + 1) ** 5.5
    return lam

def sym2_local_factory(ainv, P):
    """local lam(p^k) (k = 0..K) of Sym^2 E (analytic normalisation); bad p from PARI lfuneuler."""
    E = P.ellinit(ainv)
    L = P.lfunsympow(E, 2)
    N = int(P.ellglobalred(E)[0])
    cache = {}
    def poly(p):
        if p in cache: return cache[p]
        if N % p == 0:
            co = ulib.local_poly(P, L, p, 3.0)
        else:
            lp = int(P.ellap(E, p)) / math.sqrt(p)
            A = lp * lp - 1
            co = [1.0, -A, A, -1.0]
        cache[p] = co
        return co
    return poly, N

def lam_from_polys(Nmax, poly):
    def local(p):
        K = int(math.log(Nmax) / math.log(p)) + 1
        return series_inv(poly(p), K)
    return multiplicative_from_local(Nmax, local)

def lam_S(ainv, N, aps=None):
    P = pari()
    poly, cond = sym2_local_factory(ainv, P)
    return lam_from_polys(N, poly), poly

def Lambda_from_lam(lam, N):
    """Lambda_L(p^k) for all prime powers <= N from the Dirichlet coefficients (k lam(p^k) ln p = sum_j Lambda(p^j) lam(p^{k-j}))."""
    ns, ps, ks = prime_powers(N)
    Lam = {}
    out = np.zeros(len(ns), complex)
    for i, (n, p, k) in enumerate(zip(ns, ps, ks)):
        n, p, k = int(n), int(p), int(k)
        v = k * lam[n] * math.log(p)
        for j in range(1, k):
            v -= Lam[p ** j] * lam[p ** (k - j)]
        Lam[n] = v
        out[i] = v
    return ns, out

# ------------------------------------------------------------------ section A: phase, growth, Stirling
def X_inv(s, Q, mu):
    s = mp.mpc(s)
    v = mp.power(Q, s - 0.5)
    for m in mu:
        m = mp.mpc(m)
        GR = lambda z: mp.power(mp.pi, -z / 2) * mp.gamma(z / 2)
        v *= GR(s + m) / GR(1 - s + mp.conj(m))
    return v

def section_A():
    mp.mp.dps = 30
    sets = [('S11a', 121.0, [1.0, 1.0, 2.0]), ('G249', 249.0, [0.5, 0.5, 1.5, 1.5]),
            ('fake-complex', 37.0, [0.3 + 1.7j, 0.8, 1.2 - 0.4j])]
    # A1 phase identity X^{-1}(1/2+it) = e^{2 i theta(t)}
    err = 0
    for name, Q, mu in sets:
        for t in (0.7, 13.3, 57.1, 211.9):
            a = complex(X_inv(mp.mpc(0.5, t), Q, mu))
            b = cmath.exp(2j * float(theta(np.array([t]), Q, mu)[0]))
            err = max(err, abs(a - b))
    check('A1 phase X^-1(1/2+it)=e^{2i theta} (real and complex mu)', err < 1e-12, 'max err %.2e' % err)
    # A1b functional equation L(s) = eps X(s) Ltilde(1-s) with PARI, X5 (non-self-dual) and S11a
    P = pari()
    errs = []
    for key, dkey in (('X5', 'X5d'), ('S11a', 'S11a')):
        L = ulib.build(P, spec_of(key)); Ld = ulib.build(P, spec_of(dkey))
        N, k, mu = ulib.params(P, L)
        eps = complex(*consts(key)['eps'])
        sh = (k - 1) / 2
        for s in (complex(0.3, 7.1), complex(0.81, 23.4), complex(-0.2, 3.3)):
            Ls = complex(P.lfun(L, P(str(s.real + sh) + '+I*' + str(s.imag))))
            Lt = complex(P.lfun(Ld, P(str(1 - s.real + sh) + '+I*' + str(-s.imag))))
            Xs = 1 / complex(X_inv(mp.mpc(s.real, s.imag), N, mu))
            errs.append(abs(Ls - eps * Xs * Lt) / max(1e-300, abs(Ls)))
        # wrong root number fails
        wrong = abs(Ls - (-eps) * Xs * Lt) / abs(Ls)
        check('A1b FE with eps, %s; flipped eps fails' % key, max(errs) < 1e-12 and wrong > 0.5,
              'rel err %.2e, flipped %.2f' % (max(errs), wrong))
    # A2 growth |X^{-1}(sigma+it)| = (Q prod |t+Im mu_j|/2pi)^{sigma-1/2} (1 + O(1/t))
    rows = []
    for name, Q, mu in sets:
        for sig in (-0.25, 1.25):
            r = []
            for t in (50.0, 500.0, 5000.0):
                a = abs(complex(X_inv(mp.mpc(sig, t), Q, mu)))
                b = (Q * np.prod([abs(t + complex(m).imag) / (2 * math.pi) for m in mu])) ** (sig - 0.5)
                r.append(abs(a / b - 1) * t)
            rows.append(max(r))
    check('A2 growth of |X^-1|: t*|ratio-1| bounded (t = 50, 500, 5000)', max(rows) < 20, 'max t|ratio-1| = %.2f' % max(rows))
    # A3 Stirling: 2theta = t ln(Q(t/2pi)^d) - d t + b ln(t/2pi) + phi0 + c1/t + O(t^-2) (O(t^-3) for real mu)
    out = []
    for name, Q, mu in sets:
        d = len(mu); mu = [complex(m) for m in mu]
        b = sum(m.imag for m in mu)
        phi0 = math.pi / 2 * sum(m.real for m in mu) - d * math.pi / 4
        c1 = d / 24 + 0.5 * sum((m - m * m).real for m in mu)
        real = all(abs(m.imag) < 1e-15 for m in mu)
        mp.mp.dps = 40
        vals = []
        for t in (100.0, 400.0, 1600.0):
            th = mp.mpf(t) / 2 * mp.log(Q)
            for m in mu:
                th += mp.im(mp.loggamma((mp.mpf(0.5) + mp.mpc(m.real, m.imag) + mp.mpc(0, t)) / 2)) - (t + m.imag) / 2 * mp.log(mp.pi)
            approx = t * math.log(Q * (t / (2 * math.pi)) ** d) - d * t + b * math.log(t / (2 * math.pi)) + phi0 + c1 / t
            vals.append(float(2 * th - approx) * (t ** 3 if real else t ** 2))
        out.append((name, vals))
        print('   A3 %s: residual x t^%d at t = 100, 400, 1600: %s' % (name, 3 if real else 2, ', '.join('%.4f' % v for v in vals)))
    ok = all(max(abs(v) for v in vals) < 50 and abs(vals[-1] - vals[-2]) < 0.2 * max(1, abs(vals[-1])) for _, vals in out)
    check('A3 Stirling expansion with phi0 = (pi/2)sum Re mu - d pi/4, c1 = d/24 + (1/2)sum Re(mu - mu^2), b ln(t/2pi)', ok,
          'scaled residuals converge (see lines above)')

# ------------------------------------------------------------------ section B: stationary phase for one n (G2 constant)
def section_B():
    sets = [('S11a', 121.0, [1.0, 1.0, 2.0]), ('G249', 249.0, [0.5, 0.5, 1.5, 1.5]),
            ('D12', 1.0, [5.5, 6.5]), ('X5', 5.0, [1.0]), ('fake-complex', 37.0, [0.3 + 1.7j, 0.8, 1.2 - 0.4j]),
            ('fake-complex-b0', 37.0, [0.4 + 3j, 0.4 - 3j, 1.0])]
    allok = True
    for name, Q, mu in sets:
        d = len(mu); muc = [complex(m) for m in mu]
        b = sum(m.imag for m in muc)
        phi0 = math.pi / 2 * sum(m.real for m in muc) - d * math.pi / 4
        c1 = d / 24 + 0.5 * sum((m - m * m).real for m in muc)
        kap = c1 - 1 / (24 * d)
        A = cmath.exp(1j * (phi0 + math.pi / 4)) / math.sqrt(d)
        e0, e1 = [], []
        for tc in (100.0, 200.0, 400.0, 800.0):
            # smooth wide window: box [tc/2, 3tc/2] with Gaussian edges of width tc/8 (H' , H'' tiny at tc)
            T1, T2, D = tc / 2, 3 * tc / 2, tc / 8
            tg, h = grid(T1, T2, D, h=0.02)
            f = H_box(tg, T1, T2, D) * np.exp(2j * theta(tg, Q, muc))
            n = Q * (tc / (2 * math.pi)) ** d * 1.0001  # stationary point near tc
            I = osc_sum(np.array([math.log(n)]), tg, f, h)[0]
            ts = 2 * math.pi * (n / Q) ** (1 / d)
            lead = 2 * math.pi * A * (n / Q) ** (1 / (2 * d)) * cmath.exp(-2j * math.pi * d * (n / Q) ** (1 / d)) \
                * (n / Q) ** (1j * b / d) * H_box(ts - b / d, T1, T2, D)
            first = lead * (1 + 1j * kap / ts)
            e0.append(abs(I / lead - 1)); e1.append(abs(I / first - 1))
        # leading: error ~ |kap|/t (complex mu: also O(1/t) from b); first order: ~ 1/t^2 (real mu)
        real = abs(b) < 1e-12
        r0 = [e * t for e, t in zip(e0, (100, 200, 400, 800))]
        r1 = [e * t * t for e, t in zip(e1, (100, 200, 400, 800))]
        print('   B %s: kappa = %.6f; |I/lead-1|*t = %s; |I/first-1|*t^2 = %s' % (
            name, kap, ', '.join('%.4f' % x for x in r0), ', '.join('%.3f' % x for x in r1)))
        if real:
            ok = abs(r0[-1] - abs(kap)) < 0.05 * max(abs(kap), 0.02) + 2e-3 and max(r1) < 20 * (1 + kap * kap)
        else:
            ok = max(r0) < 10 and r0[-1] < 1.2 * r0[0] + 1e-3
        allok &= ok
        check('B leading constant e^{i(phi0+pi/4)}/sqrt(d) and first-order factor (1 + i kappa/t), %s' % name, ok,
              'leading rel err %.2e at t=800, first-order %.2e' % (e0[-1], e1[-1]))

# ------------------------------------------------------------------ section C: G1 exact identity on zeros
def g1_case(key, Q, mu, lam, Nmax, T1, T2, D, r=1.0, negzeros=None, conjtest=False, h=0.05):
    t0 = time.time()
    tg, h = grid(T1, T2, D, h=h)
    f = H_box(tg, T1, T2, D) * np.exp(2j * theta(tg, Q, mu))
    g = zeros(key)
    Z = np.sum(H_box(g, T1, T2, D) * np.exp(2j * theta(g, Q, mu)) * np.exp(-1j * g * math.log(r)))
    if negzeros is not None:
        gn = -negzeros
        Z += np.sum(H_box(gn, T1, T2, D) * np.exp(2j * theta(gn, Q, mu)) * np.exp(-1j * gn * math.log(r)))
    ns, Lam = Lambda_from_lam(lam, Nmax)
    I = osc_sum(np.log(ns) + math.log(r), tg, f, h)
    terms = Lam * ns ** -0.5 * I
    Pside = -np.sum(terms) / (2 * math.pi)
    tail = abs(np.sum(terms[ns > 0.8 * Nmax])) / (2 * math.pi)
    # E0 (form (ii)): (1/pi) int H theta' e^{2i theta} r^{-it} - (1/2pi) sum_m conj(Lambda(m)) m^{-1/2} int H e^{i(2theta + t ln(m/r))}
    bnd = h * np.sum(f * theta1(tg, Q, mu) * np.exp(-1j * tg * math.log(r))) / math.pi
    M = 2000
    sel = ns <= M
    Im = osc_sum(np.log(ns[sel]) - math.log(r), tg, f, h, sign=+1)
    E0 = bnd - np.sum(np.conj(Lam[sel]) * ns[sel] ** -0.5 * Im) / (2 * math.pi)
    res = abs(Z - Pside - E0)
    info = dict(Z=Z, P=Pside, E0=E0, res=res, tail=tail, nzeros=int(np.sum((g > T1 - 6 * D) & (g < T2 + 6 * D))),
                npp=len(ns), secs=time.time() - t0)
    if conjtest:
        Pc = -np.sum(np.conj(Lam) * ns ** -0.5 * I) / (2 * math.pi)
        info['conj_res'] = abs(Z - Pc - E0)
    return info

def section_C():
    P = pari()
    cases = []
    # X5: complex character mod 5 (degree 1, non-self-dual)
    N = 10000; lam = lam_X5(20000)
    neg = zeros('X5d')
    info = g1_case('X5', 5.0, [1.0], lam, N, 400.0, 440.0, 4.0, negzeros=neg, conjtest=True)
    cases.append(('X5 (d=1, chi mod 5, box [400,440], Delta=4)', info))
    info = g1_case('X5', 5.0, [1.0], lam, 20000, 400.0, 440.0, 4.0, r=1.5, negzeros=neg)
    cases.append(('X5 twisted r=3/2', info))
    # D12: Delta (degree 2, level 1)
    N = 30000; lam = lam_D12(N)
    zD = zeros('D12')
    info = g1_case('D12', 1.0, [5.5, 6.5], lam, N, 200.0, 260.0, 5.0, negzeros=zD)
    cases.append(('D12 (d=2, Delta, box [200,260], Delta=5)', info))
    # S11a: Sym^2 of 11a (degree 3)
    N = 1300000
    t0 = time.time(); lam, poly = lam_S([0, -1, 1, -10, -20], N); print('   (S11a coefficients to %d: %.1f s)' % (N, time.time() - t0))
    zS = zeros('S11a')
    info = g1_case('S11a', 121.0, [1.0, 1.0, 2.0], lam, N, 40.0, 50.0, 4.0, negzeros=zS)
    cases.append(('S11a (d=3, Sym^2 11a, box [40,50], Delta=4)', info))
    info = g1_case('S11a', 121.0, [1.0, 1.0, 2.0], lam, N, 40.0, 50.0, 4.0, r=2.0, negzeros=zS)
    cases.append(('S11a twisted r=2', info))
    for name, inf in cases:
        print('   C %s: zeros %s, primes %s, E0 %.2e, |zeros-primes-E0| = %.2e, tail(n>0.8N) %.1e, #zeros %d, #prime powers %d, %.1f s' % (
            name, fmtc(inf['Z']), fmtc(inf['P']), abs(inf['E0']), inf['res'], inf['tail'], inf['nzeros'], inf['npp'], inf['secs']))
        check('C G1 exact identity, ' + name, inf['res'] < 1e-8 and inf['tail'] < 1e-9, '|zeros - primes - E0| = %.2e' % inf['res'])
        if 'conj_res' in inf:
            check('C G1 conjugation test (conj Lambda in the prime sum fails), ' + name, inf['conj_res'] > 1e-2,
                  'with conj(Lambda): %.3f' % inf['conj_res'])
    # cross-check of the S11a coefficients against PARI lfunan of the object whose zeros were computed
    L = ulib.build(P, spec_of('S11a'))
    an = P.lfunan(L, 3000)
    lp = np.array([ulib.to_c(an[i]) for i in range(3000)]) / np.arange(1, 3001) ** 1.0
    check('C S11a coefficients = PARI lfunan (n <= 3000)', np.max(np.abs(lp - lam[1:3001])) < 1e-12,
          'max diff %.1e' % np.max(np.abs(lp - lam[1:3001])))
    L = ulib.build(P, spec_of('D12')); an = P.lfunan(L, 3000)
    lD = lam_D12(3000)
    lp = np.array([ulib.to_c(an[i]) for i in range(3000)]) / np.arange(1, 3001) ** 5.5
    check('C D12 coefficients = PARI lfunan (n <= 3000)', np.max(np.abs(lp - lD[1:])) < 1e-12, 'max diff %.1e' % np.max(np.abs(lp - lD[1:])))

def fmtc(z):
    return '%.10f%+.10fi' % (z.real, z.imag)

# ------------------------------------------------------------------ section D: G3 exact dual identity, and the T_q identity
def g3_lhs(Q, mu, lam, N, T1, T2, D, r, weights=None, h=0.05):
    tg, h = grid(T1, T2, D, h=h)
    f = H_box(tg, T1, T2, D) * np.exp(2j * theta(tg, Q, mu))
    n = np.arange(1, N + 1)
    I = osc_sum(np.log(n) + math.log(r), tg, f, h)
    c = lam[1:N + 1] if weights is None else lam[1:N + 1] * weights[1:N + 1]
    terms = c * n ** -0.5 * I
    return np.sum(terms), abs(np.sum(terms[n > 0.8 * N])), I

def g3_rhs(lam, eps, T1, T2, D, r, conj=True, M=None):
    M = M or int(200 * r) + 200
    m = np.arange(1, M + 1)
    c = np.conj(lam[1:M + 1]) if conj else lam[1:M + 1]
    return eps * np.sum(c * m ** -0.5 * Hhat_box(np.log(r / m), T1, T2, D))

def section_D():
    P = pari()
    T1, T2, D = 36.0, 44.0, 4.0
    N = 800000
    t0 = time.time(); lam, poly = lam_S([0, -1, 1, -10, -20], N); print('   (S11a coefficients to %d: %.1f s)' % (N, time.time() - t0))
    eps = complex(*consts('S11a')['eps'])
    Ivals = {}
    for r in (1.0, 2.0, 1.5):
        t0 = time.time()
        lhs, tail, I = g3_lhs(121.0, [1.0, 1.0, 2.0], lam, N, T1, T2, D, r)
        Ivals[r] = I
        rhs = g3_rhs(lam, eps, T1, T2, D, r)
        res = abs(lhs - rhs); flip = abs(lhs + rhs)
        print('   D S11a r=%g: lhs %s, eps*dual %s, |diff| %.2e, flipped eps |diff| %.3f, tail %.1e, %.1f s' % (
            r, fmtc(lhs), fmtc(rhs), res, flip, tail, time.time() - t0))
        check('D G3 exact identity S11a (d=3) r=%g, right root number' % r, res < 1e-8 and tail < 1e-9, '|lhs - eps dual| = %.2e' % res)
        check('D G3 S11a r=%g, flipped root number fails' % r, flip > 1e-2, '|lhs + eps dual| = %.3f' % flip)
    # G3' (log-weighted, exact, needs min Re mu > 0): sum_b lam(b) ln b b^{-1/2} I_b(a)
    #     = eps sum_m conj lam(m) m^{-1/2} int H(t) (2 theta'(t) - ln m) e^{-i t ln(a/m)} dt
    tg, h = grid(T1, T2, D, h=0.05)
    Ht = H_box(tg, T1, T2, D); th1 = 2 * theta1(tg, 121.0, [1.0, 1.0, 2.0])
    nn = np.arange(1, N + 1)
    for a in (2.0, 1.5):
        lhs = np.sum(lam[1:] * np.log(nn) * nn ** -0.5 * Ivals[a])
        M = 600
        rhs = 0j
        for m in range(1, M + 1):
            rhs += np.conj(lam[m]) * m ** -0.5 * h * np.sum(Ht * (th1 - math.log(m)) * np.exp(-1j * tg * math.log(a / m)))
        rhs *= eps
        reso = eps * np.conj(lam[int(a)]) * a ** -0.5 * h * np.sum(Ht * (th1 - math.log(a))) if a == int(a) else 0
        print('   D G3prime S11a a=%g: lhs %s, dual %s, |diff| %.2e; resonant term (m = a) %s' % (a, fmtc(lhs), fmtc(rhs), abs(lhs - rhs), fmtc(complex(reso))))
        check("D G3' exact log-weighted dual identity S11a a=%g" % a, abs(lhs - rhs) < 1e-8, '%.2e' % abs(lhs - rhs))
    # X5 (non-self-dual): conj on the dual side and the complex root number
    N5 = 20000; lam5 = lam_X5(N5); eps5 = complex(*consts('X5')['eps'])
    for (a, b, Dl) in ((400.0, 440.0, 4.0),):
        for r in (2.0, 1.0):
            lhs, tail, _ = g3_lhs(5.0, [1.0], lam5, N5, a, b, Dl, r)
            right = g3_rhs(lam5, eps5, a, b, Dl, r)
            noconj = g3_rhs(lam5, eps5, a, b, Dl, r, conj=False)
            epsbar = g3_rhs(lam5, eps5.conjugate(), a, b, Dl, r)
            print('   D X5 r=%g: lhs %s, eps*dual %s; |diff| %.2e; without conj %.3f; with conj(eps) %.3f' % (
                r, fmtc(lhs), fmtc(right), abs(lhs - right), abs(lhs - noconj), abs(lhs - epsbar)))
            check('D G3 exact identity X5 (d=1, complex eps) r=%g' % r, abs(lhs - right) < 1e-8 and tail < 1e-9, '%.2e' % abs(lhs - right))
            check('D G3 X5 r=%g: unconjugated dual coefficients or conj(eps) fail' % r,
                  abs(lhs - noconj) > 1e-2 and abs(lhs - epsbar) > 1e-2, '%.3f, %.3f' % (abs(lhs - noconj), abs(lhs - epsbar)))
    # T_q identity (singular series, exact form): sum_n c_q(n) lam(n) n^{-1/2} I_n
    #   = eps sum_{e|q} e mu(q/e) sum_f r_e(f) (ef)^{-1/2} sum_m conj lam(m) m^{-1/2} Hhat(ln(ef/m));
    # resonant part eps int H prod_{p|q} (p G(p) - 1), G(p) = 1 - F_p(1/p).
    I = Ivals[1]
    n = np.arange(0, N + 1)
    intH = T2 - T1
    for q in (2, 3, 6, 5):
        cq = ramanujan_c(q, n)
        lhs = np.sum(cq[1:] * lam[1:N + 1] * np.arange(1, N + 1) ** -0.5 * I)
        rhs, reso = tq_rhs(q, poly, lam, eps, T1, T2, D)
        print('   D T_q q=%d: lhs %s, exact dual %s, |diff| %.2e; resonant part eps*intH*prod(pG(p)-1) = %s' % (
            q, fmtc(lhs), fmtc(rhs), abs(lhs - rhs), fmtc(reso)))
        check('D T_q exact identity S11a q=%d' % q, abs(lhs - rhs) < 1e-8, '%.2e' % abs(lhs - rhs))

def ramanujan_c(q, n):
    out = np.zeros(len(n))
    for e in range(1, q + 1):
        if q % e == 0:
            mu = mobius(q // e)
            if mu: out += e * mu * (n % e == 0)
    return out

def mobius(n):
    r = 1; p = 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            if n % p == 0: return 0
            r = -r
        p += 1
    if n > 1: r = -r
    return r

def factor(n):
    f = {}; p = 2
    while p * p <= n:
        while n % p == 0: f[p] = f.get(p, 0) + 1; n //= p
        p += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def R_poly(Pp, lamp, a):
    """R_{p,a}(X) = P_p(X) * sum_j lam(p^{a+j}) X^j, a polynomial of degree <= d-1."""
    d = len(Pp) - 1
    out = []
    for k in range(d):
        out.append(sum(Pp[h] * lamp[a + k - h] for h in range(0, min(k, d) + 1)))
    return out

def tq_rhs(q, poly, lam, eps, T1, T2, D):
    total = 0j; reso_total = 0j
    for e in range(1, q + 1):
        if q % e: continue
        mu = mobius(q // e)
        if mu == 0: continue
        # r_e(f): product over p^a || e of R_{p,a}
        terms = {1: 1.0 + 0j}
        for p, a in factor(e).items():
            Pp = poly(p); lamp = series_inv(Pp, a + 10)
            R = R_poly(Pp, lamp, a)
            new = {}
            for f0, c0 in terms.items():
                for k, ck in enumerate(R):
                    new[f0 * p ** k] = new.get(f0 * p ** k, 0) + c0 * ck
            terms = new
        for f0, cf in terms.items():
            rr = e * f0
            total += e * mu * cf * rr ** -0.5 * g3_rhs(lam, eps, T1, T2, D, rr)
            reso_total += e * mu * cf * rr ** -1.0 * eps * np.conj(lam[rr]) * (T2 - T1)
    return total, reso_total

# ------------------------------------------------------------------ section E: singular series
def F_local(Pp, X, J=60):
    """F_p(X) = sum_h conj(lam(p^h)) mu(p^h) X^h, mu(p^h) = coefficients of P_p (1/L_p)."""
    lamp = series_inv(Pp, J)
    return sum(np.conj(lamp[h]) * Pp[h] * X ** h for h in range(len(Pp)))

def section_E():
    rng = np.random.default_rng(12345)
    # E1 Satake identities for random unitary parameters, d = 1..5
    e1 = e2 = e3 = e4 = 0.0
    for d in range(1, 6):
        for trial in range(200):
            al = np.exp(2j * math.pi * rng.random(d))
            Pp = np.poly(al)[::-1]  # coefficients of prod(1 - a X) in increasing powers? np.poly gives prod(x - a)
            # prod_j (1 - a_j X) = sum_k (-1)^k e_k X^k ; np.poly(al) = [1, -e1, e2, ...] for prod(x - a_j)
            Pp = np.poly(al)
            lamp = series_inv(Pp, 40)
            # (a) F_p(1) = 0
            F1 = sum(np.conj(lamp[h]) * Pp[h] for h in range(d + 1))
            e1 = max(e1, abs(F1))
            # (b) r_{p,1,k} = -mu(p^{k+1}) (R_{p,1}(X) = (1 - P_p(X))/X)
            R = R_poly(Pp, lamp, 1)
            e2 = max(e2, max(abs(R[k] + Pp[k + 1]) for k in range(d)))
            # (c) HL line factor A_p(v) = recipe factor S_p(v)/(1 - 1/p), v = -3..3, p = 7
            p = 7.0; X = 1 / p
            for v in range(-3, 4):
                # recipe: S_p(v) = sum_{m-h=v} conj lam(p^m) mu(p^h) p^{-(m+h)/2}
                S = sum(np.conj(lamp[h + v]) * Pp[h] * p ** (-(2 * h + v) / 2) for h in range(d + 1) if h + v >= 0)
                rec = S / (1 - X)
                # HL: q_p = 0 term + q_p = 1 term (mu(p)/phi(p)) * [ -(e=1 term) + p * sum_j r_{p,1,j} conj lam(p^{1+j+v}) p^{-1-j-v/2} ]
                t0 = (np.conj(lamp[v]) * p ** (-v / 2)) if v >= 0 else 0
                t1 = p * sum(R[j] * np.conj(lamp[1 + j + v]) * p ** (-1 - j - v / 2) for j in range(d) if 1 + j + v >= 0)
                hl = t0 + (-1 / (p - 1)) * (-t0 + t1)
                e3 = max(e3, abs(hl - rec))
            # (d) d = 3 closed form
            if d == 3:
                a = lamp[1]; om = np.prod(al); X = 0.37
                G = F_local(Pp, X) / (1 - X)
                Gc = 1 - (abs(a) ** 2 - 1) * X + (1 - 2 * abs(a) ** 2 + om * np.conj(a) ** 3) * X ** 2
                e4 = max(e4, abs(G - Gc))
    check('E1a F_p(1) = 0 for unitary Satake parameters (d = 1..5)', e1 < 1e-11, 'max |F_p(1)| = %.1e' % e1)
    check('E1b R_{p,1}(X) = (1 - P_p(X))/X, i.e. r_{p,1,k} = -mu(p^{k+1})', e2 < 1e-11, '%.1e' % e2)
    check('E1c Hardy-Littlewood line factor = recipe line factor rho_L(p^v)c_p/p^{|v|/2} (v = -3..3)', e3 < 1e-11, '%.1e' % e3)
    check('E1d d = 3: G_p(X) = 1 - (|a|^2-1)X + (1 - 2|a|^2 + omega conj(a)^3)X^2', e4 < 1e-11, '%.1e' % e4)
    # E2 partial Euler products prod_{p<=P} G_p(1/p) and partial q-sums, against the constants files
    P = pari()
    for key, ainv in (('S11a', [0, -1, 1, -10, -20]), ('S15a', [1, 1, 1, -10, -10]), ('D12', None), ('X5', None)):
        c = consts(key)['c']
        if key in ('S11a', 'S15a'):
            poly, cond = sym2_local_factory(ainv, P)
        elif key == 'D12':
            tauv = P('mfcoefs(mfDelta(),300000)')
            def poly(p, tauv=tauv):
                tp = float(tauv[p]) / p ** 5.5
                return [1.0, -tp, 1.0]
        else:
            ch = lam_X5(5)
            def poly(p):
                return [1.0, -ch[p % 5 if p % 5 else 5] if p != 5 else 0.0]
        Pmax = 1000000 if key != 'D12' else 300000
        pr = primes_upto(Pmax)
        logprod = 0j; marks = {}
        g = {}
        for p in pr:
            p = int(p)
            Pp = poly(p)
            Gp = F_local(Pp, 1.0 / p) / (1 - 1.0 / p)
            g[p] = 1 - Gp
            logprod += cmath.log(Gp)
            if p in (997, 9973, 99991, 299993, 999983):
                marks[p] = cmath.exp(logprod)
        line = ', '.join('P=%d: %.5f' % (k, v.real) for k, v in sorted(marks.items()))
        print('   E2 %s: c (constants file) = %.5f; prod_{p<=P} G_p(1/p): %s' % (key, c, line))
        last = marks[max(marks)]
        check('E2 Euler product of the exact local factor G_p(1/p) = F_p(1/p)/(1-1/p) tends to c_L, ' + key,
              abs(last - c) < 0.01 * max(1, abs(c)), 'last %.5f vs %.5f' % (last.real, c))
        # q-sum over squarefree q <= R of mu(q) prod_{p|q} g(p)   (Hardy-Littlewood form)
        R = 200000
        spf = np.zeros(R + 1, np.int64)
        for p in primes_upto(R)[::-1]:
            spf[p::p] = p
        val = np.zeros(R + 1, complex); val[1] = 1
        s = 1 + 0j; sums = {}
        for qq in range(2, R + 1):
            p = int(spf[qq]); m = qq // p
            if m % p == 0 or val[m] == 0:
                val[qq] = 0
            else:
                val[qq] = -val[m] * g[p]
            s += val[qq]
            if qq in (1000, 10000, 100000, 200000):
                sums[qq] = s
        print('   E2 %s: sum_{q<=R, squarefree} mu(q) prod_{p|q}(1 - G_p(1/p)): %s' % (
            key, ', '.join('R=%d: %.5f' % (k, v.real) for k, v in sorted(sums.items()))))

# ------------------------------------------------------------------ section F: G2 summed (illustration, smooth window)
def hann_conv(t, a, L, k=0, D=1.0):
    """k-th derivative (k <= 2) of (W * g_D)(t), W(u) = sin^2(pi (u - a)/L) on [a, a+L]: W^(k) * g_D
    (W and W' are continuous, so no boundary terms), by local quadrature over |u - t| <= 9D."""
    t = np.atleast_1d(np.asarray(t, float))
    z = np.linspace(-9 * D, 9 * D, 1801); dz = z[1] - z[0]
    gz = np.exp(-(z / D) ** 2) / (math.sqrt(math.pi) * D)
    out = np.empty(t.shape)
    for i0 in range(0, len(t), 2000):
        u = t[i0:i0 + 2000][:, None] - z[None, :]
        w = (u - a) / L
        inside = (u >= a) & (u <= a + L)
        if k == 0: W = np.sin(np.pi * w) ** 2
        elif k == 1: W = (np.pi / L) * np.sin(2 * np.pi * w)
        else: W = 2 * (np.pi / L) ** 2 * np.cos(2 * np.pi * w)
        W = np.where(inside, W, 0.0)
        out[i0:i0 + 2000] = dz * (W * gz[None, :]).sum(1)
    return out

def section_F():
    rows = []
    cfg = [('X5', 5.0, [1.0], None, 200.0, 400.0), ('D12', 1.0, [5.5, 6.5], None, 100.0, 300.0),
           ('S11a', 121.0, [1.0, 1.0, 2.0], [0, -1, 1, -10, -20], 20.0, 60.0)]
    for key, Q, mu, ainv, a, b in cfg:
        L = b - a; d = len(mu)
        phi0 = math.pi / 2 * sum(mu) - d * math.pi / 4
        c1 = d / 24 + 0.5 * sum(m - m * m for m in mu)
        kap = c1 - 1 / (24 * d)
        A = cmath.exp(1j * (phi0 + math.pi / 4)) / math.sqrt(d)
        # zero side
        g = zeros(key)
        Z = np.sum(hann_conv(g, a, L) * np.exp(2j * theta(g, Q, mu)))
        # prime powers with t_n* in [a - 8, b + 8]
        nmax = int(Q * ((b + 8) / (2 * math.pi)) ** d) + 2
        if key == 'X5': lam = lam_X5(nmax)
        elif key == 'D12': lam = lam_D12(nmax)
        else: lam, _ = lam_S(ainv, nmax)
        ns, Lam = Lambda_from_lam(lam, nmax)
        ts = 2 * math.pi * (ns / Q) ** (1 / d)
        sel = ts > a - 8
        ns, Lam, ts = ns[sel], Lam[sel], ts[sel]
        H0 = hann_conv(ts, a, L); H1 = hann_conv(ts, a, L, 1); H2 = hann_conv(ts, a, L, 2)
        ker = Lam * ns ** (-0.5 + 1 / (2 * d)) * np.exp(-2j * math.pi * d * (ns / Q) ** (1 / d))
        Plead = -A * Q ** (-1 / (2 * d)) * np.sum(ker * H0)
        Pfirst = -A * Q ** (-1 / (2 * d)) * np.sum(ker * (H0 * (1 + 1j * kap / ts) + 1j / (2 * d) * H1 + 1j * ts / (2 * d) * H2))
        intH = L / 2
        print('   F %s (d=%d, Hann [%g,%g] * g_1): zeros %s; leading %s (|diff| %.3f); first order %s (|diff| %.4f); kappa %.4f; int H/2pi = %.3f' % (
            key, d, a, b, fmtc(Z), fmtc(Plead), abs(Z - Plead), fmtc(Pfirst), abs(Z - Pfirst), kap, intH / (2 * math.pi)))
        rows.append((key, abs(Z - Plead), abs(Z - Pfirst)))
    check('F G2 illustration: first-order prime side closer to the zeros than the leading one (d = 2, 3)',
          all(r[2] < r[1] for r in rows if r[0] != 'X5'), ', '.join('%s %.3f -> %.4f' % r for r in rows))

def main():
    secs = sys.argv[1] if len(sys.argv) > 1 else 'ABCDEF'
    T0 = time.time()
    for S in secs:
        t0 = time.time()
        print('=== section %s' % S, flush=True)
        {'A': section_A, 'B': section_B, 'C': section_C, 'D': section_D, 'E': section_E, 'F': section_F}[S]()
        print('    (section %s: %.1f s)' % (S, time.time() - t0), flush=True)
    npass = sum(ok for _, ok in RESULTS)
    print('TOTAL %d/%d PASS, %.1f s' % (npass, len(RESULTS), time.time() - T0))

if __name__ == '__main__':
    main()
