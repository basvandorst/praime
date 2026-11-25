"""Hardy-Littlewood sums for the counts of pairs (n, n+h) in one admissible class b mod Q over [L0, L0+W),
with the split into the linear-response part (own code; definitions from STUDY.md rounds 31, 39, 42 and round 155).

Points of the class: n_k = b + Q k, k = 0..N-1, N = W // Q, profile m_k = L0 + Q/2 + Q k.
  r1(n) = (Q/phi(Q))/ln n          chance that n is prime, given that n is a unit mod Q
  r2(n) = c_Q(h)/(ln n ln(n+h))    chance that n and n+h are both prime, b admissible
  c_Q(h) = prod_{p|Q} (1-1/p)^-2 * prod_{p not | Q} (1 - nu_p(h)/p)/(1-1/p)^2,  nu_p(h) = #{0, h} mod p
  kappa(n) = r2/r1 = S(h) phi(Q)/(A_Q(h) ln(n+h))     (round 155, class form, in count units)
All singular series are class-conditioned: the factor at p | Q is 1.  For d = Q k:
  rho(k)  = S'({0,h,d,d+h})/S'({0,h})^2        four-point series   (round 42)
  s0(k)   = S'({0,d}),  s+(k) = S'({0,d+h}),  s-(k) = S'({0,d-h})     pair series of single primes
  tau1(k) = S'({0,h,d})/S'({0,h}),  tau2(k) = S'({0,h,d+h})/S'({0,h})  three-point series (twin and one prime)
With w(i) = sum_k r2_k r2_{k+i} (autocorrelation of the profile):
  Var X            = sum r2(1-r2)                         + 2 sum_{i>=1} w(i) (rho(i) - 1)
  Var(kappa Z)     = sum [r2^2/r1 (1-r1) + r2^2/r1'(1-r1') + 2 r2^2 (r2 - r1 r1')/(r1 r1')]
                                                           + 2 sum w(i) [2(s0-1) + (s+-1) + (s--1)]
  Cov(X, kappa Z)  = sum r2 [kappa (1-r1) + kappa'(1-r1')] + 2 sum w(i) [(tau1-1) + (tau2-1)]
Z_b is the class-b plus class-(b+h) prime count weighted by kappa; R = X - kappa Z.
Every local factor has mean exactly 1 over k mod p, so the infinite products are kept in full: the base sums run over
all p not dividing Q (to PMAX, tail by the exponential integral), and every special prime, also above PMAX, enters
with its ratio to the generic factor.
"""
import math
import numpy as np
from scipy.special import exp1

PMAX = 20_000_000


def small_primes(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0].astype(np.int64)


PR = small_primes(PMAX)
PRF = PR.astype(np.float64)
TAIL1 = float(exp1(math.log(PMAX)))           # sum_{p > PMAX} 1/p^2 ~ int_PMAX^inf dt/(t^2 ln t) = E1(ln PMAX)


def qfactors(Q):
    return [int(p) for p in PR[:12] if Q % p == 0]


def phi(Q):
    return int(round(Q * math.prod(1 - 1 / p for p in qfactors(Q))))


def admissible(Q, h):
    return [b for b in range(Q) if math.gcd(b, Q) == 1 and math.gcd(b + h, Q) == 1]


def S_pair(h):
    """Hardy-Littlewood constant S(h) of the pair {0, h} (h even)."""
    P = PRF[1:]
    nu = np.where(h % PR[1:] == 0, 1.0, 2.0)
    return 2.0 * math.exp(float(np.sum(np.log1p(-nu / P) - 2 * np.log1p(-1 / P))) - TAIL1)


def constants(Q, h):
    qf = qfactors(Q)
    keep = (Q % PR) != 0
    P = PRF[keep]; Pi = PR[keep]
    ph = (h % Pi) == 0                              # p | h, p not | Q
    nu2 = np.where(ph, 1.0, 2.0)
    logc = sum(-2 * math.log1p(-1 / p) for p in qf) + float(np.sum(np.log1p(-nu2 / P) - 2 * np.log1p(-1 / P))) - TAIL1
    # generic four-point factor: p not | h: (1-4/p)/(1-2/p)^2 ; p | h: (1-2/p)/(1-1/p)^2   (tail ~ -4/p^2)
    g4 = np.where(ph, np.log1p(-2 / P) - 2 * np.log1p(-1 / P), np.log1p(-4 / P) - 2 * np.log1p(-2 / P))
    base4 = float(np.sum(g4)) - 4 * TAIL1
    # pair series of single primes: generic (1-2/p)/(1-1/p)^2 = 1 - 1/(p-1)^2       (tail ~ -1/p^2)
    logC = float(np.sum(np.log1p(-1 / (P - 1) ** 2))) - TAIL1
    # three-point: generic p not | h: (1-3/p)/((1-2/p)(1-1/p)); p | h: (1-2/p)/(1-1/p)^2   (tail ~ -1/p^2)
    g3 = np.where(ph, np.log1p(-2 / P) - 2 * np.log1p(-1 / P), np.log1p(-3 / P) - np.log1p(-2 / P) - np.log1p(-1 / P))
    base3 = float(np.sum(g3)) - 2 * TAIL1                 # log t_p = -2/p^2 + O(p^-3) for p > PMAX
    A = len(admissible(Q, h))
    return dict(qf=qf, c=math.exp(logc), base4=base4, logC=logC, base3=base3, A=A, phi=phi(Q), S=S_pair(h))


def gen4(p, h):
    return math.log1p(-2 / p) - 2 * math.log1p(-1 / p) if h % p == 0 else math.log1p(-4 / p) - 2 * math.log1p(-2 / p)


def gen3(p, h):
    return math.log1p(-2 / p) - 2 * math.log1p(-1 / p) if h % p == 0 else math.log1p(-3 / p) - math.log1p(-2 / p) - math.log1p(-1 / p)


def support(Ap, Bp, k0, k1, Q, h, kind):
    """For v(k) = Ap*k + Bp, k in [k0, k1) (gcd(Ap, Bp) = 1), the additive contributions of the primes p not | Q
    dividing v(k), for the series named in kind: 'k' (v = k), 'plus'/'minus' (v = (Qk +- h)/gcd(Q, h)).
    Returns dict of arrays: rho, s, tau (log contributions relative to the generic factors)."""
    n = k1 - k0
    k = np.arange(k0, k1, dtype=np.int64)
    v = Ap * k + Bp
    vmax = int(v.max())
    out = {'rho': np.zeros(n), 's': np.zeros(n), 'tau': np.zeros(n)}
    Ls = np.zeros(n)
    lim = max(math.isqrt(vmax) + 1, 64)              # >= 64 > every prime of Q and of h
    for p in PR[PR <= lim].tolist():
        if Ap % p == 0:
            continue
        inv = pow(Ap, -1, p)
        r = (-Bp * inv) % p
        i0 = (r - k0) % p
        lp = math.log(p)
        if i0 < n:
            Ls[i0::p] += lp
            if Q % p:                                 # p not | Q: a special prime
                contrib(out, slice(i0, None, p), p, h, kind)
        pe = p * p
        while pe <= vmax:
            inv_e = pow(Ap, -1, pe)
            re = (-Bp * inv_e) % pe
            ie = (re - k0) % pe
            if ie < n:
                Ls[ie::pe] += lp
            pe *= p
    cof = np.rint(v / np.exp(Ls)).astype(np.int64)
    big = np.nonzero(cof > 1)[0]
    if big.size:
        P = cof[big].astype(np.float64)
        contrib_vec(out, big, P, h, kind, Q)
    return out


def _vals(p, h, kind):
    """(rho, s, tau) log-ratio contributions of one special prime p (not | Q) for the given series."""
    if kind == 'k':
        if h % p == 0:      # p | h and p | k: four-point nu 1, three-point nu3 1
            rho = -math.log1p(-1 / p) - gen4(p, h)
            tau = -math.log1p(-1 / p) - gen3(p, h)
        else:               # p | d: four-point nu 2 -> 1/(1-2/p); three-point nu3 2 -> 1/(1-1/p)
            rho = -math.log1p(-2 / p) - gen4(p, h)
            tau = -math.log1p(-1 / p) - gen3(p, h)
        s = math.log((p - 1) / (p - 2))
    else:                   # p | Qk +- h, p not | h (p | h primes never reach here: they divide k instead)
        rho = math.log1p(-3 / p) - 2 * math.log1p(-2 / p) - gen4(p, h)
        tau = -math.log1p(-1 / p) - gen3(p, h)
        s = math.log((p - 1) / (p - 2))
    return rho, s, tau


def contrib(out, sl, p, h, kind):
    if kind != 'k' and h % p == 0:
        return
    rho, s, tau = _vals(p, h, kind)
    out['rho'][sl] += rho; out['s'][sl] += s; out['tau'][sl] += tau


def contrib_vec(out, idx, P, h, kind, Q):
    # cofactor primes are > sqrt(vmax) >= 30 > h and do not divide Q or h
    if kind == 'k':
        rho = -np.log1p(-2 / P) - (np.log1p(-4 / P) - 2 * np.log1p(-2 / P))
        tau = -np.log1p(-1 / P) - (np.log1p(-3 / P) - np.log1p(-2 / P) - np.log1p(-1 / P))
    else:
        rho = np.log1p(-3 / P) - 2 * np.log1p(-2 / P) - (np.log1p(-4 / P) - 2 * np.log1p(-2 / P))
        tau = -np.log1p(-1 / P) - (np.log1p(-3 / P) - np.log1p(-2 / P) - np.log1p(-1 / P))
    s = np.log((P - 1) / (P - 2))
    out['rho'][idx] += rho; out['s'][idx] += s; out['tau'][idx] += tau


def profile_w(L0, Q, h, c, N, ngrid=4001):
    """w(i) = sum_{k=0}^{N-1-i} r2_k r2_{k+i} for i = 0..N-1, by Gauss-Legendre on a grid of i plus Euler-Maclaurin
    end terms, then cubic interpolation (validated against the exact FFT in check_w)."""
    r = lambda t: c / (np.log(L0 + Q / 2 + Q * t) * np.log(L0 + Q / 2 + Q * t + h))
    xg, wg = np.polynomial.legendre.leggauss(64)
    ig = np.linspace(0, N - 1, ngrid)
    M = (N - 1) - ig
    t = (M[:, None] / 2) * (xg[None, :] + 1)
    integ = np.sum(wg[None, :] * r(t) * r(t + ig[:, None]), axis=1) * M / 2
    fa = r(np.zeros_like(ig)) * r(ig); fb = r(M) * r(M + ig)
    # Euler-Maclaurin: sum_{k=0}^{M} f(k) = int + (f(0)+f(M))/2 + (f'(M)-f'(0))/12 ; f' is O(f/(x ln x)): dropped
    wgrid = integ + (fa + fb) / 2
    from scipy.interpolate import CubicSpline
    return CubicSpline(ig, wgrid)


def decompose(Q, h, L0, W, chunk=1 << 21, wfun=None):
    C = constants(Q, h)
    N = W // Q
    g = math.gcd(Q, h)
    Ap, Bp = Q // g, h // g
    if wfun is None:
        wfun = profile_w(L0, Q, h, C['c'], N)
    # diagonal sums (exact over the N points)
    mX = diagX = diagZ = diagXZ = sk = 0.0
    for j0 in range(0, N, chunk):
        m = L0 + Q / 2 + Q * np.arange(j0, min(N, j0 + chunk), dtype=np.float64)
        r1 = (Q / C['phi']) / np.log(m); r1b = (Q / C['phi']) / np.log(m + h)
        r2 = C['c'] / (np.log(m) * np.log(m + h))
        kap = r2 / r1; kapb = r2 / r1b
        mX += float(r2.sum())
        diagX += float(np.sum(r2 * (1 - r2)))
        diagZ += float(np.sum(kap ** 2 * r1 * (1 - r1) + kapb ** 2 * r1b * (1 - r1b) + 2 * kap * kapb * (r2 - r1 * r1b)))
        diagXZ += float(np.sum(r2 * (kap * (1 - r1) + kapb * (1 - r1b))))
        sk += float(np.sum(r2 * kap))
    kbar = sk / mX
    del m, r1, r1b, r2, kap, kapb
    off4 = offlin = offXZ = 0.0
    off_s0 = off_spm = 0.0
    sum_w = 0.0
    for k0 in range(1, N, chunk):
        k1 = min(N, k0 + chunk)
        w = wfun(np.arange(k0, k1, dtype=np.float64))
        rho, s0, sp, sm, t1, t2 = series(Q, h, k0, k1, C)
        off4 += 2 * float(np.sum(w * (rho - 1)))
        offlin += 2 * float(np.sum(w * (2 * (s0 - 1) + (sp - 1) + (sm - 1))))
        off_s0 += 2 * float(np.sum(w * 2 * (s0 - 1)))
        off_spm += 2 * float(np.sum(w * ((sp - 1) + (sm - 1))))
        offXZ += 2 * float(np.sum(w * ((t1 - 1) + (t2 - 1))))
        sum_w += float(np.sum(w))
    VX = diagX + off4; VZ = diagZ + offlin; CXZ = diagXZ + offXZ
    return dict(Q=Q, h=h, L0=L0, W=W, N=N, A=C['A'], c=C['c'], S=C['S'], kbar=kbar, mX=mX,
                vX=VX / mX, vZ=VZ / mX, cXZ=CXZ / mX, vR=(VX - 2 * CXZ + VZ) / mX, beta=CXZ / VZ,
                diagX=diagX / mX, off4=off4 / mX, diagZ=diagZ / mX, offlin=offlin / mX, off_same=off_s0 / mX,
                off_cross=off_spm / mX, diagXZ=diagXZ / mX, offXZ=offXZ / mX)


def series(Q, h, k0, k1, C=None):
    """rho, s0, s+, s-, tau1, tau2 for k in [k0, k1) (see the module docstring)."""
    if C is None:
        C = constants(Q, h)
    g = math.gcd(Q, h)
    Ap, Bp = Q // g, h // g
    Kk = support(1, 0, k0, k1, Q, h, 'k')
    Kp = support(Ap, Bp, k0, k1, Q, h, 'plus')
    Km = support(Ap, -Bp, k0, k1, Q, h, 'minus')
    rho = np.exp(C['base4'] + Kk['rho'] + Kp['rho'] + Km['rho'])
    s0 = np.exp(C['logC'] + Kk['s'])
    kk = np.arange(k0, k1, dtype=np.int64)
    extra = np.zeros(k1 - k0)                        # p | gcd(k, h), p not | Q, also divides d +- h
    for p in PR[PR <= h].tolist():
        if h % p == 0 and Q % p:
            extra[kk % p == 0] += math.log((p - 1) / (p - 2))
    sp = np.exp(C['logC'] + Kp['s'] + extra); sm = np.exp(C['logC'] + Km['s'] + extra)
    t1 = np.exp(C['base3'] + Kk['tau'] + Km['tau']); t2 = np.exp(C['base3'] + Kk['tau'] + Kp['tau'])
    return rho, s0, sp, sm, t1, t2


def direct(Q, h, d, kind):
    """Direct local products over all p <= PMAX (p not | Q) plus the tail, for one distance d (check only)."""
    keep = (Q % PR) != 0
    P = PR[keep]; Pf = P.astype(np.float64)
    r0 = 0; rh = h % P; rd = d % P; rdh = (d + h) % P
    if kind == 'rho':
        nu4 = 1 + (rh != 0) + ((rd != 0) & (rd != rh)) + ((rdh != 0) & (rdh != rh) & (rdh != rd))
        nu2 = 1 + (rh != 0)
        return math.exp(float(np.sum(np.log1p(-nu4 / Pf) - 2 * np.log1p(-nu2 / Pf))) - 4 * TAIL1)
    if kind == 'pair':      # S'({0, d})
        nu = 1 + (rd != 0)
        return math.exp(float(np.sum(np.log1p(-nu / Pf) - 2 * np.log1p(-1 / Pf))) - TAIL1)
    if kind == 'tau':       # S'({0,h,d})/S'({0,h})
        nu3 = 1 + (rh != 0) + ((rd != 0) & (rd != rh))
        nu2 = 1 + (rh != 0)
        return math.exp(float(np.sum(np.log1p(-nu3 / Pf) - np.log1p(-nu2 / Pf) - np.log1p(-1 / Pf))) - TAIL1)


def rho_direct(Q, h, k, P=20000):
    """Direct product for one k (check): prod over p < P, p not | Q, of (1 - nu4/p)/(1 - nu2/p)^2 (tail ignored)."""
    d = Q * k
    out = 0.0
    for p in PR[PR < P].tolist():
        if Q % p == 0:
            continue
        nu4 = len({0, h % p, d % p, (d + h) % p}); nu2 = len({0, h % p})
        out += math.log1p(-nu4 / p) - 2 * math.log1p(-nu2 / p)
    return math.exp(out)


def prime_class_v(Q, L0, W, chunk=1 << 21):
    """Var/mean of the prime count of one unit class mod Q over [L0, L0+W) from the Hardy-Littlewood pairs
    (class-conditioned pair series s0 at distances Qk), the statistic of round 36."""
    ph = phi(Q)
    N = W // Q
    mY = diag = 0.0
    for j0 in range(0, N, chunk):
        r1 = (Q / ph) / np.log(L0 + Q / 2 + Q * np.arange(j0, min(N, j0 + chunk), dtype=np.float64))
        mY += float(r1.sum()); diag += float(np.sum(r1 * (1 - r1)))
    del r1
    c1 = Q / ph
    xg, wg = np.polynomial.legendre.leggauss(64)
    r = lambda t: c1 / np.log(L0 + Q / 2 + Q * t)
    ig = np.linspace(0, N - 1, 4001); M = (N - 1) - ig
    t = (M[:, None] / 2) * (xg[None, :] + 1)
    wgr = np.sum(wg[None, :] * r(t) * r(t + ig[:, None]), axis=1) * M / 2 + (r(0 * ig) * r(ig) + r(M) * r(M + ig)) / 2
    from scipy.interpolate import CubicSpline
    wf = CubicSpline(ig, wgr)
    logC = constants(Q, 2)['logC']
    off = 0.0
    for k0 in range(1, N, chunk):
        k1 = min(N, k0 + chunk)
        Kk = support(1, 0, k0, k1, Q, 2, 'k')
        off += 2 * float(np.sum(wf(np.arange(k0, k1, dtype=np.float64)) * (np.exp(logC + Kk['s']) - 1)))
    return (diag + off) / mY, mY
