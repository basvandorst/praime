"""spectro_lib.py: primes, Hecke eigenvalues on all integers (segments), resonance weights rho, main terms.

Labels (degree 2, normalised so that |lambda(p)| <= 2):
  EC:<name>  elliptic curve (a-invariants in CURVES), lambda(n) = a_n / sqrt(n), multiplicative, from PARI ellap
  A<m>       Q(i): lambda(n) = sum_{a >= 1, b >= 0, a^2 + b^2 = n} cos(m atan2(b, a))   (m = 0 mod 4), N = 4, chi_{-4}
  C23        (r_{x^2+xy+6y^2}(n) - r_{2x^2+xy+3y^2}(n)) / 2, N = 23, chi_{-23}
Chirp: e(-2 sqrt(r n / N)); main term of the smooth integer sum
  sum_n lambda(n) n^{-1/4} W(sqrt n) e(-2 sqrt(r n/N)) ~ c_eff * (3/2) K1 r^{-1/4} int W(v) (1 - i kappa sqrt(N)/(2 pi v sqrt r)) dv
  with K1 = eps e^{-i(phi0 + pi/4)} (2 sqrt2/3) N^{-1/4} (the constant-1 law), c_eff from the resonance weights (spectro.md).
"""
import math, os, json
import numpy as np

TWO_PI = 2 * math.pi
HERE = os.path.dirname(os.path.abspath(__file__))
CURVES = {'11a1': [0, -1, 1, -10, -20], '15a1': [1, 1, 1, -10, -10], '14a1': [1, 0, 1, 4, -6], '37a1': [0, 0, 1, -1, 0]}

# ------------------------------------------------------------------ primes
def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0].astype(np.int64)

_SP = {}
def small_primes(r):
    if 'p' not in _SP or _SP['p'][-1] < r: _SP['p'] = primes_upto(max(r + 10, 100000))
    return _SP['p']

def is_prime_seg(L, U):
    """bool array s[i] = (L + i is prime), L <= n < U."""
    s = np.ones(U - L, bool)
    for n in range(L, min(U, 2)): s[n - L] = False
    r = math.isqrt(U) + 1
    for q in small_primes(r):
        q = int(q)
        if q > r: break
        st = max(q * q, -(-L // q) * q)
        if st < U: s[st - L::q] = False
    return s

# ------------------------------------------------------------------ elliptic curves: a_p table (cached), lambda on segments
_PARI = [None]
def pari():
    if _PARI[0] is None:
        import cypari2
        _PARI[0] = cypari2.Pari(); _PARI[0].allocatemem(400 * 10 ** 6)
    return _PARI[0]

def ap_table(name, X, cache=True):
    fn = os.path.join(HERE, 'aptab_%s_%d.npz' % (name, X))
    if cache and os.path.exists(fn):
        z = np.load(fn); return z['p'], z['a']
    P = pari(); ai = CURVES[name]
    ps, aa = [], []
    step = 2 * 10 ** 6
    for lo in range(0, X, step):
        hi = min(X, lo + step)
        v = P('my(E=ellinit(%s), Pr=primes([%d,%d])); [Vecsmall(Pr), Vecsmall(apply(p->ellap(E,p), Pr))]' % (ai, lo + 1, hi))
        ps.append(np.array(list(v[0]), np.int64)); aa.append(np.array(list(v[1]), np.int64))
    p = np.concatenate(ps); a = np.concatenate(aa)
    if cache: np.savez(fn, p=p, a=a)
    return p, a

def ap_at(name, plist):
    """a_p for an arbitrary list of primes (PARI), in chunks."""
    P = pari(); ai = CURVES[name]; out = []
    for i in range(0, len(plist), 100000):
        ch = plist[i:i + 100000]
        v = P('my(E=ellinit(%s)); Vecsmall(apply(p->ellap(E,p), %s))' % (ai, '[' + ','.join(map(str, map(int, ch))) + ']'))
        out.append(np.array(list(v), np.int64))
    return np.concatenate(out) if out else np.zeros(0, np.int64)

def conductor(name):
    return int(pari().ellglobalred(pari().ellinit(CURVES[name]))[0])

def lam_pow_table(lp, chi, kmax=40):
    """lambda(p^k), k = 0..kmax, from lambda(p) and chi(p) (chi = 0 at bad p): l_{k+1} = l l_k - chi l_{k-1}."""
    t = np.zeros(kmax + 1); t[0] = 1.0
    if kmax >= 1: t[1] = lp
    for k in range(1, kmax): t[k + 1] = lp * t[k] - chi * t[k - 1]
    return t

def lam_ec_seg(L, U, ptab, atab, N):
    """lambda(n) = a_n/sqrt(n) for L <= n < U (n >= 1)."""
    n = np.arange(L, U, dtype=np.int64); rem = n.copy(); lam = np.ones(U - L)
    if L == 0: lam[0] = 0.0; rem[0] = 1
    r = math.isqrt(U - 1)
    j = 0
    while j < len(ptab) and ptab[j] <= r:
        q = int(ptab[j]); a = int(atab[j]); bad = (N % q == 0)
        lp = a / math.sqrt(q)
        tab = lam_pow_table(lp, 0.0 if bad else 1.0, 64)
        st = (-L) % q
        idx = np.arange(st, U - L, q)
        if L == 0 and st == 0: idx = idx[1:]
        rr = rem[idx]; v = np.zeros(len(idx), np.int64)
        while True:
            m = (rr % q == 0)
            if not m.any(): break
            rr[m] //= q; v[m] += 1
        rem[idx] = rr; lam[idx] *= tab[v]
        j += 1
    m = rem > 1
    if m.any():
        k = np.searchsorted(ptab, rem[m]); assert np.all(ptab[k] == rem[m])
        lam[m] *= atab[k] / np.sqrt(rem[m].astype(np.float64))
    return lam

# ------------------------------------------------------------------ lattice labels
def isqrt_vec(x):
    x = np.asarray(x, np.int64); r = np.floor(np.sqrt(np.maximum(x, 0).astype(np.float64))).astype(np.int64)
    r[(r + 1) * (r + 1) <= x] += 1; r[r * r > x] -= 1
    return r

def _expand(lo, hi, step, other):
    cnt = np.maximum(0, (hi - lo) // step + 1); tot = int(cnt.sum())
    if tot == 0: return np.zeros(0, np.int64), np.zeros(0, np.int64)
    start = np.repeat(np.cumsum(cnt) - cnt, cnt)
    vals = np.repeat(lo, cnt) + step * (np.arange(tot, dtype=np.int64) - start)
    return vals, np.repeat(other, cnt)

def lam_gauss_seg(L, U, m, sub=4_000_000):
    """lambda(n) for A<m> (m = 0 mod 4), L <= n < U: sum over a >= 1, b >= 0 with a^2 + b^2 = n of cos(m atan2(b, a))."""
    lam = np.zeros(U - L)
    for L2 in range(L, U, sub):
        U2 = min(U, L2 + sub)
        bmax = math.isqrt(U2 - 1)
        b = np.arange(0, bmax + 1, dtype=np.int64)
        lo = np.maximum(1, isqrt_vec(L2 - b * b - 1) + 1) if L2 > 0 else np.ones_like(b)
        lo = np.where(L2 - b * b > 0, lo, 1)
        hi = isqrt_vec(U2 - 1 - b * b); hi = np.where(U2 - 1 - b * b >= 1, hi, 0)
        a, bb = _expand(lo, hi, 1, b)
        n = a * a + bb * bb; k = (n >= L2) & (n < U2); a, bb, n = a[k], bb[k], n[k]
        w = np.cos(m * np.arctan2(bb.astype(np.float64), a.astype(np.float64)))
        lam[L2 - L:U2 - L] += np.bincount(n - L2, weights=w, minlength=U2 - L2)
    return lam

def _form_count(L, U, A, B, C):
    """r_Q(n) for Q = A x^2 + B x y + C y^2 (positive definite), all (x, y) in Z^2, L <= n < U."""
    D = B * B - 4 * A * C
    ymax = math.isqrt(4 * A * (U - 1) // (-D)) + 1
    y = np.arange(-ymax, ymax + 1, dtype=np.int64)
    # 4A n = (2A x + B y)^2 - D y^2  =>  s = 2A x + B y in [sqrt(4A L + D y^2), sqrt(4A (U-1) + D y^2)]
    tU = 4 * A * (U - 1) + D * y * y; tL = 4 * A * L + D * y * y
    ok = tU >= 0; y = y[ok]; tU = tU[ok]; tL = tL[ok]
    sU = isqrt_vec(tU); sL = np.where(tL > 0, isqrt_vec(np.maximum(tL - 1, 0)) + 1, 0)
    cnt = np.zeros(U - L)
    for sg in (1, -1):
        # x = (sg*s - B y) / (2A), s in [sL, sU]
        if sg == 1: xlo = -(-(sL - B * y) // (2 * A)); xhi = (sU - B * y) // (2 * A)
        else: xlo = -(-(-sU - B * y) // (2 * A)); xhi = (-sL - B * y) // (2 * A)
        x, yy = _expand(xlo, xhi, 1, y)
        n = A * x * x + B * x * yy + C * yy * yy
        s = 2 * A * x + B * yy
        k = (n >= L) & (n < U) & ((s >= 0) if sg == 1 else (s < 0))
        cnt += np.bincount(n[k] - L, minlength=U - L)
    # s = 0 case (only n = 0 possible for definite forms except when ... ) ignored for n >= 1
    return cnt

def lam_c23_seg(L, U, sub=4_000_000):
    lam = np.zeros(U - L)
    for L2 in range(L, U, sub):
        U2 = min(U, L2 + sub)
        lam[L2 - L:U2 - L] = 0.5 * (_form_count(L2, U2, 1, 1, 6) - _form_count(L2, U2, 2, 1, 3))
    if L == 0: lam[0] = 0.0
    return lam

# ------------------------------------------------------------------ label metadata
def label_info(label):
    """N, eps, mu (analytic), chi(p) function, c (single-prime constant), rho_p function."""
    if label.startswith('EC:'):
        name = label[3:]
        cj = json.load(open(os.path.join(HERE, 'curves_c.json')))[name]
        return dict(N=cj['N'], eps=cj['eps'], mu=[0.5, 1.5], c=cj['c'], bad=cj['bad'])
    if label.startswith('A'):
        m = int(label[1:])
        cc = {4: 2.026446079853554, 8: 1.2660906768600704, 12: 1.2208738121184273, 16: 1.4705470091280404, 20: 0.8212970234372582}[m]
        return dict(N=4, eps=1, mu=[m / 2.0, m / 2.0 + 1], c=cc, bad=[2])
    if label == 'C23':
        return dict(N=23, eps=1, mu=[0.0, 1.0], c=1.8533516816725817, bad=[23])
    raise ValueError(label)

def K1_kappa(label):
    """K1 (the law with c = 1) and kappa of the first-order factor, degree 2."""
    I = label_info(label); d = 2; mu = I['mu']; Q = I['N']
    phi0 = (math.pi / 2) * sum(mu) - d * math.pi / 4
    K1 = complex(I['eps']) * complex(math.cos(phi0 + math.pi / 4), -math.sin(phi0 + math.pi / 4)) * (2 * math.sqrt(d) / (d + 1)) * Q ** (-1 / (2 * d))
    kappa = (d * d - 1) / (24 * d) + 0.5 * sum(m - m * m for m in mu)
    return K1, kappa

def hann(v, v0, v1):
    u = (v - v0) / (v1 - v0); return np.where((u > 0) & (u < 1), np.sin(np.pi * u) ** 2, 0.0)

def main_term_unit(label, v0, v1, r=1.0, npts=20001):
    """(3/2) K1 r^{-1/4} int W(v) (1 - i kappa sqrt(N) / (2 pi v sqrt r)) dv over the Hann window [v0, v1]."""
    K1, kappa = K1_kappa(label); N = label_info(label)['N']
    v = np.linspace(v0, v1, npts); f = hann(v, v0, v1) * (1 - 1j * kappa * math.sqrt(N) / (TWO_PI * v * math.sqrt(r)))
    h = v[1] - v[0]; s = f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum()
    return 1.5 * K1 * r ** -0.25 * s * h / 3

# ------------------------------------------------------------------ local data: lambda(p^k), chi(p), rho_p
def local_lams(label, pmax=10 ** 6):
    """dict p -> (lambda(p), chi(p)) for p <= pmax (chi = 0 at bad p; for A, C: nebentypus)."""
    ps = primes_upto(pmax); out = {}
    if label.startswith('EC:'):
        name = label[3:]; N = label_info(label)['N']
        p, a = ap_table(name, max(pmax, 2), cache=False) if pmax < 3 * 10 ** 7 else ap_table(name, pmax)
        for pp, aa in zip(p, a): out[int(pp)] = (aa / math.sqrt(pp), 0.0 if N % pp == 0 else 1.0)
        return out
    if label.startswith('A'):
        m = int(label[1:])
        lam = lam_gauss_seg(0, pmax + 1, m)
        for pp in ps:
            pp = int(pp); chi = 0.0 if pp == 2 else (1.0 if pp % 4 == 1 else -1.0)
            out[pp] = (lam[pp], chi)
        return out
    if label == 'C23':
        lam = lam_c23_seg(0, pmax + 1)
        leg = np.full(23, -1.0); leg[(np.arange(1, 23) ** 2) % 23] = 1.0; leg[0] = 0.0
        for pp in ps:
            pp = int(pp); out[pp] = (lam[pp], leg[pp % 23] if pp != 2 else 1.0)
        return out
    raise ValueError(label)

def rho_p(label, p, lp, chi):
    """resonance weight of the twist e(an/p), a != 0, relative to the main line (spectro.md, section 1):
    good p: (|lambda(p)|^2 - 1)/p = chi(p) lambda(p^2)/p; bad p: from c_p = 1 - rho_p (multiplicative: c_p = 1 + 1/p;
    ramified CM/dihedral (one Satake parameter of modulus 1): c_p = 1)."""
    if chi != 0.0: return (lp * lp - 1.0) / p
    if label.startswith('EC:'):
        return -1.0 / p if abs(lp) > 1e-12 else -1.0 / (p - 1)
    return 0.0

def R_class(label, l, lp, chi):
    """per-class resonance R(c), c mod l, normalised so that (1/l) sum_c R(c) = 1 - rho + rho = 1 ... i.e. the resonance of
    sum_{n = c mod l} lambda(n) chirp, relative to (1/l) x the full integer resonance.  Good l and multiplicative bad l:
    R(0) = 1 + (l - 1) rho, R(c != 0) = 1 - rho (the twisted resonances do not depend on a).  Ramified prime of a form with
    nebentypus chi_D (C23 at 23): R(0) = 1, R(c) = 1 + chi_D(c) (lambda vanishes on chi_D(n) = -1); A<m> at 2: R = 1.
    (checked on integers to 1e-4: out_validate_*.txt, class rows)."""
    if label == 'C23' and l == 23:
        leg = np.full(23, -1.0); leg[(np.arange(1, 23) ** 2) % 23] = 1.0; leg[0] = 0.0
        R = 1.0 + leg; R[0] = 1.0; return R
    rho = rho_p(label, l, lp, chi)
    R = np.full(l, 1.0 - rho); R[0] = 1.0 + (l - 1) * rho; return R

def local_factor(label, l, lp, chi, f):
    """L_l = (1/l) sum_c f(c) R(c) for a local weight table f (length l)."""
    return float(np.dot(f, R_class(label, l, lp, chi)) / l)
