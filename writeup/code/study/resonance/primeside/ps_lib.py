"""ps_lib.py: primes and their Hecke eigenvalues for the prime-side test of the one-point law (fact 2055).

Everything here is elementary: a segmented sieve over odd numbers, lattice enumerations of binary quadratic forms,
and the eigenvalue lambda(p) of each case read off from the representation of p.

Cases (series): d = degree, Q = conductor, mu = archimedean parameters (analytic normalisation), eps = root number
(filled in by ps_constants.py from PARI), lam(p) = normalised Hecke eigenvalue at p.
  A<m>  (m = 4, 8, 12, 16, 20): Q(i), Hecke character (z/|z|)^m; p = a^2 + b^2: lam = 2 cos(m*atan(b/a)); inert 0; lam(2) = cos(m pi/4)*... (see SPECIAL)
  B<m>  (m = 6, 12): Q(sqrt-3), (z/|z|)^m; p = a^2 + 3b^2: lam = 2 cos(m*atan(sqrt3 b/a)); p = 2 mod 3: 0
  C23, C31: S3 Artin (x^3 - x - 1, x^3 + x - 1): lam = 2 (p = x^2+xy+6y^2 resp. x^2+xy+8y^2), -1 (other split), 0 (inert)
  C47a, C47b: D5 (class group of Q(sqrt-47), order 5): lam = 2 (principal), 2cos(2pi/5) / 2cos(4pi/5) on the class of
        2x^2+xy+6y^2, 2cos(4pi/5) / 2cos(2pi/5) on the class of 3x^2+xy+4y^2, 0 inert
  D11a, D37a: elliptic curves, lam = a_p / sqrt p (PARI ellap)
  E11a: Sym^2 of 11a (degree 3, Q = 121): lam = (a_p^2/p) - 1 (p != 11), 1/11 at p = 11
  DD: Ramanujan Delta, lam = tau(p)/p^(11/2)
"""
import numpy as np, math

TWO_PI = 2 * math.pi

# ---------------------------------------------------------------- primes
def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]

_SP = [None]
def sieve_odd(L, U):
    """L < U even, L >= 0.  ok[i] is True iff L + 1 + 2i is prime (i < (U - L)/2)."""
    assert L % 2 == 0 and U % 2 == 0 and 0 <= L < U
    r = math.isqrt(U) + 2
    if _SP[0] is None or _SP[0][-1] < r: _SP[0] = primes_upto(max(r, 50000))
    ok = np.ones((U - L) // 2, bool)
    if L == 0: ok[0] = False                                  # 1 is not prime
    for q in _SP[0][1:]:
        q = int(q)
        if q * q >= U: break
        st = max(q * q, -(-(L + 1) // q) * q)
        if st % 2 == 0: st += q
        if st < U: ok[(st - L - 1) // 2::q] = False
    return ok

def isqrt_vec(x):
    """floor(sqrt(x)) for an int64 array x >= 0 (exact)."""
    x = np.asarray(x, np.int64); r = np.floor(np.sqrt(np.maximum(x, 0).astype(np.float64))).astype(np.int64)
    r[(r + 1) * (r + 1) <= x] += 1; r[r * r > x] -= 1
    return r

def _expand(lo, hi, step, other):
    """for each k: values lo[k], lo[k]+step, ..., <= hi[k]; returns (vals, others repeated)."""
    cnt = np.maximum(0, (hi - lo) // step + 1); tot = int(cnt.sum())
    if tot == 0: return np.zeros(0, np.int64), np.zeros(0, np.int64)
    start = np.repeat(np.cumsum(cnt) - cnt, cnt)
    vals = np.repeat(lo, cnt) + step * (np.arange(tot, dtype=np.int64) - start)
    return vals, np.repeat(other, cnt)

# ---------------------------------------------------------------- Q(i): a > b > 0, a + b odd, L < a^2 + b^2 < U
def gaussian_octant(L, U, ok):
    """primes p = a^2 + b^2 in (L, U) with a > b >= 1 (one per prime p = 1 mod 4); returns p, a, b."""
    bmax = math.isqrt(U // 2) + 1
    b = np.arange(1, bmax + 1, dtype=np.int64)
    lo = np.maximum(b + 1, isqrt_vec(L - b * b) - 1); hi = isqrt_vec(U - 1 - b * b)
    lo = lo + ((lo + b) % 2 == 0)                                 # a + b odd
    a, bb = _expand(lo, hi, 2, b)
    n = a * a + bb * bb; m = (n > L) & (n < U) & (a > bb); a, bb, n = a[m], bb[m], n[m]
    m = ok[(n - L - 1) // 2]
    return n[m], a[m], bb[m]

# ---------------------------------------------------------------- Q(sqrt-3): a, b >= 1, a + b odd, L < a^2 + 3b^2 < U
def eisenstein_rep(L, U, ok):
    bmax = math.isqrt(U // 3) + 1
    b = np.arange(1, bmax + 1, dtype=np.int64)
    lo = np.maximum(1, isqrt_vec(L - 3 * b * b) - 1); hi = isqrt_vec(U - 1 - 3 * b * b)
    lo = lo + ((lo + b) % 2 == 0)
    a, bb = _expand(lo, hi, 2, b)
    n = a * a + 3 * bb * bb; m = (n > L) & (n < U); a, bb, n = a[m], bb[m], n[m]
    m = ok[(n - L - 1) // 2]
    return n[m], a[m], bb[m]

# ---------------------------------------------------------------- binary forms A x^2 + B xy + C y^2, y >= 1, given parities
def mark_form(L, U, A, B, C, xpar, ypar):
    """boolean array over odd n in (L, U): n = A x^2 + B x y + C y^2 for some x (parity xpar) and y >= 1 (parity ypar)."""
    D = B * B - 4 * A * C; rep = np.zeros((U - L) // 2, bool)
    ymax = math.isqrt(4 * A * U // (-D)) + 2
    y = np.arange(1 if ypar == 1 else 2, ymax + 1, 2, dtype=np.int64)
    # x = (-B y +- sqrt(4 A n + D y^2)) / (2A), n in (L, U)
    sU = isqrt_vec(4 * A * U + D * y * y) + 2
    sL2 = 4 * A * L + D * y * y; sL = np.where(sL2 > 0, isqrt_vec(np.maximum(sL2, 0)) - 2, -1)
    sL = np.maximum(sL, 0)
    for sgn in (1, -1):
        if sgn == 1: lo = (-B * y + sL) // (2 * A) - 1; hi = (-B * y + sU) // (2 * A) + 1
        else: lo = (-B * y - sU) // (2 * A) - 1; hi = (-B * y - sL) // (2 * A) + 1
        lo = lo + ((lo - xpar) % 2 != 0)
        x, yy = _expand(lo, hi, 2, y)
        n = A * x * x + B * x * yy + C * yy * yy; m = (n > L) & (n < U) & (n % 2 == 1)
        rep[(n[m] - L - 1) // 2] = True
    return rep

def legendre_table(q):
    t = np.full(q, -1, np.int8); t[(np.arange(1, q) ** 2) % q] = 1; t[0] = 0
    return t

# ---------------------------------------------------------------- the cases
C5a, C5b = 2 * math.cos(2 * math.pi / 5), 2 * math.cos(4 * math.pi / 5)
SERIES = {
    # key: (group, d, Q, mu list (analytic), weight k)
    'A4': ('A', 2, 4, [2.0, 3.0], 5), 'A8': ('A', 2, 4, [4.0, 5.0], 9), 'A12': ('A', 2, 4, [6.0, 7.0], 13),
    'A16': ('A', 2, 4, [8.0, 9.0], 17), 'A20': ('A', 2, 4, [10.0, 11.0], 21),
    'B6': ('B', 2, 3, [3.0, 4.0], 7), 'B12': ('B', 2, 3, [6.0, 7.0], 13),
    'C23': ('C23', 2, 23, [0.0, 1.0], 1), 'C31': ('C31', 2, 31, [0.0, 1.0], 1),
    'C47a': ('C47', 2, 47, [0.0, 1.0], 1), 'C47b': ('C47', 2, 47, [0.0, 1.0], 1),
    'D11a': ('D11a', 2, 11, [0.5, 1.5], 2), 'D37a': ('D37a', 2, 37, [0.5, 1.5], 2),
    'E11a': ('D11a', 3, 121, [1.0, 1.0, 2.0], 3),
    'DD': ('DD', 2, 1, [5.5, 6.5], 12),
}
GROUP_DQ = {'A': (2, 4), 'B': (2, 3), 'C23': (2, 23), 'C31': (2, 31), 'C47': (2, 47), 'D11a': None, 'D37a': (2, 37), 'DD': (2, 1)}
# nebentypus at p (for lambda(p^k) recursion), and the special small primes (p = 2, ramified)
def special_primes(key):
    """dict p -> lambda(p) for primes not produced by the odd-norm enumerations (p = 2 and the ramified prime)."""
    if key[0] == 'A':
        m = int(key[1:]); return {2: math.cos(m * math.pi / 4) * 1.0}    # ((1+i)/sqrt2)^m = e^{i m pi/4} = (-1)^(m/4) (single Satake parameter)
    if key[0] == 'B':
        m = int(key[1:]); return {2: 0.0, 3: math.cos(m * math.pi / 2)}  # 2 inert; (sqrt-3/|.|)^m = e^{i m pi/2} = (-1)^(m/6)
    if key in ('C23', 'C31'): return {2: -1.0, int(key[1:]): 1.0}
    if key == 'C47a': return {2: C5a, 47: 1.0}
    if key == 'C47b': return {2: C5b, 47: 1.0}
    return {}

def chi_neben(key, p):
    """nebentypus chi(p) (for degree 2): chi_{-4}, chi_{-3}, chi_D for the weight-1 forms, 1 for D11a, D37a, DD (0 at bad p)."""
    p = np.asarray(p)
    if key[0] == 'A': return np.where(p % 2 == 0, 0, np.where(p % 4 == 1, 1, -1))
    if key[0] == 'B': return np.where(p % 3 == 0, 0, np.where(p % 3 == 1, 1, -1))
    if key[0] == 'C':
        q = int(key[1:3]); t = legendre_table(q).astype(int)
        return np.where(p == 2, 1, t[p % q])                  # (D/2) = 1 for D = 1 mod 8 (all three)
    if key in ('D11a', 'D37a'): N = int(key[1:3]); return np.where(p % N == 0, 0, 1)
    return np.ones_like(p)

def lam_series(group, L, U, ok, keys):
    """primes in (L, U) (odd) with lambda for each key of the group: returns p (int64), {key: lam (float64)}, extra."""
    out = {}
    if group == 'A':
        p, a, b = gaussian_octant(L, U, ok); th = np.arctan2(b, a)
        for k in keys: out[k] = 2 * np.cos(int(k[1:]) * th)
        return p, out, {'a': a, 'b': b, 'th': th}
    if group == 'B':
        p, a, b = eisenstein_rep(L, U, ok); th = np.arctan2(math.sqrt(3) * b, a)
        for k in keys: out[k] = 2 * np.cos(int(k[1:]) * th)
        return p, out, {}
    if group in ('C23', 'C31', 'C47'):
        q = int(group[1:]); idx = np.nonzero(ok)[0]; p = L + 1 + 2 * idx
        leg = legendre_table(q)[p % q]
        cc = {23: 6, 31: 8, 47: 12}[q]
        prin = mark_form(L, U, 1, 1, cc, 1, 0)[idx]
        if q == 47:
            f1 = mark_form(L, U, 2, 1, 6, 1, 1)[idx]
            out['C47a'] = np.where(leg == 1, np.where(prin, 2.0, np.where(f1, C5a, C5b)), 0.0)
            out['C47b'] = np.where(leg == 1, np.where(prin, 2.0, np.where(f1, C5b, C5a)), 0.0)
        else:
            out[group] = np.where(leg == 1, np.where(prin, 2.0, -1.0), 0.0)
        for k in out: out[k][p == q] = 1.0
        return p, out, {}
    raise ValueError(group)

# ---------------------------------------------------------------- the law
def law_constants(d, Q, mu, eps, c):
    """K with sum_{n<=x} Lambda(n) e(-d (n/Q)^(1/d)) ~ K x^((d+1)/(2d)); kappa of the first-order factor; phi0."""
    phi0 = (math.pi / 2) * sum(mu) - d * math.pi / 4
    K = complex(eps) * c * complex(math.cos(phi0 + math.pi / 4), -math.sin(phi0 + math.pi / 4)) * (2 * math.sqrt(d) / (d + 1)) * Q ** (-1 / (2 * d))
    kappa = (d * d - 1) / (24 * d) + 0.5 * sum(m - m * m for m in mu)
    return K, kappa, phi0

def main_term(K, kappa, d, Q, x0, x1, weight=None, npts=20001):
    """integral of weight(x) dM(x) over (x0, x1], dM = ((d+1)/2d) K x^((1-d)/2d) (1 - i kappa/t(x)) dx, t = 2 pi (x/Q)^(1/d).
    weight: function of v = x^(1/d) (vectorised) or None (sharp, weight 1). Integration in v (dx = d v^(d-1) dv)."""
    v = np.linspace(x0 ** (1 / d), x1 ** (1 / d), npts); x = v ** d
    t = TWO_PI * (x / Q) ** (1 / d)
    dens = ((d + 1) / (2 * d)) * K * x ** ((1 - d) / (2 * d)) * (1 - 1j * kappa / t) * d * v ** (d - 1)
    if weight is not None: dens = dens * weight(v)
    # Simpson
    h = v[1] - v[0]; s = dens[0] + dens[-1] + 4 * dens[1:-1:2].sum() + 2 * dens[2:-1:2].sum()
    return s * h / 3

# ---------------------------------------------------------------- windows (in v = p^(1/d)) and statistics
def hann(v, v0, v1):
    u = (v - v0) / (v1 - v0); return np.where((u > 0) & (u < 1), np.sin(np.pi * u) ** 2, 0.0)

def windows(v0, v1, nsub=8):
    """main Hann window on [v0, v1] and nsub half-overlapping Hann sub-windows."""
    h = (v1 - v0) / (nsub + 1)
    W = [(v0, v1)] + [(v0 + j * h, v0 + (j + 2) * h) for j in range(nsub)]
    return W

CTRL = [0.6713, 0.7919, 0.8837, 0.9431, 1.0617, 1.1237, 1.2389, 1.3711]   # control frequency factors (alpha^2 far from a/b, a,b <= 12)

def seg_stats(p, lams, d, Q, W, ctrl=CTRL):
    """per-segment sums: sharp (resonant + controls), smooth for each window (resonant + controls), noise sums."""
    pf = p.astype(np.float64); lp = np.log(pf); v = pf ** (1.0 / d)
    ph = TWO_PI * d * (pf / Q) ** (1.0 / d)
    E = [np.exp(-1j * ph)] + [np.exp(-1j * a * ph) for a in ctrl]
    g = pf ** (-0.5 + 0.5 / d)
    Ws = [hann(v, a, b) for (a, b) in W]
    res = {}
    for k, lam in lams.items():
        wl = lam * lp
        r = {'sharp': np.array([np.dot(wl, e) for e in E]),
             'noise_sharp': float(np.dot(wl, wl))}
        sm = np.zeros((len(W), len(E)), complex); nz = np.zeros(len(W))
        for i, w in enumerate(Ws):
            ww = wl * g * w
            nz[i] = float(np.dot(ww, ww))
            for j, e in enumerate(E): sm[i, j] = np.dot(ww, e)
        r['smooth'] = sm; r['noise_smooth'] = nz; r['count'] = int((lam != 0).sum())
        res[k] = r
    return res

# ---------------------------------------------------------------- Pythagorean rays (Q(i))
def pyth_directions(H):
    """primitive Pythagorean directions w = (u, v), u > 0, v > 0, |w| = h <= H; weight A_h/(4/pi) = prod_{p|h} (-1/(p-1)) if h squarefree else 0."""
    out = []
    for s in range(2, int(math.isqrt(H)) + 2):
        for t in range(1, s):
            if (s - t) % 2 == 1 and math.gcd(s, t) == 1:
                h = s * s + t * t
                if h > H: continue
                u, v = s * s - t * t, 2 * s * t
                # weight
                hh, wgt, sqf = h, 1.0, True; q = 2
                while q * q <= hh:
                    if hh % q == 0:
                        e = 0
                        while hh % q == 0: hh //= q; e += 1
                        if e > 1: sqf = False
                        wgt *= -1.0 / (q - 1)
                    q += 1
                if hh > 1: wgt *= -1.0 / (hh - 1)
                if not sqf: wgt = 0.0
                for (x, y) in ((u, v), (v, u)): out.append((math.atan2(y, x), h, wgt, x, y))
    return sorted(out)

BETA = 2.0
def zone_F(delta, beta=BETA, half=False):
    """F(delta) = int w(u + delta) e^{-i pi u^2} du, w(u) = exp(-pi u^2/beta^2) (u over R, or u > 0 if half)."""
    if not half:
        a = 1 / beta ** 2 + 1j
        return a ** -0.5 * np.exp(-1j * np.pi * delta ** 2 / (1 + 1j * beta ** 2))
    u = np.linspace(0, 12 * beta, 24001); f = np.exp(-np.pi * (u + delta) ** 2 / beta ** 2 - 1j * np.pi * u * u)
    return np.trapezoid(f, u) if hasattr(np, 'trapezoid') else np.trapz(f, u)

def zone_sums(pf, a, b, phis, ms=(0, 4), beta=BETA, weight=None):
    """Gaussian points z = a + bi with a > b > 0 (octant) and their mirrors b + ai (together the open first quadrant):
    Z_m(phi) = sum_{z in first quadrant} wt(z) e^{i m arg z} e^{-2 pi i |z|} w(s/sqrt t), t = <z, e_phi>, s = <z, i e_phi>,
    w(u) = exp(-pi u^2/beta^2), wt = ln N(z) (primes) or the given weight.  Mirror points have arg = pi/2 - theta,
    e^{i m (pi/2 - theta)} = e^{-i m theta} for m = 0 mod 4."""
    a = a.astype(np.float64); b = b.astype(np.float64)
    r = np.sqrt(pf); lp = np.log(pf) if weight is None else weight; th = np.arctan2(b, a); base = lp * np.exp(-1j * TWO_PI * r)
    out = {}
    for ph in phis:
        cs, sn = math.cos(ph), math.sin(ph)
        acc = {mm: 0j for mm in ms}; nz = 0.0
        for (x, y, sg) in ((a, b, 1), (b, a, -1)):
            t = x * cs + y * sn; s = -x * sn + y * cs
            u = s / np.sqrt(np.maximum(t, 1e-9)); m = (t > 0) & (np.abs(u) < 6 * beta)
            if not m.any(): continue
            ww = np.exp(-np.pi * u[m] ** 2 / beta ** 2); w = ww * base[m]
            nz += float(np.sum((ww * lp[m]) ** 2)) if np.ndim(lp) else float(np.sum(ww * ww)) * lp * lp
            for mm in ms: acc[mm] += complex(np.sum(w * np.exp(1j * sg * mm * th[m])))
        for mm in ms: out[(ph, mm)] = acc[mm]
        out[(ph, 'noise')] = nz
    return out

def zone_prediction(phi_c, ms, r0, r1, dirs, beta=BETA, nr=401):
    """ray model: Z_m(phi_c) = sum_rays A e^{i m phi_h} int_{r0}^{r1} sqrt r F(sqrt r sin(phi_h - phi_c)) dr;
    rays: the two half axes (A = 4/pi, half zones) and the Pythagorean directions in dirs (A = (4/pi) * weight)."""
    rr = np.linspace(r0, r1, nr); out = {mm: 0j for mm in ms}
    rays = [(0.0, 4 / math.pi, 'half+'), (math.pi / 2, 4 / math.pi, 'half-')] + [(ph, 4 / math.pi * wgt, 'full') for (ph, h, wgt, x, y) in dirs if wgt != 0]
    for (ph, A, kind) in rays:
        dl = np.sqrt(rr) * math.sin(ph - phi_c)
        if np.min(np.abs(dl)) > 6 * beta: continue
        if kind == 'full': F = zone_F(dl, beta)
        else:
            F = np.array([zone_F(x if kind == 'half+' else -x, beta, half=True) for x in dl])
        f = np.sqrt(rr) * F; J = np.trapezoid(f, rr) if hasattr(np, 'trapezoid') else np.trapz(f, rr)
        for mm in ms: out[mm] += A * complex(math.cos(mm * ph), math.sin(mm * ph)) * J
    return out
