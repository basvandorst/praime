"""consec_lib.py: primes in a segment and their Frobenius-type labels (consecutive-primes study).

Self-contained (the sieve and lattice enumeration follow writeup/code/study/resonance/primeside/ps_lib.py, copied here).
Every label is read off from the representation of p by binary quadratic forms, from the angle of p = a^2 + b^2, or
from a_p of the elliptic curve 37a1 (PARI ellap).  Codes (int8), -1 = excluded (p = 2, ramified p):

  C23, C31  (S3: x^3 - x - 1, x^3 + x - 1):  0 inert ((p/q) = -1, one root), 1 principal (3 roots), 2 other split (0 roots)
  C47       (D5, class group Z/5):           0 inert, 1 principal, 2 class of (2, +-1, 6), 3 class of (3, +-1, 4)
  C39       (D4, class group Z/4):           0 inert, 1 principal, 2 class of (3, 3, 4) (same genus), 3 class of (2, +-1, 5)
  G         (Gaussian, octant angle theta of p = a^2 + b^2, a > b > 0, theta in (0, pi/4)):
            0 for p = 3 mod 4, 1..4 = sector floor(theta / (pi/16)) + 1;  theta kept as float
  E37       a_p of y^2 + y = x^3 - x (int);  derived labels in consec_tables.py

Non-abelian functions (mean 0, mean square 1 over the split primes under Chebotarev / Sato-Tate):
  dihedral: f = a_psi(p)/sqrt2 with a_psi(p) = psi(P) + psi(Pbar) (psi a class-group character of order 3, 5 (two), 4)
  G: f_m = sqrt2 cos(4 m theta), m = 1, 2;   E37: U_1 = a_p/sqrt p, U_2 = U_1^2 - 1, sign(a_p).
"""
import numpy as np, math

# ---------------------------------------------------------------- primes (odd numbers in (L, U), L, U even)
def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]

_SP = [None]
def sieve_odd(L, U):
    """ok[i] True iff L + 1 + 2i is prime, 0 <= i < (U - L)/2."""
    assert L % 2 == 0 and U % 2 == 0 and 0 <= L < U
    r = math.isqrt(U) + 2
    if _SP[0] is None or _SP[0][-1] < r: _SP[0] = primes_upto(max(r, 50000))
    ok = np.ones((U - L) // 2, bool)
    if L == 0: ok[0] = False
    for q in _SP[0][1:]:
        q = int(q)
        if q * q >= U: break
        st = max(q * q, -(-(L + 1) // q) * q)
        if st % 2 == 0: st += q
        if st < U: ok[(st - L - 1) // 2::q] = False
    return ok

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

def mark_form(L, U, A, B, C):
    """boolean over odd n in (L, U): n = A x^2 + B x y + C y^2 for some integers x, y with y >= 1.
    (y = 0 gives A x^2, never an odd prime > A for A in our forms.)  Enumerates y in chunks to bound memory."""
    D = B * B - 4 * A * C; rep = np.zeros((U - L) // 2, bool)
    ymax = math.isqrt(4 * A * U // (-D)) + 2
    for y0 in range(1, ymax + 1, 500):
        y = np.arange(y0, min(ymax, y0 + 499) + 1, dtype=np.int64)
        sU = isqrt_vec(4 * A * U + D * y * y) + 2
        sL2 = 4 * A * L + D * y * y; sL = np.where(sL2 > 0, isqrt_vec(np.maximum(sL2, 0)) - 2, 0); sL = np.maximum(sL, 0)
        for sgn in (1, -1):
            if sgn == 1: lo = (-B * y + sL) // (2 * A) - 1; hi = (-B * y + sU) // (2 * A) + 1
            else: lo = (-B * y - sU) // (2 * A) - 1; hi = (-B * y - sL) // (2 * A) + 1
            x, yy = _expand(lo, hi, 1, y)
            if len(x) == 0: continue
            n = A * x * x + B * x * yy + C * yy * yy; m = (n > L) & (n < U) & (n % 2 == 1)
            rep[(n[m] - L - 1) // 2] = True
    return rep

def legendre_table(q):
    t = np.full(q, -1, np.int8); t[(np.arange(1, q) ** 2) % q] = 1; t[0] = 0
    return t

# ---------------------------------------------------------------- the dihedral cases
FORMS = {  # disc: (principal, [(label, form), ...] for the other classes; the last class = the rest of the split primes)
    23: ((1, 1, 6), []),            # rest (2, +-1, 3) -> label 2
    31: ((1, 1, 8), []),            # rest (2, +-1, 4) -> label 2
    47: ((1, 1, 12), [(2, (2, 1, 6))]),   # rest (3, +-1, 4) -> label 3
    39: ((1, 1, 10), [(2, (3, 3, 4))]),   # rest (2, +-1, 5) -> label 3
}
def dihedral_labels(D, L, U, ok, p):
    idx = (p - L - 1) // 2
    # chi_{-D}(p) = (p/D) for D = 23, 31, 47 (D = 3 mod 4 prime); chi_{-39}(p) = (-3/p)(13/p) = (p/3)(p/13)
    if D == 39: leg = (legendre_table(3)[p % 3] * legendre_table(13)[p % 13]).astype(np.int8)
    else: leg = legendre_table(D)[p % D]
    prin_f, others = FORMS[D]
    lab = np.where(leg == 1, np.int8(2 if D in (23, 31) else 3), np.int8(0)).astype(np.int8)
    for code, f in others:
        m = mark_form(L, U, *f)[idx]; lab[(leg == 1) & m] = code
    m = mark_form(L, U, *prin_f)[idx]; lab[(leg == 1) & m] = 1
    lab[(p % D == 0) | (p == 2)] = -1
    if D == 39: lab[(p == 3) | (p == 13)] = -1
    return lab

# a_psi / sqrt2 per label code (index = label code; label -1 maps to nan)
C5a, C5b = 2 * math.cos(2 * math.pi / 5), 2 * math.cos(4 * math.pi / 5)
F_DIH = {
    23: [np.array([0.0, 2.0, -1.0]) / math.sqrt(2)],
    31: [np.array([0.0, 2.0, -1.0]) / math.sqrt(2)],
    47: [np.array([0.0, 2.0, C5a, C5b]) / math.sqrt(2), np.array([0.0, 2.0, C5b, C5a]) / math.sqrt(2)],
    39: [np.array([0.0, 2.0, -2.0, 0.0]) / math.sqrt(2)],
}
CHEB = {  # Chebotarev probability of each label code
    23: [1 / 2, 1 / 6, 1 / 3], 31: [1 / 2, 1 / 6, 1 / 3], 47: [1 / 2, 1 / 10, 1 / 5, 1 / 5], 39: [1 / 2, 1 / 8, 1 / 8, 1 / 4],
}

# ---------------------------------------------------------------- Gaussian primes
def gaussian_theta(L, U, ok, p):
    """theta in (0, pi/4) for p = 1 mod 4 (nan otherwise)."""
    bmax = math.isqrt(U // 2) + 1; th = np.full(len(p), np.nan)
    for b0 in range(1, bmax + 1, 1000):
        b = np.arange(b0, min(bmax, b0 + 999) + 1, dtype=np.int64)
        lo = np.maximum(b + 1, isqrt_vec(L - b * b) - 1); hi = isqrt_vec(U - 1 - b * b)
        lo = lo + ((lo + b) % 2 == 0)
        a, bb = _expand(lo, hi, 2, b)
        if len(a) == 0: continue
        n = a * a + bb * bb; m = (n > L) & (n < U) & (a > bb); a, bb, n = a[m], bb[m], n[m]
        m = ok[(n - L - 1) // 2]; a, bb, n = a[m], bb[m], n[m]
        pos = np.searchsorted(p, n); th[pos] = np.arctan2(bb, a)
    return th

def gaussian_labels(th, p):
    lab = np.where(p % 4 == 3, 0, 1 + np.floor(np.nan_to_num(th) / (math.pi / 16))).astype(np.int8)
    lab[(p % 4 == 1) & ~np.isfinite(th)] = -1          # must not happen
    lab[p == 2] = -1
    return lab

# ---------------------------------------------------------------- elliptic curve 37a1
_pari = [None]
def pari():
    if _pari[0] is None:
        import cypari2
        _pari[0] = cypari2.Pari(); _pari[0].allocatemem(300 * 10 ** 6)
    return _pari[0]
def ec_ap(L, U, p):
    """a_p of 37a1 for the primes p (sorted, in (L, U))."""
    P = pari(); out = np.zeros(len(p), np.int64)
    step = 2 * 10 ** 6
    for lo in range(L, U, step):
        hi = min(U, lo + step)
        v = P('my(E=ellinit([0,0,1,-1,0]), Pr=primes([%d,%d])); [Vecsmall(Pr), Vecsmall(apply(q->ellap(E,q), Pr))]' % (lo + 1, hi - 1))
        pr = np.array(v[0].list(), dtype=np.int64) if len(v[0]) else np.zeros(0, np.int64)
        ap = np.array(v[1].list(), dtype=np.int64) if len(v[1]) else np.zeros(0, np.int64)
        m = pr > 2; pr, ap = pr[m], ap[m]
        pos = np.searchsorted(p, pr); assert np.all(p[pos] == pr)
        out[pos] = ap
    return out

# ---------------------------------------------------------------- one segment
def segment(L, U, cases=('C23', 'C31', 'C47', 'C39', 'G', 'E37')):
    """odd primes in (L, U) and their labels.  Returns dict with 'p' (int64) and per case arrays."""
    ok = sieve_odd(L, U); p = L + 1 + 2 * np.nonzero(ok)[0].astype(np.int64)
    out = {'p': p}
    for c in cases:
        if c[0] == 'C': out[c] = dihedral_labels(int(c[1:]), L, U, ok, p)
        elif c == 'G':
            th = gaussian_theta(L, U, ok, p); out['G_th'] = th; out['G'] = gaussian_labels(th, p)
        elif c == 'E37': out['E37_ap'] = ec_ap(L, U, p)
    return out
