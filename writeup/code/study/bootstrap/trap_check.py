"""Independent numerical check of the traps for the double zeros (section 5g of bootstrap.md).

Nothing of the Bessel code is used.  Every function h with transform Xi(t)^2 B_h(t) is evaluated from the zero side,
    h(u) = (1/pi) int_0^oo Xi(t)^2 B_h(t) cos(tu) dt,   h'(u) = -(1/pi) int_0^oo Xi^2 B_h t sin(tu) dt,  h'' likewise,
by Gauss-Legendre quadrature (120 nodes on each [k, k+1], k < 400, 520 bits), as in cert_check.py.

Input: the optimum of half-degree d' (file theory/convex/certs/polish_A_d<d'>.json: coefficients b with g(0) = 1,
double zeros u_k, weights mu_k, cost lam).

For each node k:
  h_k in V_d' with h_k(u_j) = h_k'(u_j) = 0 (j != k), h_k(u_k) = 1, h_k'(u_k) = 0, top coefficient 0;
  W(h_k) = 2 sum Lambda(n) n^(-1/2) h_k(log n)  and the test  W(h_k) = lam h_k(0) + mu_k  (the KKT identity);
  eta_0 = lam / W(h_k);  eta = FACTOR * eta_0;  condition (i): eta B_h <= B_g on the real line (grid and roots);
  the set S = {u in (0, 3.8] : eta h_k(u) > g(u)} on a grid (step 1/200) and its endpoints by bisection;
  trap radius = largest distance from log n_k to a point of S.
Beyond u = 3.8 the values are below the working precision; there w = eta h - g has the sign of -b_top(g) < 0.

Run: (ulimit -v 1500000; timeout 2400 nice python trap_check.py [half-degree d' = 24] [FACTOR = 2] [cache.pkl])
"""
import sys, os, json, ast, pickle
import numpy as np
import gmpy2
from gmpy2 import mpfr
import mpmath as mp

PREC = 520
mp.mp.prec = PREC
gmpy2.get_context().precision = PREC
T, NODES = 400, 120
HERE = os.path.dirname(os.path.abspath(__file__))
CERT = os.path.join(HERE, 'theory', 'convex', 'certs')
DH = int(sys.argv[1]) if len(sys.argv) > 1 else 24
FACTOR = sys.argv[2] if len(sys.argv) > 2 else '2'
CACHE = sys.argv[3] if len(sys.argv) > 3 else None


def gl_nodes(n):
    xs, ws = [], []
    for i in range(1, n + 1):
        x = mp.cos(mp.pi * (i - mp.mpf(1) / 4) / (n + mp.mpf(1) / 2))
        for _ in range(60):
            p0, p1 = mp.mpf(1), x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < mp.mpf(2) ** (-PREC + 20):
                break
        p0, p1 = mp.mpf(1), x
        for k in range(2, n + 1):
            p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
        dp = n * (x * p1 - p0) / (x * x - 1)
        xs.append(x)
        ws.append(2 / ((1 - x * x) * dp * dp))
    return xs, ws


def Xi(t):
    s = mp.mpf(1) / 2 + 1j * t
    return (s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)).real


def G(x):
    return mpfr(mp.nstr(x, 160, strip_zeros=False))


def build():
    if CACHE and os.path.exists(CACHE):
        with open(CACHE, 'rb') as f:
            t, w = pickle.load(f)
        return np.array([mpfr(x) for x in t], dtype=object), np.array([mpfr(x) for x in w], dtype=object)
    xs, ws = gl_nodes(NODES)
    t, w = [], []
    for k in range(T):
        for x, wt in zip(xs, ws):
            tt = k + (x + 1) / 2
            t.append(G(tt))
            w.append(G(wt / 2 * Xi(tt) ** 2 / mp.pi))
    if CACHE:
        with open(CACHE, 'wb') as f:
            pickle.dump(([str(x) for x in t], [str(x) for x in w]), f)
    return np.array(t, dtype=object), np.array(w, dtype=object)


t, w = build()
TAU = mpfr(30)                         # basis (t/TAU)^(2j), to keep the linear system tame
s2 = t * t / (TAU * TAU)
pw = [np.array([mpfr(1)] * len(t), dtype=object)]
for j in range(DH):
    pw.append(pw[-1] * s2)
wpw = [w * p for p in pw]              # weights times basis
ZERO = mpfr(0)


def tot(a):
    return gmpy2.fsum(list(a))


def basis_at(u, der=0):
    """values of the der-th derivative of the DH+1 basis functions at u"""
    u = mpfr(u)
    if der == 0:
        c = np.array([gmpy2.cos(x * u) for x in t], dtype=object)
    elif der == 1:
        c = -np.array([gmpy2.sin(x * u) for x in t], dtype=object) * t
    else:
        c = -np.array([gmpy2.cos(x * u) for x in t], dtype=object) * t * t
    return [tot(wp * c) for wp in wpw]


def val(coef, bas):
    return sum(c * b for c, b in zip(coef, bas))


def mangoldt(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0:
                m //= p
            return gmpy2.log(mpfr(p)) if m == 1 else ZERO
    return ZERO


def e(x, d=6):
    return mp.nstr(mp.mpf(str(x)), d)


c = json.load(open(os.path.join(CERT, 'polish_A_d%d.json' % DH)))
braw = c['b'] if not isinstance(c['b'], str) else ast.literal_eval(c['b'])
bg = [mpfr(x) * TAU ** (2 * j) for j, x in enumerate(braw)]          # coefficients of g in the scaled basis
U = [mpfr(x) for x in c['touch_points']]
MU = [mpfr(x) for x in c['mu']]
LAM = mpfr(c['lam'])
K = len(U)
PPOW = [n for n in range(2, 200) if mangoldt(n) != 0]
print('half-degree d\' = %d (degree %d), K = %d nodes, lambda = %s, FACTOR = %s' % (DH, 2 * DH, K, e(LAM, 10), FACTOR))

# prime-side functional on the basis, and the basis at 0, at the nodes
NMAX = 47
Wb = [ZERO] * (DH + 1)
for n in PPOW:
    if n > NMAX:
        break
    bas = basis_at(gmpy2.log(mpfr(n)))
    f = 2 * mangoldt(n) / gmpy2.sqrt(mpfr(n))
    Wb = [a + f * b for a, b in zip(Wb, bas)]
B0 = basis_at(0)
BU = [basis_at(u) for u in U]
DU = [basis_at(u, 1) for u in U]
D2U = [basis_at(u, 2) for u in U]
print('checks of the stored optimum: g(0) = %s, W(g) = %s, largest |g(u_k)| = %s, largest |g\'(u_k)| = %s'
      % (e(val(bg, B0), 12), e(val(bg, Wb), 10), e(max(abs(val(bg, b)) for b in BU), 3), e(max(abs(val(bg, b)) for b in DU), 3)))

# grid for the set S
GRID = [mpfr(i) / 200 for i in range(1, 761)]
BG = [basis_at(u) for u in GRID]
gG = [val(bg, b) for b in BG]

print('  k |  n | u_k - log n | mu_k       | 2L/sqrt n  | g\'\'(u_k)   | W(h_k)     | KKT defect | h_k(0)     | eta        | eta_1/eta  | S (grid) | a - log n   | b - log n   | radius     | sqrt(2 eta/g\'\') | weight >=')
for k in range(K):
    n = PPOW[k]
    rows, rhs = [], []
    for j in range(K):
        rows.append(BU[j]); rhs.append(mpfr(1) if j == k else ZERO)
        rows.append(DU[j]); rhs.append(ZERO)
    top = [ZERO] * DH + [mpfr(1)]
    rows.append(top); rhs.append(ZERO)
    A = mp.matrix([[mp.mpf(str(x)) for x in r] for r in rows])
    rr = [mp.mpf(str(x)) for x in rhs]
    for i in range(A.rows):
        m = max(abs(A[i, j]) for j in range(A.cols))
        for j in range(A.cols):
            A[i, j] /= m
        rr[i] /= m
    mp.mp.prec = 3 * PREC
    sol = mp.lu_solve(A, mp.matrix(rr))
    mp.mp.prec = PREC
    ch = [G(x) for x in sol]
    Wh = val(ch, Wb)
    h0 = val(ch, B0)
    defect = (Wh - LAM * h0 - MU[k]) / MU[k]
    eta0 = LAM / Wh
    eta = mpfr(FACTOR) * eta0
    g2 = val(bg, D2U[k])
    # condition (i): eta B_h <= B_g on the real line; ratio B_g/B_h where B_h > 0, on a grid in s = t^2/TAU^2
    mp.mp.prec = 2000
    try:
        diff = [mp.mpf(str(a)) - mp.mpf(str(eta)) * mp.mpf(str(b)) for a, b in zip(bg, ch)]
        while abs(diff[-1]) == 0:
            diff.pop()
        r = mp.polyroots(list(reversed(diff)), maxsteps=3000, extraprec=6000)
        nreal = sum(1 for z in r if abs(z.imag) < abs(z) * mp.mpf(10) ** (-30) and z.real >= 0)
    finally:
        mp.mp.prec = PREC
    eta1 = None
    for i in range(0, 4001):
        s = mpfr(i) / 100
        bh = sum(cc * s ** j for j, cc in enumerate(ch))
        if bh > 0:
            q = sum(cc * s ** j for j, cc in enumerate(bg)) / bh
            eta1 = q if eta1 is None or q < eta1 else eta1
    # the set S on the grid
    wG = [eta * val(ch, b) - g for b, g in zip(BG, gG)]
    pos = [i for i, v in enumerate(wG) if v > 0 and GRID[i] <= mpfr('3.8')]
    comp = 0
    for idx, i in enumerate(pos):
        if idx == 0 or pos[idx - 1] != i - 1:
            comp += 1
    # endpoints around u_k by bisection on w(u) = eta h(u) - g(u), starting from u_k (where w = eta > 0)
    def wv(u):
        bas = basis_at(u)
        return eta * val(ch, bas) - val(bg, bas)
    rad0 = gmpy2.sqrt(2 * eta / g2)
    ends = []
    for sgn in (-1, 1):
        lo, hi = U[k], U[k] + sgn * 4 * rad0
        if wv(hi) > 0:
            ends.append(None)
            continue
        for _ in range(60):
            mid = (lo + hi) / 2
            if wv(mid) > 0:
                lo = mid
            else:
                hi = mid
        ends.append(hi)
    ln = gmpy2.log(mpfr(n))
    if None in ends:
        a_s = b_s = rad_s = 'n/a'
    else:
        a_s, b_s = e(ends[0] - ln, 4), e(ends[1] - ln, 4)
        rad_s = e(max(abs(ends[0] - ln), abs(ends[1] - ln)), 4)
    # grid points in S that are not inside (a, b)
    far = [float(GRID[i]) for i in pos if None in ends or not (ends[0] <= GRID[i] <= ends[1])]
    wlow = (eta * Wh - LAM * (1 + max(eta * h0 - 1, 0))) / eta           # sup_S w <= eta (h_k <= 1 near u_k up to o(1))
    print('%3d | %2d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %d pts, %d comp, outside (a,b): %d | %s | %s | %s | %s | %s (roots of B_g - eta B_h with s >= 0: %d)'
          % (k + 1, n, e(U[k] - ln, 4), e(MU[k], 8), e(2 * mangoldt(n) / gmpy2.sqrt(mpfr(n)), 8), e(g2, 4), e(Wh, 8), e(defect, 2),
             e(h0, 4), e(eta, 4), e(eta1 / eta, 3) if eta1 is not None else 'inf', len(pos), comp, len(far), a_s, b_s, rad_s, e(rad0, 4), e(wlow, 6), nreal), flush=True)
