"""Feature factory: groups of features from outside elementary arithmetic.

Each group: (group_name, domain, mode, fn). fn(ctx, idx) -> list of (name, kind, values)
  idx  : int64 indices into ctx.P (the sorted primes); k = ctx.P[idx]
  kind : 'bin' (0/1) or 'cont' (threshold chosen on discovery deciles); NaN = undefined
  mode : 'full' (every validation prime) or 'sample' (importance sample, HT-weighted)
"""
import math, hashlib
import numpy as np
import mpmath, gmpy2

LN2 = math.log(2.0)
U64 = np.uint64


class Ctx:
    def __init__(self, P, NMAX, zeros):
        self.P = P.astype(np.int64)
        self.NMAX = int(NMAX)          # largest k for which sequences must exist
        self.zeros = zeros             # float64 array of gamma
        self._cache = {}
        self.lnP = np.log(self.P.astype(np.float64))

    def cache(self, key, maker):
        if key not in self._cache:
            if len(self._cache) > 3:
                self._cache.pop(next(iter(self._cache)))
            self._cache[key] = maker()
        return self._cache[key]


# ---------------- exact fractional parts frac(k*alpha) ----------------
def alpha_words(alpha_mp):
    """frac(alpha) as two 64-bit words A (top) and A2 (next)."""
    with mpmath.workdps(60):
        f = alpha_mp - mpmath.floor(alpha_mp)
        x = int(mpmath.floor(f * mpmath.mpf(2) ** 128))
    return x >> 64, (x & ((1 << 64) - 1))


def frac_mul(k, alpha_mp):
    """frac(k*alpha) for int64 k < 2^53, error < 1e-13 (k^2 range) / 1e-18 (k < 2^31)."""
    A, A2 = alpha_words(alpha_mp)
    x = np.asarray(k, dtype=np.uint64)
    A2hi = U64(A2 >> 32)
    x1 = x >> U64(32); x0 = x & U64(0xFFFFFFFF)
    with np.errstate(over="ignore"):
        lo = x * U64(A)                                   # wraps mod 2^64
        hi = x1 * A2hi + ((x0 * A2hi) >> U64(32))         # floor(x*A2hi/2^32)
        u = lo + hi
    return (u >> U64(11)).astype(np.float64) / float(1 << 53)


def frac(x):
    return x - np.floor(x)


def popcount(x):
    return np.bitwise_count(np.asarray(x, dtype=np.uint64)).astype(np.int64)


def _make_consts():
    with mpmath.workdps(60):
        CONSTS = {
            "pi": mpmath.pi, "e": mpmath.e, "phi": (1 + mpmath.sqrt(5)) / 2, "sqrt2": mpmath.sqrt(2),
            "sqrt3": mpmath.sqrt(3), "sqrt5": mpmath.sqrt(5), "ln2": mpmath.log(2), "ln3": mpmath.log(3),
            "ln10": mpmath.log(10), "zeta3": mpmath.zeta(3), "catalan": mpmath.catalan,
            "eulergamma": mpmath.euler, "feigdelta": mpmath.mpf("4.669201609102990671853203820466"),
            "feigalpha": mpmath.mpf("2.502907875095892822283902873218"), "khinchin": mpmath.khinchin,
            "glaisher": mpmath.glaisher, "mertens": mpmath.mertens, "twinprime": mpmath.twinprime,
            "inv2pi": 1 / (2 * mpmath.pi), "inv_year": 1 / mpmath.mpf("365.24219"),
            "inv_synmonth": 1 / mpmath.mpf("29.530589"), "inv_alphaEM": 1 / mpmath.mpf("137.035999084"),
            "plastic": mpmath.mpf("1.324717957244746025960908854478"), "silver": 1 + mpmath.sqrt(2),
            "e_pi": mpmath.exp(mpmath.pi), "pi_e": mpmath.pi ** mpmath.e, "semitone": mpmath.mpf(2) ** (mpmath.mpf(1) / 12),
            "golomb_dickman": mpmath.mpf("0.6243299885435508709929363831"), "omega": mpmath.lambertw(1).real,
            "conway": mpmath.mpf("1.303577269034296391257099112153"), "inv_day_s": 1 / mpmath.mpf(86400),
        }
        return {k: +v for k, v in CONSTS.items()}


CONSTS = _make_consts()



def g_constants(ctx, idx):
    k = ctx.P[idx]
    for nm, a in CONSTS.items():
        f = frac_mul(k, a)
        yield ((f"frac(k*{nm})", "cont", f))
        yield ((f"cos2pi(k*{nm})", "cont", np.cos(2 * np.pi * f)))


def g_quadratic(ctx, idx):
    k = ctx.P[idx]
    for nm in ("sqrt2", "pi", "phi", "e"):
        f = frac_mul(k * k, CONSTS[nm])
        yield ((f"frac(k^2*{nm})", "cont", f))


LOGB = {"2": 2.0, "e": math.e, "3": 3.0, "5": 5.0, "6": 6.0, "7": 7.0, "10": 10.0, "12": 12.0,
        "pi": math.pi, "phi": (1 + 5 ** 0.5) / 2, "silver": 1 + 2 ** 0.5, "semitone": 2 ** (1 / 12), "1.5": 1.5}


def g_logs(ctx, idx):
    lk = ctx.lnP[idx]
    for nm, b in LOGB.items():
        f = frac(lk / math.log(b))
        yield ((f"frac(log_{nm} k)", "cont", f))
        yield ((f"cos2pi(log_{nm} k)", "cont", np.cos(2 * np.pi * f)))
    # pitch class of a tone of k Hz (A4 = 440 Hz, equal temperament)
    pc = np.mod(np.rint(12 * (lk - math.log(440.0)) / LN2).astype(np.int64), 12)
    for c, nm in enumerate(["A", "A#", "B", "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#"]):
        yield ((f"pitch class {nm} (k Hz)", "bin", (pc == c).astype(float)))
    yield (("cents off pitch (k Hz)", "cont", 12 * (lk - math.log(440.0)) / LN2 - np.rint(12 * (lk - math.log(440.0)) / LN2)))


# ---------------- zeta zeros ----------------
def zero_kernel(gam, N):
    """cosine fade on the top seventh of the first N zeros (project's faded cut); 0 beyond N."""
    K = np.zeros(len(gam))
    K[:N] = 1.0
    c = int(round(N * 6 / 7))
    t = np.arange(c, N)
    K[c:N] = 0.5 * (1 + np.cos(np.pi * (t - c + 1) / (N - c + 1)))
    return K


def g_zeta_sums(ctx, idx):
    k = ctx.P[idx].astype(np.float64)
    gam = ctx.zeros
    Ns = (len(gam), 300, 50)
    Ks = np.stack([zero_kernel(gam, N) for N in Ns], axis=1)        # [nz, 3]
    Kr = zero_kernel(gam, len(gam)) * (1.0 / (0.25 + gam ** 2))        # for psi reconstruction
    xs = {"k": np.log(k), "k-1": np.log(k - 1), "k+1": np.log(k + 1), "2k+1": np.log(2 * k + 1),
          "2k-1": np.log(2 * k - 1), "M_k": k * LN2}
    res = {}
    CH = 1500
    for xn, L in xs.items():
        C = np.zeros((len(k), 3)); S = np.zeros((len(k), 3)); R = np.zeros(len(k))
        for s in range(0, len(k), CH):
            ph = np.multiply.outer(L[s:s + CH], gam)
            c = np.cos(ph); sn = np.sin(ph)
            C[s:s + CH] = c @ Ks; S[s:s + CH] = sn @ Ks
            if xn == "k":
                R[s:s + CH] = (0.5 * c + sn * gam) @ Kr
        res[xn] = (C, S, R)
    sumK = Ks.sum(axis=0)
    for xn, (C, S, R) in res.items():
        sx = np.sqrt(k if xn in ("k", "k-1", "k+1") else (2 * k if xn != "M_k" else np.ones_like(k)))
        for j, N in enumerate(Ns):
            if xn != "M_k":
                yield ((f"zeta det (2/sqrt x)sum cos, x={xn}, N={N}", "cont", 2 * C[:, j] / sx))
            yield ((f"zeta det normalized -sum cos/sumK, x={xn}, N={N}", "cont", -C[:, j] / sumK[j]))
        if xn in ("k", "M_k", "2k+1"):
            for j, N in ((0, Ns[0]), (2, 50)):
                ang = frac(np.arctan2(S[:, j], C[:, j]) / (2 * np.pi))
                yield ((f"zeta phase frac, x={xn}, N={N}", "cont", ang))
                yield ((f"zeta phase cos, x={xn}, N={N}", "cont", np.cos(2 * np.pi * ang)))
        if xn in ("k", "M_k"):
            for j, N in ((0, Ns[0]), (1, 300)):
                yield ((f"zeta |sum x^(i g)|/sumK, x={xn}, N={N}", "cont", np.hypot(C[:, j], S[:, j]) / sumK[j]))
        if xn == "k":
            psi_zero = -2 * np.sqrt(k) * R              # truncated zero sum ~ psi(k) - k
            yield (("zeta sum -2Re(k^rho/rho)/sqrt k (N=all)", "cont", psi_zero / np.sqrt(k)))
            psi = ctx.cache("psi", lambda: psi_at_primes(ctx))[idx]
            yield (("psi(k)-k minus zero sum, /sqrt k", "cont", (psi - k - psi_zero + math.log(2 * math.pi)) / np.sqrt(k)))


def psi_at_primes(ctx):
    """psi(p) = sum of log q over prime powers q^m <= p, at every prime p in ctx.P."""
    P = ctx.P
    lp = ctx.lnP
    theta = np.cumsum(lp)
    # add prime powers q^m <= p (m >= 2): for each such power, add log q to all P >= q^m
    extra = np.zeros(len(P))
    qs = P[P <= int(math.isqrt(int(P[-1]))) + 1]
    pw_pos, pw_val = [], []
    for q in qs:
        q = int(q); v = q * q
        while v <= P[-1]:
            pw_pos.append(np.searchsorted(P, v)); pw_val.append(math.log(q)); v *= q
    add = np.zeros(len(P) + 1)
    np.add.at(add, np.array(pw_pos), np.array(pw_val))
    return theta + np.cumsum(add)[:len(P)]


def g_zeta_phases(ctx, idx):
    k = ctx.P[idx].astype(np.float64)
    lk = ctx.lnP[idx]
    for j in range(6):
        g = ctx.zeros[j]
        f = frac(g * lk / (2 * np.pi))
        yield ((f"frac(gamma_{j+1} ln k/2pi)", "cont", f))
        if j < 2:
            yield ((f"cos(gamma_{j+1} ln k)", "cont", np.cos(2 * np.pi * f)))
        with mpmath.workdps(40):
            a = mpmath.mpf(ctx.zeros_str[j]) * mpmath.log(2) / (2 * mpmath.pi)
        f2 = frac_mul(ctx.P[idx], a)
        yield ((f"frac(gamma_{j+1} k ln2/2pi) [phase at 2^k]", "cont", f2))
        if j < 2:
            yield ((f"cos(gamma_{j+1} k ln2) [phase at 2^k]", "cont", np.cos(2 * np.pi * f2)))


def li_arr(x):
    from scipy.special import expi
    return expi(np.log(x))


def riemann_R(x):
    from scipy.special import expi
    lx = np.log(x)
    mob = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1, 0, -1, 0, -1, 0, 1, 1, -1, 0, 0, 1, 0, 0, -1, -1]
    R = np.zeros_like(x)
    for n in range(1, 31):
        if mob[n] == 0:
            continue
        y = lx / n
        ok = y > math.log(1.5)
        R[ok] += mob[n] / n * expi(y[ok])
    return R


def g_residuals(ctx, idx):
    P = ctx.P
    k = P[idx].astype(np.float64)
    lk = ctx.lnP[idx]
    pik = (idx + 1).astype(np.float64)
    nrm = np.sqrt(k) / lk
    psi = ctx.cache("psi", lambda: psi_at_primes(ctx))[idx]
    yield (("(psi(k)-k)/sqrt k", "cont", (psi - k) / np.sqrt(k)))
    theta = np.cumsum(ctx.lnP)[idx]
    yield (("(theta(k)-k)/sqrt k", "cont", (theta - k) / np.sqrt(k)))
    yield (("(pi(k)-li(k))/(sqrt k/ln k)", "cont", (pik - li_arr(k)) / nrm))
    yield (("(pi(k)-R(k))/(sqrt k/ln k)", "cont", (pik - riemann_R(k)) / nrm))
    for m, A, Bc in ((4, [3], [1]), (3, [2], [1]), (8, [3], [1]), (8, [5], [1]), (8, [7], [1]),
                     (5, [2, 3], [1, 4]), (12, [5, 7, 11], [1])):
        r = np.mod(P, m)
        ca = np.cumsum(np.isin(r, A)); cb = np.cumsum(np.isin(r, Bc))
        lead = (ca[idx] - cb[idx] * (len(A) / len(Bc))) / nrm
        yield ((f"Chebyshev race mod {m}: {A} - {Bc}", "cont", lead))
        if m in (4, 3):
            yield ((f"Chebyshev race mod {m} lead > 0", "bin", (lead > 0).astype(float)))


# ---------------- automata, words, dynamics ----------------
def kolakoski(n):
    K = np.array([1, 2, 2], dtype=np.int8)
    while len(K) < n:
        m = min(len(K), n)  # use all known run lengths
        sym = np.where(np.arange(m) % 2 == 0, 1, 2).astype(np.int8)
        K = np.repeat(sym, K[:m])
    return K[:n]


def tribonacci_word(n):
    w = np.array([0], dtype=np.int8)
    while len(w) < n:
        # morphism 0->01, 1->02, 2->0
        lens = np.where(w == 2, 1, 2)
        out = np.zeros(lens.sum(), dtype=np.int8)
        pos = np.concatenate([[0], np.cumsum(lens)[:-1]])
        out[pos] = 0
        two = lens == 2
        out[pos[two] + 1] = np.where(w[two] == 0, 1, 2)
        w = out
    return w[:n]


def g_words(ctx, idx):
    k = ctx.P[idx]
    ku = k.astype(np.uint64)
    yield (("Thue-Morse t(k)", "bin", (popcount(ku) & 1).astype(float)))
    yield (("Rudin-Shapiro r(k)", "bin", (popcount(ku & (ku >> U64(1))) & 1).astype(float)))
    # Baum-Sweet: 1 iff no block of zeros of odd length in binary k
    x = k.copy(); run = np.zeros_like(x); bad = np.zeros(len(x), bool)
    nb = int(k.max()).bit_length()
    for b in range(nb + 1):
        bit = (x >> b) & 1
        top = (x >> b) == 0
        end = (bit == 1) | top
        bad |= end & (run % 2 == 1)
        run = np.where(bit == 0, run + 1, 0)
        run = np.where(top, 0, run)
    yield (("Baum-Sweet b(k)", "bin", (~bad).astype(float)))
    phi = (1 + 5 ** 0.5) / 2
    fphi = frac_mul(k, CONSTS["phi"])
    yield (("Fibonacci word bit s(k)", "bin", (fphi < 2 - phi).astype(float)))
    N = ctx.NMAX + 2
    K = ctx.cache("kol", lambda: kolakoski(N + 1001))
    kk = K[k - 1]                                   # 1-indexed: K(1)=1
    yield (("Kolakoski K(k)=2", "bin", (kk == 2).astype(float)))
    cs = ctx.cache("kolcs", lambda: np.concatenate([[0], np.cumsum(K == 2, dtype=np.int32)]))
    lo = np.maximum(k - 1 - 500, 0)
    yield (("Kolakoski density of 2s in k+-500", "cont", (cs[k - 1 + 501] - cs[lo]) / (k - 1 + 501 - lo)))
    yield (("Kolakoski discrepancy #2-#1 up to k /sqrt k", "cont", (2 * cs[k] - k) / np.sqrt(k)))
    T = ctx.cache("trib", lambda: tribonacci_word(N + 2))
    tk = T[k - 1]
    for c in range(3):
        yield ((f"Tribonacci word letter {c}", "bin", (tk == c).astype(float)))
    # Golay-Rudin-Shapiro partial sum up to k, /sqrt k (exact, chunked)
    rs_ps = ctx.cache("rsps", lambda: rs_partial(N + 2))
    yield (("Rudin-Shapiro partial sum /sqrt k", "cont", rs_ps[k] / np.sqrt(k)))


def rs_partial(n):
    out = np.zeros(n + 1, dtype=np.int32)
    acc = 0; CH = 1 << 22
    for s in range(0, n, CH):
        m = np.arange(s, min(s + CH, n), dtype=np.uint64)
        r = 1 - 2 * (popcount(m & (m >> U64(1))) & 1)
        c = np.cumsum(r) + acc
        out[s + 1:s + 1 + len(m)] = c
        acc = int(c[-1])
    return out                                      # out[k] = sum_{m<k} r(m)


def ca_raster(rule, W, n):
    rows = n // W + 2
    row = np.zeros(W, dtype=np.uint8); row[W // 2] = 1
    tab = np.array([(rule >> i) & 1 for i in range(8)], dtype=np.uint8)
    out = np.empty(rows * W, dtype=np.uint8)
    for t in range(rows):
        out[t * W:(t + 1) * W] = row
        row = tab[(np.roll(row, 1) << 2) | (row << 1) | np.roll(row, -1)]
    return out


def g_ca(ctx, idx):
    k = ctx.P[idx]
    for W in (1009, 3001):
        for rule in (30, 110, 90, 150, 45, 54, 73, 105, 22, 126):
            R = ctx.cache(("ca", rule, W), lambda: ca_raster(rule, W, ctx.NMAX + 3 * W))
            yield ((f"Rule {rule} raster bit, width {W}", "bin", R[k].astype(float)))
            if W == 1009 and rule in (30, 110):
                # density in the 5x5 space-time block around cell k
                acc = np.zeros(len(k))
                for dr in (-2, -1, 0, 1, 2):
                    for dc in (-2, -1, 0, 1, 2):
                        j = np.clip(k + dr * W + dc, 0, len(R) - 1)
                        acc += R[j]
                yield ((f"Rule {rule} 5x5 density, width {W}", "cont", acc / 25))


def g_dynamics(ctx, idx):
    k = ctx.P[idx]
    for nm in ("pi", "sqrt2", "e"):
        x = frac_mul(k, CONSTS[nm])
        cnt = np.full(len(k), 30.0); done = np.zeros(len(k), bool)
        for t in range(30):
            hit = (~done) & (x >= 0.4) & (x < 0.6)
            cnt[hit] = t; done |= hit
            x = 4 * x * (1 - x)
        yield ((f"logistic r=4 steps to [0.4,0.6) from frac(k*{nm})", "cont", cnt))
    for nm in ("pi", "e", "sqrt2", "phi"):
        f = frac_mul(k, CONSTS[nm])
        yield ((f"Gauss map first CF digit of frac(k*{nm})", "cont", np.floor(1 / np.maximum(f, 1e-18))))
    # Collatz 3x+1
    n = k.copy(); steps = np.zeros(len(k)); glide = np.full(len(k), -1.0); odd = np.zeros(len(k))
    mx = n.copy(); act = np.ones(len(k), bool); t = 0
    while act.any():
        t += 1
        a = np.flatnonzero(act)
        v = n[a]
        isodd = (v & 1) == 1
        v = np.where(isodd, 3 * v + 1, v >> 1)
        n[a] = v; odd[a] += isodd; mx[a] = np.maximum(mx[a], v)
        g0 = (glide[a] < 0) & (v < k[a]); glide[a[g0]] = t
        fin = v == 1
        steps[a[fin]] = t; act[a[fin]] = False
    yield (("Collatz total stopping time", "cont", steps))
    yield (("Collatz glide", "cont", glide))
    yield (("Collatz log2(max/k)", "cont", np.log2(mx / k)))
    yield (("Collatz odd-step fraction", "cont", odd / steps))
    yield (("Collatz total steps odd", "bin", (steps % 2).astype(float)))
    # 3x-1: which of the three known cycles (min element 1, 5, 17)
    n = k.copy()
    for t in range(1500):
        n = np.where(n & 1, 3 * n - 1, n >> 1)
    mn = n.copy()
    for t in range(25):
        n = np.where(n & 1, 3 * n - 1, n >> 1); mn = np.minimum(mn, n)
    for c in (1, 5, 17):
        yield ((f"3x-1 map ends in cycle {c}", "bin", (mn == c).astype(float)))
    # Mandelbrot escape count
    cre = -2 + 2.5 * frac_mul(k, CONSTS["sqrt2"]); cim = -1.25 + 2.5 * frac_mul(k, CONSTS["sqrt3"])
    zr = np.zeros(len(k)); zi = np.zeros(len(k)); esc = np.full(len(k), 60.0); live = np.ones(len(k), bool)
    for t in range(60):
        zr, zi = zr * zr - zi * zi + cre, 2 * zr * zi + cim
        e = live & (zr * zr + zi * zi > 4)
        esc[e] = t; live &= ~e
        zr[~live] = 0; zi[~live] = 0
    yield (("Mandelbrot escape at c(k)", "cont", esc))
    # Henon map escape
    x = 2 * frac_mul(k, CONSTS["pi"]) - 1; y = 0.4 * frac_mul(k, CONSTS["e"]) - 0.2
    esc = np.full(len(k), 50.0); live = np.ones(len(k), bool)
    for t in range(50):
        x, y = 1 - 1.4 * x * x + y, 0.3 * x
        e = live & (np.abs(x) > 10)
        esc[e] = t; live &= ~e
        x[~live] = 0; y[~live] = 0
    yield (("Henon escape from (frac k pi, frac k e)", "cont", esc))


# ---------------- geometry ----------------
def ulam(k):
    k = k.astype(np.int64)
    m = np.ceil((np.sqrt(k.astype(np.float64)) - 1) / 2).astype(np.int64)
    m = np.where((2 * m - 1) ** 2 >= k, m - 1, m); m = np.where((2 * m + 1) ** 2 < k, m + 1, m)
    t = k - (2 * m - 1) ** 2 - 1
    side = t // (2 * m); pos = t % (2 * m)
    x = np.select([side == 0, side == 1, side == 2], [m, m - 1 - pos, -m], -m + 1 + pos)
    y = np.select([side == 0, side == 1, side == 2], [-m + 1 + pos, m, m - 1 - pos], -m)
    return m, side, pos, x, y


def hilbert_d2xy(order, d):
    d = d.astype(np.int64).copy()
    x = np.zeros_like(d); y = np.zeros_like(d)
    s = 1
    while s < (1 << order):
        rx = 1 & (d // 2)
        ry = 1 & (d ^ rx)
        swap = ry == 0
        flip = swap & (rx == 1)
        x = np.where(flip, s - 1 - x, x); y = np.where(flip, s - 1 - y, y)
        x, y = np.where(swap, y, x), np.where(swap, x, y)
        x += s * rx; y += s * ry
        d //= 4; s *= 2
    return x, y


def powmod(b, e, m):
    b = np.asarray(b, dtype=np.int64) % m; e = np.asarray(e, dtype=np.int64).copy()
    r = np.ones_like(m)
    while np.any(e > 0):
        odd = (e & 1) == 1
        r = np.where(odd, (r * b) % m, r)
        b = (b * b) % m
        e >>= 1
    return r


def cornacchia(p, r, d):
    """x^2 + d y^2 = p from r^2 = -d mod p (vectorized). Returns x, y (NaN where it fails)."""
    p = p.astype(np.int64); r = np.where(r > p // 2, r, p - r)  # r in (p/2, p)
    a, b = p.copy(), r.copy()
    lim = np.floor(np.sqrt(p.astype(np.float64))).astype(np.int64)
    for _ in range(64):
        go = b > lim
        if not go.any():
            break
        a2 = np.where(go, b, a); b2 = np.where(go, a % np.maximum(b, 1), b)
        a, b = a2, b2
    x = b
    rem = p - x * x
    y2 = rem // d
    y = np.rint(np.sqrt(np.maximum(y2, 0).astype(np.float64))).astype(np.int64)
    ok = (rem % d == 0) & (y * y == y2)
    return np.where(ok, x, -1), np.where(ok, y, -1)


def sqrt_neg1(p):
    """c with c^2 = -1 mod p for p = 1 mod 4."""
    c = np.zeros_like(p); done = np.zeros(len(p), bool)
    for n in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71):
        cc = powmod(np.full_like(p, n), (p - 1) // 4, p)
        good = (~done) & ((cc * cc) % p == p - 1)
        c[good] = cc[good]; done |= good
        if done.all():
            break
    return c, done


def g_geometry(ctx, idx):
    k = ctx.P[idx]
    m, side, pos, x, y = ulam(k)
    yield (("Ulam ring index (monotone control)", "cont", m.astype(float)))
    yield (("Ulam ring parity", "bin", (m % 2).astype(float)))
    for s in range(4):
        yield ((f"Ulam side {s}", "bin", (side == s).astype(float)))
    yield (("Ulam position along side", "cont", pos / (2 * m)))
    ang = frac(np.arctan2(y, x) / (2 * np.pi))
    yield (("Ulam angle", "cont", ang))
    yield (("Ulam distance to a diagonal / ring", "cont", np.minimum(np.abs(x - y), np.abs(x + y)) / m))
    yield (("Ulam on a diagonal", "bin", (np.abs(x) == np.abs(y)).astype(float)))
    yield (("Ulam x parity", "bin", (x % 2).astype(float)))
    sq = np.sqrt(k.astype(np.float64))
    yield (("Sacks spiral angle frac(sqrt k)", "cont", frac(sq)))
    yield (("Sacks spiral cos(2pi sqrt k)", "cont", np.cos(2 * np.pi * sq)))
    # hexagonal spiral
    hm = np.ceil((-3 + np.sqrt(9 + 12 * (k - 1.0))) / 6).astype(np.int64)
    hm = np.where(3 * hm * (hm + 1) + 1 < k, hm + 1, hm); hm = np.where(3 * (hm - 1) * hm + 1 >= k, hm - 1, hm)
    t = k - (3 * hm * (hm - 1) + 1) - 1
    hs = t // hm
    for s in range(6):
        yield ((f"hex spiral side {s}", "bin", (hs == s).astype(float)))
    yield (("hex spiral position along side", "cont", (t % hm) / hm))
    n = np.floor((np.sqrt(8 * k + 1.0) - 1) / 2).astype(np.int64)
    n = np.where(n * (n + 1) // 2 > k, n - 1, n)
    yield (("triangle row position j/(n+1)", "cont", (k - n * (n + 1) // 2) / (n + 1)))
    hx, hy = hilbert_d2xy(14, k)
    yield (("Hilbert x/sqrt k", "cont", hx / sq)); yield (("Hilbert y/sqrt k", "cont", hy / sq))
    yield (("Hilbert angle", "cont", np.arctan2(hy, hx)))
    yield (("Hilbert x parity", "bin", (hx % 2).astype(float))); yield (("Hilbert y parity", "bin", (hy % 2).astype(float)))
    # Gaussian / Eisenstein / x^2+2y^2 angles
    p = k
    g1 = np.full(len(k), np.nan); g2 = np.full(len(k), np.nan)
    sel = np.flatnonzero(p % 4 == 1)
    if len(sel):
        c, ok = sqrt_neg1(p[sel])
        a, b = cornacchia(p[sel], c, 1)
        a = a.astype(float); b = b.astype(float)
        aa = np.where(a % 2 == 1, a, b); bb = np.where(a % 2 == 1, b, a)
        good = ok & (a > 0)
        g1[sel[good]] = (np.arctan2(bb, aa) / (np.pi / 2))[good]
        g2[sel[good]] = (np.arctan2(np.minimum(a, b), np.maximum(a, b)) / (np.pi / 4))[good]
    yield (("Gaussian prime angle (a odd)", "cont", g1))
    yield (("Gaussian prime octant angle", "cont", g2))
    e1 = np.full(len(k), np.nan)
    sel = np.flatnonzero((p % 3 == 1))
    if len(sel):
        ps = p[sel]; om = np.zeros_like(ps); done = np.zeros(len(ps), bool)
        for nn in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53):
            w = powmod(np.full_like(ps, nn), (ps - 1) // 3, ps)
            g = (~done) & (w != 1)
            om[g] = w[g]; done |= g
        r = (2 * om + 1) % ps
        x3, y3 = cornacchia(ps, r, 3)
        good = done & (x3 > 0)
        e1[sel[good]] = (np.arctan2(np.sqrt(3) * y3, x3) / (np.pi / 2))[good]
    yield (("Eisenstein angle (x^2+3y^2)", "cont", e1))
    q1 = np.full(len(k), np.nan)
    sel = np.flatnonzero(np.isin(p % 8, [1, 3]))
    if len(sel):
        ps = p[sel]; r = np.zeros_like(ps); done = np.zeros(len(ps), bool)
        m3 = ps % 8 == 3
        r[m3] = powmod(ps[m3] - 2, (ps[m3] + 1) // 4, ps[m3]); done[m3] = True
        m1 = np.flatnonzero(ps % 8 == 1)
        for nn in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67):
            if len(m1) == 0:
                break
            pp = ps[m1]
            z = powmod(np.full_like(pp, nn), (pp - 1) // 8, pp)
            z4 = powmod(z, np.full_like(pp, 4), pp)
            g = z4 == pp - 1
            rr = (z + powmod(z, np.full_like(pp, 3), pp)) % pp
            r[m1[g]] = rr[g]; done[m1[g]] = True
            m1 = m1[~g]
        x2, y2 = cornacchia(ps, r, 2)
        good = done & (x2 > 0)
        q1[sel[good]] = (np.arctan2(np.sqrt(2) * y2, x2) / (np.pi / 2))[good]
    yield (("x^2+2y^2 angle", "cont", q1))
    # Stern diatomic / Calkin-Wilf / Stern-Brocot
    a = np.ones(len(k), dtype=np.int64); b = np.ones(len(k), dtype=np.int64)
    nb = np.floor(np.log2(k.astype(np.float64))).astype(np.int64)
    for bit in range(int(nb.max()) - 1, -1, -1):
        act = nb > bit
        bt = (k >> bit) & 1
        na = np.where(bt == 0, a, a + b); nbb = np.where(bt == 0, a + b, b)
        a = np.where(act, na, a); b = np.where(act, nbb, b)
    yield (("Calkin-Wilf rational log(s(k)/s(k+1))", "cont", np.log(a / b)))
    yield (("Stern s(k): log s(k)/bits", "cont", np.log(a) / (nb + 1)))
    yield (("Stern s(k) odd", "bin", (a % 2).astype(float)))
    yield (("Stern s(k) = 0 mod 3", "bin", (a % 3 == 0).astype(float)))
    # Stern-Brocot BFS order: reverse the bits below the leading one
    rev = np.ones(len(k), dtype=np.int64)
    for bit in range(int(nb.max())):
        act = nb > bit
        rev = np.where(act, (rev << 1) | ((k >> bit) & 1), rev)
    a2 = np.ones(len(k), dtype=np.int64); b2 = np.ones(len(k), dtype=np.int64)
    for bit in range(int(nb.max()) - 1, -1, -1):
        act = nb > bit
        bt = (rev >> bit) & 1
        na = np.where(bt == 0, a2, a2 + b2); nbb = np.where(bt == 0, a2 + b2, b2)
        a2 = np.where(act, na, a2); b2 = np.where(act, nbb, b2)
    yield (("Stern-Brocot BFS rational log ratio", "cont", np.log(a2 / b2)))
    # Zeckendorf
    F = [1, 2]
    while F[-1] < 2 ** 40:
        F.append(F[-1] + F[-2])
    r = k.copy(); cnt = np.zeros(len(k)); low = np.zeros(len(k))
    for i in range(len(F) - 1, -1, -1):
        take = r >= F[i]
        r = np.where(take, r - F[i], r); cnt += take; low = np.where(take, i, low)
    yield (("Zeckendorf number of terms /log k", "cont", cnt / ctx.lnP[idx]))
    yield (("Zeckendorf terms odd", "bin", (cnt % 2).astype(float)))
    yield (("Zeckendorf lowest index even", "bin", (low % 2 == 0).astype(float)))


# ---------------- statistics: local primes ----------------
def g_gaps(ctx, idx):
    P = ctx.P; k = P[idx]; lk = ctx.lnP[idx]
    gm = k - P[idx - 1]; gp = P[idx + 1] - k
    gm2 = k - P[idx - 2]; gp2 = P[idx + 2] - k
    yield from [("gap before /ln k", "cont", gm / lk), ("gap after /ln k", "cont", gp / lk),
            ("gap after - gap before /ln k", "cont", (gp - gm) / lk),
            ("max neighbour gap /ln k", "cont", np.maximum(gm, gp) / lk),
            ("min neighbour gap /ln k", "cont", np.minimum(gm, gp) / lk),
            ("two gaps before /ln k", "cont", gm2 / lk), ("two gaps after /ln k", "cont", gp2 / lk),
            ("balanced prime (equal gaps)", "bin", (gm == gp).astype(float)),
            ("gap before = 0 mod 6", "bin", (gm % 6 == 0).astype(float)),
            ("gap after = 0 mod 6", "bin", (gp % 6 == 0).astype(float)),
            ("gap before = 4", "bin", (gm == 4).astype(float)), ("gap after = 4", "bin", (gp == 4).astype(float)),
            ("gap before = 6", "bin", (gm == 6).astype(float)), ("gap after = 6", "bin", (gp == 6).astype(float)),
            ("good prime p^2 > prev*next", "bin", (k * k > P[idx - 1] * P[idx + 1]).astype(float))]
    isp = lambda q: P[np.minimum(np.searchsorted(P, q), len(P) - 1)] == q
    yield from [("k-6 prime", "bin", isp(k - 6).astype(float)), ("k+6 prime", "bin", isp(k + 6).astype(float)),
            ("k-4 prime", "bin", isp(k - 4).astype(float)), ("k+4 prime", "bin", isp(k + 4).astype(float)),
            ("safe prime (k-1)/2 prime", "bin", isp((k - 1) // 2).astype(float))]
    for W in (100, 1000, 10000):
        c = np.searchsorted(P, k + W, side="right") - np.searchsorted(P, k - W, side="left")
        ex = li_arr((k + W).astype(float)) - li_arr(np.maximum(k - W, 2).astype(float))
        yield ((f"local prime count in k+-{W}, z-score", "cont", (c - ex) / np.sqrt(ex)))
    # rank of the gap after k among the 41 gaps around it
    G = np.diff(P)
    win = np.stack([G[np.clip(idx + j, 0, len(G) - 1)] for j in range(-20, 21)], axis=1)
    yield (("rank of gap after among 41 neighbours", "cont", (win < gp[:, None]).sum(1) / 41.0))
    pik = idx + 1
    yield from [("prime index pi(k) odd", "bin", (pik % 2).astype(float)),
            ("prime index pi(k) prime (super-prime)", "bin", isp(pik).astype(float)),
            ("prime index = 0 mod 3", "bin", (pik % 3 == 0).astype(float)),
            ("frac(pi(k)*phi)", "cont", frac_mul(pik, CONSTS["phi"]))]


# ---------------- misc: digits, constants' digits, hashes, words ----------------
WORDS1 = "zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen fifteen sixteen seventeen eighteen nineteen".split()
TENS = "  twenty thirty forty fifty sixty seventy eighty ninety".split(" ")


def say(n):
    def h(n):
        s = []
        if n >= 100:
            s.append(WORDS1[n // 100] + "hundred"); n %= 100
        if n >= 20:
            s.append(TENS[n // 10]); n %= 10
            if n:
                s.append(WORDS1[n])
        elif n > 0:
            s.append(WORDS1[n])
        return "".join(s)
    parts = []
    for div, nm in ((10**9, "billion"), (10**6, "million"), (1000, "thousand")):
        if n >= div:
            parts.append(h(n // div) + nm); n %= div
    if n:
        parts.append(h(n))
    return "".join(parts)


def const_bits(name, nbits):
    gmpy2.get_context().precision = nbits + 64
    if name == "pi":
        v = gmpy2.const_pi()
    elif name == "e":
        v = gmpy2.exp(1)
    elif name == "sqrt2":
        v = gmpy2.sqrt(2)
    elif name == "ln2":
        v = gmpy2.const_log2()
    elif name == "phi":
        v = (1 + gmpy2.sqrt(5)) / 2
    z = int(gmpy2.floor(gmpy2.frac(v) * gmpy2.mpz(2) ** nbits))
    by = z.to_bytes((nbits + 7) // 8, "big")
    bits = np.unpackbits(np.frombuffer(by, dtype=np.uint8))[-nbits:]
    gmpy2.get_context().precision = 53
    return bits                                       # bits[j] = (j+1)-th bit after the point


def pi_decimal(nd):
    gmpy2.get_context().precision = int(nd * 3.33) + 64
    z = gmpy2.floor(gmpy2.const_pi() * gmpy2.mpz(10) ** nd)
    s = gmpy2.mpz(z).digits(10)
    gmpy2.get_context().precision = 53
    return np.frombuffer(s[1:].encode(), dtype=np.uint8) - 48     # digit j = (j+1)-th decimal


def champernowne_digit(k):
    k = k.astype(np.int64) - 1                       # 0-indexed position in "123456789101112..."
    d = np.ones_like(k); cnt = 9 * np.ones_like(k); start = np.ones_like(k)
    for _ in range(12):
        over = k >= d * cnt
        k = np.where(over, k - d * cnt, k); d = np.where(over, d + 1, d)
        cnt = np.where(over, cnt * 10, cnt); start = np.where(over, start * 10, start)
    num = start + k // d; off = k % d
    return (num // (10 ** (d - 1 - off))) % 10


def ce_digits(ctx, nd):
    P = ctx.P
    s = "".join(map(str, P[:min(len(P), nd)].tolist()))
    assert len(s) >= nd, "need more primes for Copeland-Erdos"
    return np.frombuffer(s[:nd].encode(), dtype=np.uint8) - 48


def ce_bits(ctx, nb):
    P = ctx.P
    parts = []; tot = 0; i = 0
    while tot < nb:
        s = "".join(format(int(q), "b") for q in P[i:i + 200000]); parts.append(s); tot += len(s); i += 200000
    s = "".join(parts)[:nb]
    return np.frombuffer(s.encode(), dtype=np.uint8) - 48


SEG7 = np.array([6, 2, 5, 5, 4, 5, 6, 3, 7, 6])
DAH = np.array([5, 4, 3, 2, 1, 0, 1, 2, 3, 4])


def g_misc(ctx, idx):
    k = ctx.P[idx]
    for b in (3, 7, 10, 12):
        x = k.copy(); ds = np.zeros(len(k)); nd = np.zeros(len(k))
        while np.any(x > 0):
            ds += x % b; nd += (x > 0); x //= b
        mu = (b - 1) / 2; var = (b * b - 1) / 12
        yield ((f"base-{b} digit sum z-score", "cont", (ds - mu * nd) / np.sqrt(var * nd)))
    x = k.copy(); seg = np.zeros(len(k)); dah = np.zeros(len(k)); nd = np.zeros(len(k))
    while np.any(x > 0):
        d = x % 10; live = x > 0
        seg += np.where(live, SEG7[d], 0); dah += np.where(live, DAH[d], 0); nd += live; x //= 10
    yield (("seven-segment lit segments per digit", "cont", seg / nd))
    yield (("Morse dahs per digit", "cont", dah / nd))
    cd = champernowne_digit(k)
    yield (("Champernowne digit at k", "cont", cd.astype(float)))
    yield (("Champernowne digit odd", "bin", (cd % 2).astype(float)))
    CE = ctx.cache("ce", lambda: ce_digits(ctx, ctx.NMAX + 2))
    yield (("Copeland-Erdos digit at k", "cont", CE[k - 1].astype(float)))
    yield (("Copeland-Erdos digit odd", "bin", (CE[k - 1] % 2).astype(float)))
    CB = ctx.cache("ceb", lambda: ce_bits(ctx, ctx.NMAX + 2))
    yield (("binary Copeland-Erdos bit at k", "bin", CB[k - 1].astype(float)))
    for nm in ("pi", "e", "sqrt2", "ln2", "phi"):
        bits = ctx.cache(("bits", nm), lambda: const_bits(nm, ctx.NMAX + 8))
        yield ((f"bit k of {nm}", "bin", bits[k - 1].astype(float)))
    PD = ctx.cache("pidec", lambda: pi_decimal(ctx.NMAX + 2))
    yield (("decimal digit k of pi", "cont", PD[k - 1].astype(float)))
    yield (("decimal digit k of pi odd", "bin", (PD[k - 1] % 2).astype(float)))
    h = np.array([hashlib.sha256(str(int(v)).encode()).digest()[0] for v in k])
    yield (("sha256(str k) first byte [noise control]", "cont", h.astype(float)))
    yield (("sha256(str k) low bit [noise control]", "bin", (h & 1).astype(float)))
    h2 = np.array([hashlib.md5(str(int(v)).encode()).digest()[0] & 1 for v in k])
    yield (("md5(str k) low bit [noise control]", "bin", h2.astype(float)))
    letters = np.array([len(say(int(v))) for v in k], dtype=float)
    es = np.array([say(int(v)).count("e") for v in k], dtype=float)
    yield (("letters in the English name of k / log k", "cont", letters / ctx.lnP[idx]))
    yield (("share of 'e' in the English name of k", "cont", es / letters))


GROUPS = [
    ("zeta_sums", "zeta", "sample", g_zeta_sums),
    ("zeta_phases", "zeta", "full", g_zeta_phases),
    ("residuals", "zeta", "full", g_residuals),
    ("constants", "constants", "full", g_constants),
    ("quadratic", "constants", "full", g_quadratic),
    ("logs", "constants", "full", g_logs),
    ("words", "automata", "full", g_words),
    ("ca", "automata", "full", g_ca),
    ("dynamics", "dynamics", "full", g_dynamics),
    ("geometry", "geometry", "full", g_geometry),
    ("gaps", "statistics", "full", g_gaps),
    ("misc", "misc", "full", g_misc),
]
