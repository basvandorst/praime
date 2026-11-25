"""review143 own library: no imports from the repo code. Characters, theta, Gauss sums, zeros, X matrices, regression."""
import math, os
import numpy as np

ZD = __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/data/riemann/lzeros/'
FZ = __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/swap2/'


# ------------------------------------------------------------------ arithmetic
def kron(D, n):
    """Kronecker symbol (D/n), n >= 1 (own implementation)."""
    if math.gcd(D, n) != 1:
        return 0
    r = 1
    while n % 2 == 0:
        n //= 2
        if D % 8 in (3, 5): r = -r
    # Jacobi (D/n), n odd
    a = D % n
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: r = -r
        a %= n
    return r if n == 1 else 0


def factor(n):
    f = {}; p = 2
    while p * p <= n:
        while n % p == 0: f[p] = f.get(p, 0) + 1; n //= p
        p += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f


def kappa(q):
    r = 1.0
    for p in factor(q): r *= p / (p - 1)
    return r


def conrey_table(N, m):
    """values chi_N(m, n) for n = 0..N-1 (Conrey labelling), own construction (complex; 0 if gcd > 1)."""
    out = np.ones(N, complex)
    for p, e in factor(N).items():
        pe = p ** e; ph = pe - pe // p
        mm = m % pe
        if p == 2:
            if e == 1: continue
            # n = eps * 5^b mod 2^e ; eps in {1, -1}
            if e == 2:
                vals = {1: 1.0, 3: (-1.0 if mm == 3 else 1.0)}
                for n in range(N):
                    r = n % 4
                    out[n] *= vals.get(r, 0.0) if r in (1, 3) else 0.0
                continue
            logs = {}
            x = 1
            for b in range(pe // 4):
                logs[x] = (0, b); logs[(-x) % pe] = (1, b); x = x * 5 % pe
            em, bm = logs[mm]
            for n in range(N):
                r = n % pe
                if r % 2 == 0: out[n] = 0; continue
                en, bn = logs[r]
                out[n] *= np.exp(2j * math.pi * (em * en / 2 + bm * bn / 2 ** (e - 2)))
            continue
        # Conrey generator: smallest primitive root mod p that is primitive mod p^2 (if e >= 2)
        g = 2
        while True:
            ok = all(pow(g, (p - 1) // r, p) != 1 for r in factor(p - 1))
            if ok and e >= 2: ok = pow(g, p - 1, p * p) != 1
            if ok: break
            g += 1
        lg = {}; x = 1
        for k in range(ph): lg[x] = k; x = x * g % pe
        lm = lg[mm]
        for n in range(N):
            r = n % pe
            if r % p == 0: out[n] = 0; continue
            out[n] *= np.exp(2j * math.pi * lm * lg[r] / ph)
    for n in range(N):
        if math.gcd(n, N) != 1: out[n] = 0
    return out


def real_table(D):
    q = abs(D)
    return np.array([kron(D, n) if n > 0 else 0 for n in range(q)], float)


def gauss_eps(tab, a):
    q = len(tab); n = np.arange(q)
    tau = np.sum(tab * np.exp(2j * math.pi * n / q))
    return tau / (1j ** a * math.sqrt(q))


# ------------------------------------------------------------------ theta (own Stirling series, no scipy)
def im_lgamma(z):
    z = np.asarray(z, complex)
    lg = (z - 0.5) * np.log(z) - z + 0.5 * math.log(2 * math.pi) + 1 / (12 * z) - 1 / (360 * z ** 3) + 1 / (1260 * z ** 5) - 1 / (1680 * z ** 7)
    return lg.imag


def theta(t, q, a):
    t = np.asarray(t, float)
    return t / 2 * math.log(q / math.pi) + im_lgamma(0.25 + a / 2 + 0.5j * t)


# ------------------------------------------------------------------ zeros
def load_group(files):
    Z = {}
    for f in files:
        with np.load(ZD + f) as z:
            for k in z.files: Z[k] = np.array(z[k], float)
    return Z


R120F = ['r117_part0.npz', 'r117_part1.npz', 'r120_part0.npz', 'r120_part1.npz']
CPLXF = ['r121_part%d.npz' % i for i in range(4)]


# ------------------------------------------------------------------ X by direct double sum over zero pairs (own algorithm)
def X_direct(zi, phi_i, zj, phi_j, sig, W1, W2):
    """(1/T) sum e^{2i th_i(g1)} conj(e^{2i th_j(g2)}) K(g2 - g1) w(m), w(m) = [erf((W2-m)/s) - erf((W1-m)/s)]/2, s = sig/sqrt2."""
    from scipy.special import erf
    s = sig / math.sqrt(2); cut = 7 * sig
    lo = np.searchsorted(zj, zi - cut); hi = np.searchsorted(zj, zi + cut)
    cnt = hi - lo
    ii = np.repeat(np.arange(len(zi)), cnt)
    jj = np.repeat(lo - np.cumsum(np.r_[0, cnt[:-1]]), cnt) + np.arange(cnt.sum())
    u = zj[jj] - zi[ii]; m = (zj[jj] + zi[ii]) / 2
    w = 0.5 * (erf((W2 - m) / s) - erf((W1 - m) / s))
    return np.sum(phi_i[ii] * np.conj(phi_j[jj]) * np.exp(-u * u / (2 * sig * sig)) * w) / (W2 - W1)


# ------------------------------------------------------------------ regression: OLS with delete-one-unit jackknife (own)
def ols_jack(Acols, y, ui, uj, mask=None):
    A = np.column_stack(Acols).astype(complex if np.iscomplexobj(y) or any(np.iscomplexobj(c) for c in Acols) else float)
    y = np.asarray(y)
    if mask is not None: A, y, ui, uj = A[mask], y[mask], ui[mask], uj[mask]
    cpl = np.iscomplexobj(A) or np.iscomplexobj(y)
    if cpl: A = A.astype(complex); y = y.astype(complex)
    AH = A.conj()
    MM = AH.T @ A; vv = AH.T @ y
    c = np.linalg.solve(MM, vv)
    units = np.unique(np.r_[ui, uj])
    reps = []
    for f in units:
        sel = (ui == f) | (uj == f)
        As = A[sel]; ys = y[sel]
        reps.append(np.linalg.solve(MM - As.conj().T @ As, vv - As.conj().T @ ys))
    reps = np.array(reps); n = len(reps)
    se_r = np.sqrt((n - 1) / n * ((reps.real - reps.real.mean(0)) ** 2).sum(0))
    se_i = np.sqrt((n - 1) / n * ((reps.imag - reps.imag.mean(0)) ** 2).sum(0))
    return c, se_r, se_i, reps


def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]
