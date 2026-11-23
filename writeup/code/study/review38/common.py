import numpy as np, glob, os, math
ZD = '/home/bas/praime/writeup/data/riemann/lzeros/'
PR = [2, 3, 5, 7, 11, 13]

def smoothN(t, q, const=0.0):
    t = np.asarray(t, float)
    return t / (2 * np.pi) * np.log(q * t / (2 * np.pi * np.e)) + const

def dS(l, p, L, K=60):
    k = np.arange(1, K + 1)
    return -(2 / np.pi**2) * np.sum(np.sin(np.pi * l * k * np.log(p) / L)**2 / (k**2 * float(p)**k))

def sigma2(u, lo, hi, ls, step=0.01, mode='dense', rng=None, nrand=None, center='emp'):
    """u: sorted unfolded zeros; window [x, x+l] with x in [lo, hi-l]."""
    out = []
    for l in ls:
        if mode == 'dense':
            x = np.arange(lo, hi - l, step)
        else:
            n = nrand or int((hi - l - lo) / step)
            x = rng.uniform(lo, hi - l, n)
        c = np.searchsorted(u, x + l) - np.searchsorted(u, x)
        if center == 'emp':
            out.append(c.var())
        else:
            out.append(np.mean((c - l)**2))
    return np.array(out)

def real_zeros():
    d = {}
    for f in ['r117_part0.npz', 'r117_part1.npz', 'r120_part0.npz', 'r120_part1.npz']:
        z = np.load(ZD + f)
        for k in z.keys():
            d[int(k[1:])] = np.sort(z[k])
    return d

def r122_zeros():
    d = {}
    for f in glob.glob(ZD + 'r122/D*.npy'):
        d[int(os.path.basename(f)[1:-4])] = np.sort(np.load(f))
    return d

def complex_zeros():
    d = {}
    for f in sorted(glob.glob(ZD + 'r121_part*.npz')):
        z = np.load(f)
        for k in z.keys():
            N, m = k[1:].split('_m')
            d[(int(N), int(m))] = np.sort(z[k])
    return d

def primes_of(n):
    n = abs(n); ps = []; p = 2
    while p * p <= n:
        if n % p == 0:
            ps.append(p)
            while n % p == 0: n //= p
        p += 1
    if n > 1: ps.append(n)
    return ps

def unfold_stats(g, q, tlo, thi, ls, frac=0.0, const=0.0, **kw):
    """g: positive heights; restrict range, optionally trimming frac at each end of the t range."""
    a = tlo + frac * (thi - tlo); b = thi - frac * (thi - tlo)
    gg = g[(g >= a) & (g <= b)]
    u = smoothN(gg, q, const)
    lo, hi = smoothN(a, q, const), smoothN(b, q, const)
    s = sigma2(u, lo, hi, ls, **kw)
    Lmean = np.mean(np.log(q * gg / (2 * np.pi)))
    return s, Lmean, len(gg)

def ols(X, y):
    X = np.asarray(X, float); y = np.asarray(y, float)
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    s2 = r @ r / (len(y) - X.shape[1])
    cov = s2 * np.linalg.inv(X.T @ X)
    return b, np.sqrt(np.diag(cov))
