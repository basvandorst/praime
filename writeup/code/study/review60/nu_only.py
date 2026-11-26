def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]
PR = primes_upto(2_000_000).astype(float)
u = 1.0 / ((PR - 1) * (PR ** 2 - 1))
C0 = float(np.exp(np.log1p(-u).sum()))

def nu_l(l, e):
    if e == 0: return 1 - 1.0 / ((l - 1) * (l * l - 1))
    v = 1 + 1.0 / (l * l - 1) + (l / (l + 1.0)) * sum(l ** (-j) for j in range(1, e // 2 + 1))
    if e % 2: v -= l ** (-(e - 1) / 2) / (l * l - 1.0)
    return v

def fac(m):
    f = {}; d = 2
    while d * d <= m:
        while m % d == 0: f[d] = f.get(d, 0) + 1; m //= d
        d += 1
    if m > 1: f[m] = f.get(m, 0) + 1
    return f

def nu(m):
    v = C0
    for l, e in fac(m).items(): v *= nu_l(l, e) / nu_l(l, 0)
    return v

