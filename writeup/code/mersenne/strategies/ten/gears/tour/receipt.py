"""A Mersenne test with a receipt (Pietrzak's halving proof for repeated squaring), in two groups:
one cell (x -> x^2 mod N from x = 3: the probable-prime test 3^(2^p) = 9) and two cells
(a + b sqrt 3 -> its square from 2 + sqrt 3: the Lucas-Lehmer test in the form used by machine 11)."""
import hashlib
import gmpy2
from gmpy2 import mpz

RBITS = 64


class One:
    def __init__(self, p):
        self.p = p
        self.N = (mpz(1) << p) - 1
        self.x0 = mpz(3) % self.N
        self.T = p

    def red(self, v):
        p, N = self.p, self.N
        v = (v & N) + (v >> p)
        v = (v & N) + (v >> p)
        return v - N if v >= N else v

    def sq(self, v):
        return self.red(v * v)

    def mul(self, u, v):
        return self.red(u * v)

    def one(self):
        return mpz(1)

    def ser(self, v):
        return v.digits(16)

    def verdict(self, y):
        return y == 9 % self.N


class Two(One):
    def __init__(self, p):
        One.__init__(self, p)
        self.x0 = (mpz(2), mpz(1))
        self.T = p - 2

    def sq(self, v):
        a, b = v
        return (self.red(2 * a * a + self.N - 1), self.red(2 * a * b))

    def mul(self, u, v):
        m0 = u[0] * v[0]
        m1 = u[1] * v[1]
        return (self.red(m0 + 3 * m1), self.red((u[0] + u[1]) * (v[0] + v[1]) + 2 * self.N * self.N - m0 - m1))

    def one(self):
        return (mpz(1), mpz(0))

    def ser(self, v):
        return v[0].digits(16) + ',' + v[1].digits(16)

    def verdict(self, y):
        return y[0] == 0


def power(G, v, e):
    r = G.one()
    for bit in bin(e)[2:]:
        r = G.sq(r)
        if bit == '1':
            r = G.mul(r, v)
    return r


def forge(G, k, y, saved, eps, d, nonce=0):
    """Try to build a receipt for the wrong result y * eps, where eps has small order d. Needs the true chain.
    Returns (y_wrong, receipt) or None when the hash did not fall right."""
    hs, _ = halves(G.T, k)
    ystar = G.mul(y, eps)
    x = G.x0
    T = G.T
    Y, Ystar = y, ystar
    e_eff = eps
    if T % 2:
        Y, Ystar, e_eff = G.sq(Y), G.sq(Ystar), G.sq(eps)
        T += 1
    mu = saved[hs[0]]
    pw = [G.one()]
    for _ in range(d - 1):
        pw.append(G.mul(pw[-1], e_eff))
    hit = None
    for t in range(d):
        mus = G.mul(mu, pw[t])
        r = challenge(G, x, Ystar, mus, T, nonce)
        if pw[(t * ((pow(2, hs[0], d) - r) % d)) % d] == e_eff and (Ystar != Y or t):
            hit = (t, r, mus)
            break
    if hit is None:
        return None
    t, r, mus = hit
    delta = pw[t]
    rs = [r]
    out = [mus]
    memo = {}

    def f(j, m):
        if j == 0:
            return saved[m]
        key = (j, m)
        if key not in memo:
            memo[key] = G.mul(power(G, f(j - 1, m), rs[j - 1]), f(j - 1, m + hs[j - 1]))
        return memo[key]

    x = G.mul(power(G, x, r), mus)
    yy = G.mul(power(G, mus, r), Ystar)
    T //= 2
    c = 1
    for i in range(1, k):
        if T % 2:
            yy = G.sq(yy)
            T += 1
        m = G.mul(f(i, hs[i]), power(G, delta, (c * pow(2, hs[i], d)) % d))
        r = challenge(G, x, yy, m, T, nonce)
        rs.append(r)
        out.append(m)
        x = G.mul(power(G, x, r), m)
        yy = G.mul(power(G, m, r), yy)
        c = (c * (r + pow(2, hs[i], d))) % d
        T //= 2
    return ystar, out


def halves(T, k):
    hs = []
    for _ in range(k):
        if T % 2:
            T += 1
        T //= 2
        hs.append(T)
    return hs, T


def needed(hs):
    idx = set()
    sums = [0]
    for h in hs:
        for s in sums:
            idx.add(s + h)
        sums = sums + [s + h for s in sums]
    return idx


def challenge(G, x, y, mu, T, nonce=0):
    h = hashlib.sha256(('%d|%d|%d|%s|%s|%s' % (G.p, T, nonce, G.ser(x), G.ser(y), G.ser(mu))).encode()).digest()
    return int.from_bytes(h[:RBITS // 8], 'big')


def run(G, k, flip=None):
    """The whole test: T squarings, keeping the residues a receipt of k halvings needs. Returns (y, saved, squarings)."""
    hs, _ = halves(G.T, k)
    idx = needed(hs)
    top = max(idx | {G.T})
    saved = {}
    v = G.x0
    y = None
    for i in range(1, top + 1):
        v = G.sq(v)
        if flip is not None and i == flip[0]:
            v = v ^ (mpz(1) << flip[1]) if not isinstance(v, tuple) else (v[0] ^ (mpz(1) << flip[1]), v[1])
        if i in idx:
            saved[i] = v
        if i == G.T:
            y = v
    return y, saved, top


def build(G, k, y, saved, nonce=0, count=None):
    """The receipt: one middle value per halving, each a product of saved residues."""
    hs, _ = halves(G.T, k)
    rs = []
    memo = {}

    def f(j, m):
        if j == 0:
            return saved[m]
        key = (j, m)
        if key not in memo:
            if count is not None:
                count[0] += 1
            memo[key] = G.mul(power(G, f(j - 1, m), rs[j - 1]), f(j - 1, m + hs[j - 1]))
        return memo[key]

    x = G.x0
    T = G.T
    mus = []
    for i in range(k):
        if T % 2:
            y = G.sq(y)
            T += 1
        mu = f(i, hs[i])
        r = challenge(G, x, y, mu, T, nonce)
        rs.append(r)
        mus.append(mu)
        x = G.mul(power(G, x, r), mu)
        y = G.mul(power(G, mu, r), y)
        T //= 2
    return mus


def check(G, y, mus, nonce=0, count=None):
    """The check: k halvings and then only the last T / 2^k squarings."""
    x = G.x0
    T = G.T
    for mu in mus:
        if T % 2:
            y = G.sq(y)
            T += 1
        r = challenge(G, x, y, mu, T, nonce)
        x = G.mul(power(G, x, r), mu)
        y = G.mul(power(G, mu, r), y)
        T //= 2
    for _ in range(T):
        x = G.sq(x)
    if count is not None:
        count[0] = T
    return x == y


def lucas_lehmer(p):
    if p == 2:
        return True
    N = (mpz(1) << p) - 1
    s = mpz(4)
    for _ in range(p - 2):
        s = s * s - 2
        s = (s & N) + (s >> p)
        if s >= N:
            s -= N
    return s == 0


if __name__ == '__main__':
    import sys, time
    for arg in sys.argv[1:]:
        p = int(arg)
        for cls, name in ((One, 'one cell'), (Two, 'two cells')):
            G = cls(p)
            k = 6
            t = time.time()
            y, saved, top = run(G, k)
            t1 = time.time() - t
            t = time.time()
            mus = build(G, k, y, saved)
            t2 = time.time() - t
            t = time.time()
            ok = check(G, y, mus)
            t3 = time.time() - t
            print(p, name, 'verdict', G.verdict(y), 'receipt accepted', ok, 'test %.2f s, receipt %.2f s, check %.3f s' % (t1, t2, t3))
