import json, os, random, sys, time
from math import gcd, isqrt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.setrecursionlimit(50000)
RHO = 20000
B = 100000
LIMIT = B * B
SIEVE = bytearray([1]) * (B + 1)
SIEVE[0] = SIEVE[1] = 0
for i in range(2, isqrt(B) + 1):
    if SIEVE[i]:
        SIEVE[i * i::i] = bytearray(len(SIEVE[i * i::i]))
SMALL = [i for i in range(B + 1) if SIEVE[i]]
DISCS = [(D, h, [int(c) for c in co]) for D, h, co in json.load(open(os.path.join(HERE, 'classpoly.json')))]


class Composite(Exception):
    pass


def small_prime(n):
    if n < 2:
        return False
    for p in SMALL:
        if p * p > n:
            return True
        if n % p == 0:
            return n == p
    return True


def trial(m):
    fs = []
    for p in SMALL:
        if p * p > m:
            break
        if m % p == 0:
            e = 0
            while m % p == 0:
                m //= p
                e += 1
            fs.append((p, e))
    if 1 < m < LIMIT:
        fs.append((m, 1))
        m = 1
    return fs, m


def sprp(n, a=2):
    if n < 3 or n % 2 == 0:
        return n == 2
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    x = pow(a, d, n)
    if x == 1 or x == n - 1:
        return True
    for _ in range(s - 1):
        x = x * x % n
        if x == n - 1:
            return True
    return False


def rho(n, iters=RHO):
    y = 2
    c = 1
    q = 1
    g = 1
    r = 1
    done = 0
    while g == 1 and done < iters:
        x = y
        for _ in range(r):
            y = (y * y + c) % n
        k = 0
        while k < r and g == 1 and done < iters:
            ys = y
            lim = min(100, r - k)
            for _ in range(lim):
                y = (y * y + c) % n
                q = q * abs(x - y) % n
            g = gcd(q, n)
            k += lim
            done += lim
        r *= 2
    if g == n:
        g = 1
        for _ in range(iters):
            ys = (ys * ys + c) % n
            g = gcd(abs(x - ys), n)
            if g > 1:
                break
    return g if 1 < g < n else None


def big_prime_part(c):
    for _ in range(4):
        if sprp(c):
            return c
        f = rho(c)
        if f is None:
            return None
        c = max(f, c // f)
    return c if sprp(c) else None


def jacobi(a, n):
    a %= n
    r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            r = -r
        a %= n
    return r if n == 1 else 0


def inv(x, n):
    x %= n
    g = gcd(x, n)
    if g != 1:
        raise Composite('a factor: %d' % g if g != n else 'zero divisor')
    return pow(x, -1, n)


def sqrt_mod(a, p):
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        r = pow(a, (p + 1) // 4, p)
    else:
        q = p - 1
        s = 0
        while q % 2 == 0:
            q //= 2
            s += 1
        z = 2
        while jacobi(z, p) != -1:
            z += 1
            if z > 2000:
                return None
        m = s
        c = pow(z, q, p)
        t = pow(a, q, p)
        r = pow(a, (q + 1) // 2, p)
        while t != 1:
            i = 0
            t2 = t
            while t2 != 1:
                t2 = t2 * t2 % p
                i += 1
                if i == m:
                    return None
            b = pow(c, 1 << (m - i - 1), p)
            m = i
            c = b * b % p
            t = t * c % p
            r = r * b % p
    return r if r * r % p == a else None


def cornacchia(D, N):
    if jacobi(D, N) != 1:
        return None
    x0 = sqrt_mod(D, N)
    if x0 is None:
        return None
    if (x0 - D) % 2:
        x0 = N - x0
    a, b = 2 * N, x0
    lim = isqrt(4 * N)
    while b > lim:
        a, b = b, a % b
    r = 4 * N - b * b
    if r % (-D):
        return None
    c = r // (-D)
    s = isqrt(c)
    if s * s != c:
        return None
    return b, s


def traces(D, t, s):
    if D == -3:
        return [t, -t, (t + 3 * s) // 2, -(t + 3 * s) // 2, (t - 3 * s) // 2, -(t - 3 * s) // 2]
    if D == -4:
        return [t, -t, 2 * s, -2 * s]
    return [t, -t]


def ec_add(P, Q, a, N):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % N == 0:
            return None
        if y1 != y2:
            raise Composite('two square roots')
        lam = (3 * x1 * x1 + a) * inv(2 * y1, N) % N
    else:
        lam = (y2 - y1) * inv(x2 - x1, N) % N
    x3 = (lam * lam - x1 - x2) % N
    return x3, (lam * (x1 - x3) - y1) % N


def ec_mul(k, P, a, N):
    R = None
    for bit in bin(k)[2:]:
        R = ec_add(R, R, a, N)
        if bit == '1':
            R = ec_add(R, P, a, N)
    return R


def pmulmod(u, v, H, N):
    h = len(H) - 1
    r = [0] * (2 * h - 1)
    for i, x in enumerate(u):
        if x:
            for j, y in enumerate(v):
                r[i + j] = (r[i + j] + x * y) % N
    for i in range(2 * h - 2, h - 1, -1):
        c = r[i]
        if c:
            for j in range(h + 1):
                r[i - h + j] = (r[i - h + j] - c * H[j]) % N
    return r[:h]


def ppowmod(base, e, H, N):
    h = len(H) - 1
    r = [1] + [0] * (h - 1)
    for bit in bin(e)[2:]:
        r = pmulmod(r, r, H, N)
        if bit == '1':
            r = pmulmod(r, base, H, N)
    return r


def pgcd(u, v, N):
    def trim(w):
        w = [x % N for x in w]
        while w and w[-1] == 0:
            w.pop()
        return w
    u, v = trim(u), trim(v)
    while v:
        iv = inv(v[-1], N)
        v = [x * iv % N for x in v]
        while len(u) >= len(v):
            c = u[-1]
            k = len(u) - len(v)
            for j in range(len(v)):
                u[k + j] = (u[k + j] - c * v[j]) % N
            u = trim(u)
        u, v = v, u
    return u


def class_root(co, N, rng):
    H = [c % N for c in reversed(co)]
    h = len(H) - 1
    if h == 1:
        return (-H[0]) % N
    for _ in range(40):
        if h == 1:
            return (-H[0]) * inv(H[1], N) % N
        base = [rng.randrange(N), 1] + [0] * (h - 2)
        w = ppowmod(base, (N - 1) // 2, H, N)
        w[0] = (w[0] - 1) % N
        g = pgcd(H[:], w, N)
        if 1 < len(g) < len(H):
            if 2 * (len(g) - 1) > h:
                q = []
                u = H[:]
                while len(u) >= len(g):
                    c = u[-1]
                    q.append(c)
                    k = len(u) - len(g)
                    for j in range(len(g)):
                        u[k + j] = (u[k + j] - c * g[j]) % N
                    u.pop()
                g = q[::-1]
                iv = inv(g[-1], N)
                g = [x * iv % N for x in g]
            H = g
            h = len(H) - 1
    return None


def find_curve(D, co, m, q, N, rng):
    j = None
    if D not in (-3, -4):
        j = class_root(co, N, rng)
        if j is None or j in (0, 1728 % N):
            return None
        c = j * inv(j - 1728, N) % N
        a0, b0 = (-3 * c) % N, (2 * c) % N
    for _ in range(60):
        x = rng.randrange(1, N)
        if D == -3:
            y = rng.randrange(1, N)
            a, b = 0, (y * y - x * x * x) % N
        elif D == -4:
            y = rng.randrange(1, N)
            a, b = (y * y - x * x * x) * inv(x, N) % N, 0
        else:
            z = (x * x * x + a0 * x + b0) % N
            if z == 0:
                continue
            a, b = a0 * z * z % N, b0 * z * z * z % N
            x, y = z * x % N, z * z % N
        if gcd(4 * a * a * a + 27 * b * b, N) != 1:
            continue
        P = (x, y)
        Q = ec_mul(m // q, P, a, N)
        if Q is None:
            continue
        if ec_mul(q, Q, a, N) is None:
            return {'t': 'ec', 'N': N, 'D': D, 'a': a, 'b': b, 'x': x, 'y': y, 'm': m, 'q': q}
    return None


def qmul(u, v, c, N):
    return (u[0] * v[0] + c * u[1] * v[1]) % N, (u[0] * v[1] + u[1] * v[0]) % N


def qpow(u, e, c, N):
    r = (1, 0)
    for bit in bin(e)[2:]:
        r = qmul(r, r, c, N)
        if bit == '1':
            r = qmul(r, u, c, N)
    return r


def beta(a, c, N):
    g = inv(a * a - c, N)
    u2 = qmul((a, 1), (a, 1), c, N)
    return u2[0] * g % N, u2[1] * g % N


def step_minus(N, fs):
    out = []
    for f in fs:
        for a in range(2, 300):
            if pow(a, N - 1, N) != 1:
                raise Composite('Fermat base %d' % a)
            if gcd(pow(a, (N - 1) // f, N) - 1, N) == 1:
                out.append([f, a])
                break
        else:
            return None
    return {'t': 'n-1', 'N': N, 'f': out}


def step_plus(N, fs):
    c = 2
    while jacobi(c, N) != -1:
        if jacobi(c, N) == 0 and gcd(c, N) not in (1, N):
            raise Composite('a factor')
        c += 1
        if c > 2000:
            return None
    out = []
    for f in fs:
        for a in range(1, 300):
            if gcd(a * a - c, N) != 1:
                continue
            bt = beta(a, c, N)
            if qpow(bt, N + 1, c, N) != (1, 0):
                raise Composite('two cells, start %d' % a)
            w = qpow(bt, (N + 1) // f, c, N)
            if gcd((w[0] - 1) ** 2 - c * w[1] ** 2, N) == 1:
                out.append([f, a])
                break
        else:
            return None
    return {'t': 'n+1', 'N': N, 'c': c, 'f': out}


class Prover:
    def __init__(self, hmax=4, seed=1, use_pm=True, use_ec=True, stop=LIMIT):
        self.stop = stop
        self.discs = [d for d in DISCS if d[1] <= hmax]
        self.rng = random.Random(seed)
        self.memo = {}
        self.use_pm = use_pm
        self.use_ec = use_ec
        self.nodes = 0

    def candidates(self, N, deep=False):
        cands = []
        if self.use_pm:
            for sign, kind in ((-1, 'n-1'), (1, 'n+1')):
                fs, cof = trial(N + sign)
                F = 1
                for p, e in fs:
                    F *= p ** e
                ok = F * F > N if sign == -1 else F > 1 and (F - 1) ** 2 > N
                if ok:
                    if deep:
                        continue
                    part = []
                    G = 1
                    for p, e in sorted(fs, reverse=True):
                        part.append(p)
                        G *= p ** e
                        if (G * G > N) if sign == -1 else (G - 1) ** 2 > N:
                            break
                    cands.append((0, kind, part, None, None))
                elif cof > 1:
                    if deep == sprp(cof):
                        continue
                    if deep:
                        cof = big_prime_part(cof)
                        if cof is None or cof * cof <= 4 * N:
                            continue
                    if cof * cof > N if sign == -1 else (cof - 1) ** 2 > N:
                        cands.append((cof, kind, [cof], cof, None))
                    else:
                        cands.append((cof, kind, [p for p, e in fs] + [cof], cof, None))
        if self.use_ec and N % 2 and N % 3:
            for D, h, co in self.discs:
                ts = cornacchia(D, N)
                if ts is None:
                    continue
                for tr in traces(D, *ts):
                    m = N + 1 - tr
                    fs, cof = trial(m)
                    q = cof if cof > 1 else (fs[-1][0] if fs else 1)
                    if q < 2:
                        continue
                    r = isqrt(q)
                    if r < 2 or (r - 1) ** 4 <= N:
                        continue
                    if deep:
                        if cof < 2 or sprp(cof):
                            continue
                        q = big_prime_part(cof)
                        if q is None:
                            continue
                        r = isqrt(q)
                        if (r - 1) ** 4 <= N:
                            continue
                    elif cof > 1 and not sprp(cof):
                        continue
                    cands.append((q, 'ec', (D, co, m), q, None))
        rank = {'n-1': 0, 'n+1': 1, 'ec': 2}
        cands.sort(key=lambda c: (c[0], rank[c[1]]))
        return cands

    def prove(self, N):
        if N < self.stop:
            return [{'t': 'td', 'N': N}] if small_prime(N) else None
        if N in self.memo:
            return self.memo[N]
        self.nodes += 1
        res = None
        for deep in (False, True) if sprp(N) else ():
            if res is not None:
                break
            for q0, kind, info, q, _ in self.candidates(N, deep):
                sub = []
                if q is not None and q >= self.stop:
                    sub = self.prove(q)
                    if sub is None:
                        continue
                if kind == 'n-1':
                    st = step_minus(N, info)
                elif kind == 'n+1':
                    st = step_plus(N, info)
                else:
                    D, co, m = info
                    st = find_curve(D, co, m, q, N, self.rng)
                if st is None:
                    continue
                res = [st] + sub
                break
        self.memo[N] = res
        return res


def certify(N, hmax=4, seed=1, use_pm=True, use_ec=True, stop=LIMIT):
    P = Prover(hmax, seed, use_pm, use_ec, stop)
    try:
        return P.prove(N), P.nodes
    except Composite as e:
        return 'composite: ' + str(e), P.nodes


if __name__ == '__main__':
    from chain_verify import verify
    for arg in sys.argv[1:]:
        N = eval(arg.replace('^', '**'))
        t = time.time()
        cert, nodes = certify(N)
        dt = time.time() - t
        if isinstance(cert, list):
            print(arg, ': proved,', len(cert), 'steps', [s['t'] for s in cert], '%.2f s' % dt, 'checker:', verify(N, cert))
        else:
            print(arg, ':', cert, '%.2f s' % dt)
