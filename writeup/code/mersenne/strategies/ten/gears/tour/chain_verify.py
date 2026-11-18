from math import gcd, isqrt

TD_LIMIT = 10**10


def is_small_prime(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def add(P, Q, a, N):
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
            raise ValueError('not a group')
        num, den = 3 * x1 * x1 + a, 2 * y1
    else:
        num, den = y2 - y1, x2 - x1
    den %= N
    if gcd(den, N) != 1:
        raise ValueError('not invertible')
    lam = num * pow(den, -1, N) % N
    x3 = (lam * lam - x1 - x2) % N
    return x3, (lam * (x1 - x3) - y1) % N


def mul(k, P, a, N):
    R = None
    Q = P
    while k:
        if k & 1:
            R = add(R, Q, a, N)
        Q = add(Q, Q, a, N)
        k >>= 1
    return R


def qmul(u, v, c, N):
    return (u[0] * v[0] + c * u[1] * v[1]) % N, (u[0] * v[1] + u[1] * v[0]) % N


def qpow(u, e, c, N):
    r = (1, 0)
    while e:
        if e & 1:
            r = qmul(r, u, c, N)
        u = qmul(u, u, c, N)
        e >>= 1
    return r


def verify(N, cert):
    try:
        return _verify(N, cert)
    except (ValueError, KeyError, IndexError, TypeError, ZeroDivisionError):
        return False


def _verify(N, cert):
    if not cert or cert[0]['N'] != N:
        return False
    proved_later = set()
    for st in reversed(cert):
        n = st['N']
        if n < 2:
            return False

        def known(f):
            return f in proved_later or (f < TD_LIMIT and is_small_prime(f))

        if st['t'] == 'td':
            if n >= TD_LIMIT or not is_small_prime(n):
                return False
        elif st['t'] == 'n-1':
            F = 1
            for f, a in st['f']:
                if not known(f) or (n - 1) % f:
                    return False
                if pow(a, n - 1, n) != 1 or gcd(pow(a, (n - 1) // f, n) - 1, n) != 1:
                    return False
                e = n - 1
                while e % f == 0:
                    e //= f
                    F *= f
            if len(set(f for f, a in st['f'])) != len(st['f']) or F * F <= n:
                return False
        elif st['t'] == 'n+1':
            c = st['c']
            if gcd(c, n) != 1:
                return False
            F = 1
            for f, a in st['f']:
                if not known(f) or (n + 1) % f:
                    return False
                g = (a * a - c) % n
                if gcd(g, n) != 1:
                    return False
                gi = pow(g, -1, n)
                u2 = qmul((a, 1), (a, 1), c, n)
                bt = (u2[0] * gi % n, u2[1] * gi % n)
                if qpow(bt, n + 1, c, n) != (1, 0):
                    return False
                w = qpow(bt, (n + 1) // f, c, n)
                if gcd((w[0] - 1) ** 2 - c * w[1] ** 2, n) != 1:
                    return False
                e = n + 1
                while e % f == 0:
                    e //= f
                    F *= f
            if len(set(f for f, a in st['f'])) != len(st['f']) or F < 2 or (F - 1) ** 2 <= n:
                return False
        elif st['t'] == 'ec':
            a, b, x, y, m, q = st['a'], st['b'], st['x'], st['y'], st['m'], st['q']
            if gcd(n, 6) != 1 or not known(q) or m % q:
                return False
            r = isqrt(q)
            if r < 2 or (r - 1) ** 4 <= n:
                return False
            if gcd(4 * a**3 + 27 * b * b, n) != 1:
                return False
            if (y * y - x**3 - a * x - b) % n:
                return False
            Q = mul(m // q, (x % n, y % n), a, n)
            if Q is None:
                return False
            if mul(q, Q, a, n) is not None:
                return False
        else:
            return False
        proved_later.add(n)
    return True
