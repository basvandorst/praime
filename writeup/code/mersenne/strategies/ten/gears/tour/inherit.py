"""Wheels that light each other (machine 23 of prime3.html).

Round 0: only the prime 2 is lit. In every round a dark number N looks at the lit primes that divide N - 1 or N + 1.
F1 is the part of N - 1 made of lit primes (with their full powers), F2 the same for N + 1, L = lcm(F1, F2).
If L * L > N, then N can be proved with those helpers alone:
  one cell   (Pocklington): a^(N-1) = 1 and gcd(a^((N-1)/q) - 1, N) = 1 for every prime q of F1  =>  every prime r | N is 1 mod F1
  two cells  (Morrison):    beta^(N+1) = 1 and beta^((N+1)/q) - 1 a unit for every prime q of F2  =>  every prime r | N is +1 or -1 mod F2
so a prime factor r <= sqrt N < L of N is 1 mod L (impossible) or equals x, the number below L that is 1 mod F1 and -1 mod F2:
one division by x finishes the proof.
rounds(N) = 1 + the smallest g such that the helpers lit in round g or before are enough.
Of those helpers the proof takes the largest prime powers first and stops as soon as L * L > N."""
from math import gcd, isqrt
from norm2 import jacobi


def factor(n, primes):
    f = []
    for p in primes:
        if p * p > n:
            break
        if n % p == 0:
            pe = 1
            while n % p == 0:
                n //= p
                pe *= p
            f.append((p, pe))
    if n > 1:
        f.append((n, n))
    return f


def helpers(N, f_minus, f_plus, depth, sides=2):
    """(g, F1, Q1, F2, Q2): the smallest round g whose lit primes are enough, and the parts they give."""
    items = [(depth[p], p, pe, 0) for p, pe in f_minus if p in depth]
    if sides == 2:
        items += [(depth[p], p, pe, 1) for p, pe in f_plus if p in depth]
    for g in sorted(set(i[0] for i in items)):
        F1 = F2 = L = 1
        Q1, Q2 = [], []
        for d, p, pe, side in sorted((i for i in items if i[0] <= g), key=lambda i: (-i[2], i[3])):
            G1, G2 = (F1, F2 * pe) if side else (F1 * pe, F2)
            M = G1 * G2 // gcd(G1, G2)
            if M == L:
                continue
            F1, F2, L = G1, G2, M
            (Q2 if side else Q1).append(p)
            if L * L > N:
                return g, F1, Q1, F2, Q2
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


def one_cell(N, Q1, amax=300):
    """Witnesses [(q, a)] or None. A base with a^(N-1) != 1 proves N composite: returns False."""
    out = []
    for q in Q1:
        for a in range(2, amax):
            if pow(a, N - 1, N) != 1:
                return False
            if gcd(pow(a, (N - 1) // q, N) - 1, N) == 1:
                out.append((q, a))
                break
        else:
            return None
    return out


def two_cells(N, Q2, amax=300):
    c = 2
    while True:
        j = jacobi(c, N)
        if j == -1:
            break
        if j == 0 and 1 < gcd(c, N) < N:
            return False
        c += 1
        if c > 2000:
            return None
    out = []
    for q in Q2:
        for a in range(1, amax):
            g = (a * a - c) % N
            if gcd(g, N) != 1:
                continue
            gi = pow(g, -1, N)
            u2 = qmul((a, 1), (a, 1), c, N)
            beta = (u2[0] * gi % N, u2[1] * gi % N)
            if qpow(beta, N + 1, c, N) != (1, 0):
                return False
            w = qpow(beta, (N + 1) // q, c, N)
            if gcd(((w[0] - 1) ** 2 - c * w[1] * w[1]) % N, N) == 1:
                out.append((q, a))
                break
        else:
            return None
    return c, out


def prove(N, F1, Q1, F2, Q2):
    """True when the helpers prove N prime; False when a test shows N is not prime; None when no witness was found."""
    L = F1 * F2 // gcd(F1, F2)
    assert L * L > N
    if Q1:
        r = one_cell(N, Q1)
        if not r:
            return r
    if Q2:
        r = two_cells(N, Q2)
        if not r:
            return r
    x = 1
    while (x + 1) % F2:
        x += F1
    if 1 < x < N and N % x == 0:
        return False
    return True


def run(top, sides=2, check=True):
    """Depth of every prime up to top, by the rule alone (primes come from the proofs, not from a sieve)."""
    depth = {2: 0}
    lit = [2]
    dark_passed = 0
    for N in range(3, top + 1, 2):
        h = helpers(N, factor(N - 1, lit), factor(N + 1, lit), depth, sides)
        if h is None:
            continue
        g, F1, Q1, F2, Q2 = h
        if prove(N, F1, Q1, F2, Q2) is True:
            depth[N] = g + 1
            lit.append(N)
    return depth


if __name__ == '__main__':
    import sys, json
    from sympy import isprime, primerange
    top = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    out = {}
    lines = []

    def say(s):
        lines.append(s)
        print(s)
        sys.stdout.flush()

    for sides, name in ((2, 'n - 1 and n + 1'), (1, 'n - 1 only')):
        d = run(top, sides)
        lit = sorted(d)
        ok = lit == list(primerange(2, top + 1))
        cnt = {}
        first = {}
        for N in lit:
            cnt[d[N]] = cnt.get(d[N], 0) + 1
            first.setdefault(d[N], N)
        say('%s, every odd number up to %d run through the proofs: the lit numbers are exactly the primes: %s; rounds: %s; first prime of each round: %s' % (name, top, ok, sorted(cnt.items()), sorted(first.items())))
        out[name] = dict(top=top, ok=ok, cnt=sorted(cnt.items()), first=sorted(first.items()))
        if sides == 2:
            for M in (30, 250, 5000):
                c = {}
                for N in lit:
                    if N <= M:
                        c[d[N]] = c.get(d[N], 0) + 1
                say('  up to %d: lit by round: %s' % (M, sorted(c.items())))
                out['upto_%d' % M] = sorted(c.items())
    # the theorem on composites: with all of N - 1 and N + 1 as helpers, no witness at all makes a composite pass
    d = run(3000, 2)
    tried = passed = 0
    for N in range(9, 3000, 2):
        if isprime(N):
            continue
        f1 = factor(N - 1, sorted(d))
        f2 = factor(N + 1, sorted(d))
        F1 = F2 = 1
        for p, pe in f1:
            F1 *= pe
        for p, pe in f2:
            F2 *= pe
        tried += 1
        if prove(N, F1, [p for p, pe in f1], F2, [p for p, pe in f2]) is True:
            passed += 1
    say('odd composites below 3000 with every prime factor of n - 1 and n + 1 as helper: %d tried, %d passed' % (tried, passed))
    out['composites'] = (tried, passed)
    json.dump(out, open('inherit.json', 'w'))
    open('inherit.txt', 'w').write('\n'.join(lines) + '\n')
