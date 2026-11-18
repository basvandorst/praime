"""Three machines that find primes before their own tick (machines 20, 21 and 22 of prime3.html).

ladder(n)      the two-cell wheel of machine 19, turned by squaring: one tick is "square, and step once more
               when the next binary digit of n is 1". Wheel n has decided after bit_length(n) - 1 ticks.
self_feed(top) a sieve that takes its gears from its own lamps: tick j adds the j-th prime gear; with the gears
               p_1 .. p_j every number below (p_(j+1))^2 is decided.
lucy(x)        counters instead of lamps: one counter for every value floor(x / k); prime gear p takes from every
               counter v >= p*p the numbers whose smallest factor is p. At the end counter v reads pi(v)."""
import sys
from math import isqrt
from norm2 import start


def ladder(n):
    """(ticks, verdict): verdict True (home), False (not home), or the reason no wheel turns."""
    if n < 2:
        return 0, False
    if n % 2 == 0:
        return 0, n == 2
    s = start(n)
    if s is None:
        return 0, False
    if s[0] == 'factor':
        return 0, False
    a, c = s
    x, y = a % n, 1
    ticks = 0
    for bit in bin(n)[3:]:
        x, y = (x * x + c * y * y) % n, 2 * x * y % n
        if bit == '1':
            x, y = (a * x + c * y) % n, (x + a * y) % n
        ticks += 1
    return ticks, (x == a % n and y == n - 1)


def self_feed(top):
    """Yields (tick, gear, next gear, new lamps) while the frontier is at most top. Uses nothing but its own lamps."""
    struck = bytearray(top + 1)
    lit = bytearray(top + 1)
    gear = 1
    tick = 0
    done = 2
    while True:
        gear += 1
        while struck[gear]:
            gear += 1
        tick += 1
        for m in range(gear * gear, top + 1, gear):
            struck[m] = 1
        nxt = gear + 1
        while struck[nxt]:
            nxt += 1
        if nxt * nxt > top + 1:
            return
        new = [m for m in range(done, nxt * nxt) if not struck[m]]
        for m in new:
            lit[m] = 1
        done = nxt * nxt
        yield tick, gear, nxt, new


def lucy(x):
    """Returns (V, C, moves, primes): V the counter values (descending), C the final counters (pi of each value),
    moves the number of counter changes, primes the gears whose counter rose (found by counting only)."""
    r = isqrt(x)
    V = [x // i for i in range(1, r + 1)]
    V += list(range(V[-1] - 1, 0, -1))
    C = [v - 1 for v in V]
    n = len(V)

    def idx(v):
        return n - v if v <= r else x // v - 1

    moves = 0
    primes = []
    for p in range(2, r + 1):
        if C[idx(p)] == C[idx(p - 1)]:
            continue
        primes.append(p)
        cp = C[idx(p - 1)]
        p2 = p * p
        for i in range(n):
            v = V[i]
            if v < p2:
                break
            C[i] -= C[idx(v // p)] - cp
            moves += 1
    return V, C, moves, primes


if __name__ == '__main__':
    from sympy import isprime, primepi, primerange
    import json
    out = {}
    lines = []

    def say(s):
        lines.append(s)
        print(s)
        sys.stdout.flush()

    # machine 20
    from norm2 import home
    top = 300000
    bad = badt = 0
    by_tick = {}
    for n in range(2, top):
        t, v = ladder(n)
        h = home(n)
        if (v is True) != (h is True):
            bad += 1
        if n % 2 and v is not None and start(n) and start(n)[0] != 'factor' and t != n.bit_length() - 1:
            badt += 1
        if v is True:
            by_tick[n.bit_length() - 1] = by_tick.get(n.bit_length() - 1, 0) + 1
        if (v is True) != isprime(n):
            bad += 1
    say('squaring wheel, every n from 2 to %d: verdicts that differ from the stepping wheel or from the primes: %d; wheels whose number of ticks is not (binary digits - 1): %d' % (top - 1, bad, badt))
    cum = 0
    rows = []
    for T in range(1, 18):
        cum += by_tick.get(T, 0)
        rows.append((T, 2 ** (T + 1) - 1, cum, int(primepi(2 ** (T + 1) - 1))))
    assert all(r[2] == r[3] for r in rows)
    out['ladder'] = dict(top=top - 1, bad=bad, badt=badt, rows=rows)
    say('lamps lit by tick T (all wheels below 2^(T+1)), T = 1..17: ' + ', '.join(str(r[2]) for r in rows) + ' = pi(2^(T+1) - 1) for every T')

    # machine 21
    TOP = 12000000
    isp = bytearray([1]) * (TOP + 1)
    isp[0] = isp[1] = 0
    for i in range(2, isqrt(TOP) + 1):
        if isp[i]:
            isp[i * i::i] = bytearray(len(isp[i * i::i]))
    total = 0
    sf = []
    wrong = 0
    for tick, gear, nxt, new in self_feed(TOP):
        total += len(new)
        if any(not isp[m] for m in new) or total != sum(isp[:nxt * nxt]):
            wrong += 1
        sf.append((tick, gear, nxt, nxt * nxt, len(new), total))
    assert wrong == 0
    assert [r[1] for r in sf] == list(primerange(2, sf[-1][1] + 1))
    out['self'] = dict(top=TOP, ticks=len(sf), rows=sf[:40], last=sf[-1], wrong=wrong,
                       min_new=min((r[4], r[0]) for r in sf[1:]), fewest=sorted((r[4], r[0], r[1]) for r in sf[1:])[:6])
    say('self-feeding sieve up to %d: %d ticks, every lamp a prime and every prime below the frontier lit at every tick; last tick %d: gear %d, frontier %d, lamps %d' % (TOP, len(sf), sf[-1][0], sf[-1][1], sf[-1][3], sf[-1][5]))
    say('new lamps at ticks 1..20: ' + ', '.join(str(r[4]) for r in sf[:20]))
    say('fewest new lamps at a tick after the first: %d (tick %d); the six smallest (lamps, tick, gear): %s' % (out['self']['min_new'] + (out['self']['fewest'],)))
    for N in (30, 250, 5000, 100000, 10 ** 7):
        t = next(r for r in sf if r[3] > N)
        say('all primes up to %d are lit at tick %d (gear %d, frontier %d): %d primes' % (N, t[0], t[1], t[3], int(primepi(N))))
        out['self']['all_%d' % N] = (t[0], t[1], t[3], int(primepi(N)))
    for T in (10, 30, 100):
        r = sf[T - 1]
        say('at tick %d: gear %d, every number below %d decided, %d lamps' % (T, r[1], r[3], r[5]))
        out['self']['tick_%d' % T] = r

    # machine 22
    cnt = []
    for x in (30 ** 2, 100 ** 2, 250 ** 2, 1000 ** 2, 5000 ** 2, 10 ** 8, 10 ** 9):
        V, C, moves, primes = lucy(x)
        assert primes == list(primerange(2, isqrt(x) + 1))
        okc = all(C[i] == int(primepi(V[i])) for i in range(0, len(V), max(1, len(V) // 400)))
        assert okc and C[0] == int(primepi(x))
        cnt.append((x, len(V), len(primes), moves, C[0]))
        say('counting wheels, x = %d: %d counters, %d prime gears, %d moves, pi(x) = %d (moves per prime %.3f); gears found by a risen counter = the primes up to sqrt x; sampled counters equal pi' % (x, len(V), len(primes), moves, C[0], moves / C[0]))
    out['count'] = cnt
    cross = None
    prev = None
    for N in range(2, 3000):
        V, C, moves, primes = lucy(N * N)
        if moves < C[0]:
            if prev is None:
                prev = N
        else:
            prev = None
    out['cross'] = prev
    say('counting wheels with N wheels (x = N^2): the moves stay below the number of primes counted for every N from %d on (tried up to 2999)' % prev)
    json.dump(out, open('early.json', 'w'))
    open('early.txt', 'w').write('\n'.join(lines) + '\n')
