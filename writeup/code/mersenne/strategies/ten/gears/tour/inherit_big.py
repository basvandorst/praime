"""Rounds needed by large primes on machine 23 (wheels that light each other), by the rule of inherit.py,
worked out from the top down so that only the needed neighbours are factored (PARI).
usage: inherit_big.py check            the rule from the top down against inherit.py's run from the bottom up
       inherit_big.py DIGITS COUNT [K] COUNT primes after random starts of DIGITS digits (seed 20261118 + DIGITS + 1000 K)"""
import sys, random, json, time
from math import gcd
import cypari2

pari = cypari2.Pari()
pari.allocatemem(400 * 10 ** 6)
FAC = {}
MEMO = {}


def items(N):
    if N not in FAC:
        out = []
        for side, m in ((0, N - 1), (1, N + 1)):
            F = pari.factor(m)
            for i in range(len(F[0])):
                p = int(F[0][i])
                out.append((p ** int(F[1][i]), p, side))
        out.sort(reverse=True)
        FAC[N] = out
    return FAC[N]


def atmost(N, g):
    """Is N lit by round g?"""
    if N == 2:
        return True
    if g <= 0:
        return False
    if g == 1:
        L = max((N - 1) & -(N - 1), (N + 1) & -(N + 1))
        return L * L > N
    key = (N, g)
    if key not in MEMO:
        F1 = F2 = 1
        ok = False
        for pe, p, side in items(N):
            if atmost(p, g - 1):
                if side:
                    F2 *= pe
                else:
                    F1 *= pe
                L = F1 * F2 // gcd(F1, F2)
                if L * L > N:
                    ok = True
                    break
        MEMO[key] = ok
    return MEMO[key]


def rounds(N):
    g = 1
    while not atmost(N, g):
        g += 1
    return g


if __name__ == '__main__':
    if sys.argv[1] == 'check':
        from inherit import run
        top = 200000
        d = run(top, 2)
        bad = sum(1 for N in d if N > 2 and rounds(N) != d[N])
        print('rounds from the top down against the run from the bottom up, %d primes up to %d: %d differ' % (len(d) - 1, top, bad))
        for N in (21767, 1688173, 2223990533):
            print('  %d: %d rounds; neighbours %s' % (N, rounds(N), ' and '.join(' x '.join(str(pe) for pe, p, s in sorted(items(N)) if s == side) for side in (0, 1))))
            for pe, p, s in items(N):
                if p > 1000:
                    print('     helper %d needs %d rounds' % (p, rounds(p)))
    else:
        digits, count = int(sys.argv[1]), int(sys.argv[2])
        rng = random.Random(20261118 + digits + (1000 * int(sys.argv[3]) if len(sys.argv) > 3 else 0))
        hist = {}
        t0 = time.time()
        worst = (0, 0)
        for i in range(count):
            N = int(pari.nextprime(rng.randrange(10 ** (digits - 1), 10 ** digits)))
            r = rounds(N)
            hist[r] = hist.get(r, 0) + 1
            if r > worst[0]:
                worst = (r, N)
            if len(FAC) > 200000:
                FAC.clear()
                MEMO.clear()
        n = sum(hist.values())
        mean = sum(k * v for k, v in hist.items()) / n
        print(json.dumps(dict(digits=digits, count=n, mean=round(mean, 4), hist=sorted(hist.items()), worst=worst, seconds=round(time.time() - t0, 1))))
