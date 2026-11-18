"""A test of the reason given for the rule of machine 23. The reason says: with k sources of helpers (numbers
built from N whose lit part helps), each source has to supply 1/(2k) of the digits, so one round raises the
reachable size to its (2k)-th power and rounds(N) = log_(2k)(ln N) + constant.
k = 1: N - 1.   k = 2: N - 1 and N + 1 (machine 23).   k = 3: also N^2 + 1.   k = 5: also N^2 + N + 1 and N^2 - N + 1.
Only the count of rounds is worked out here (rule: the lcm of the lit parts, squared, exceeds N; helpers are
primes below N). The proof steps for N^2 + 1 and N^2 +- N + 1 exist in the literature (Williams; Lenstra) but
are not built here.
usage: inherit_sources.py K"""
import sys, json, math, random
import cypari2

pari = cypari2.Pari()
pari.allocatemem(400 * 10 ** 6)
K = int(sys.argv[1])
POLY = [lambda N: N - 1, lambda N: N + 1, lambda N: N * N + 1, lambda N: N * N + N + 1, lambda N: N * N - N + 1][:K]
FAC = {}
MEMO = {}


def lcm(a, b):
    return a * b // math.gcd(a, b)


def items(N):
    if N not in FAC:
        out = []
        for side, P in enumerate(POLY):
            F = pari.factor(P(N))
            for i in range(len(F[0])):
                p = int(F[0][i])
                if p < N:
                    out.append((p ** int(F[1][i]), p, side))
        out.sort(reverse=True)
        FAC[N] = out
    return FAC[N]


def atmost(N, g):
    if N == 2:
        return True
    if g <= 0:
        return False
    key = (N, g)
    if key not in MEMO:
        parts = [1] * K
        ok = False
        for pe, p, side in items(N):
            if atmost(p, g - 1):
                parts[side] *= pe
                L = 1
                for x in parts:
                    L = lcm(L, x)
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


rows = []
for a in (10, 12, 14, 16, 18, 20):
    p = int(pari.nextprime(2 ** a))
    tot = n = 0
    while p < 2 ** (a + 1):
        tot += rounds(p)
        n += 1
        p = int(pari.nextprime(p + 1))
    rows.append(('2^%d to 2^%d' % (a, a + 1), math.log(math.log(1.5 * 2 ** a)), n, tot / n))
rng = random.Random(20261118 + K)
for digits, count in ((9, 2000), (12, 2000), (15, 1500), (20, 800), (25, 400), (30, 200)):
    tot = 0
    for _ in range(count):
        tot += rounds(int(pari.nextprime(rng.randrange(10 ** (digits - 1), 10 ** digits))))
        if len(FAC) > 100000:
            FAC.clear()
            MEMO.clear()
    rows.append(('%d digits' % digits, math.log(math.log(5.5 * 10 ** (digits - 1))), count, tot / count))
base = math.log(2 * K)
cs = [m - x / base for name, x, n, m in rows]
n = len(rows)
sx = sum(r[1] for r in rows); sy = sum(r[3] for r in rows); sxx = sum(r[1] ** 2 for r in rows); sxy = sum(r[1] * r[3] for r in rows)
slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
print(json.dumps(dict(k=K, rows=[(a, round(b, 4), c, round(d, 4)) for a, b, c, d in rows], slope=round(slope, 3), predicted_slope=round(1 / base, 3), const_min=round(min(cs), 3), const_max=round(max(cs), 3))))
