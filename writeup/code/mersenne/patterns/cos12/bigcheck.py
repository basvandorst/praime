"""Independent check of every number in the big-picture section (bigdata.py -> bigdata.txt), in plain Python:
a bytearray sieve, one pass over the primes, gcd for the columns, mpmath for li, sympy for spot checks.
usage: bigcheck.py"""
import math, bisect
from math import gcd
import mpmath, sympy

M = 10 ** 7
sv = bytearray([1]) * (M + 3); sv[0] = sv[1] = 0
for i in range(2, int((M + 2) ** .5) + 1):
    if sv[i]: sv[i * i::i] = bytearray(len(range(i * i, M + 3, i)))
ok = []
def check(name, got, want):
    ok.append(got == want); print(('ok  ' if got == want else 'FAIL'), name, got, '' if got == want else f'(want {want})')

UP, DOWN = {1, 11}, {5, 7}
# the pictures, width 240
w = 240
lit = [a for a in range(1, w + 1) if gcd(a, w) == 1]
check('columns of 240 that can hold primes', len(lit), 64)
lane = lambda a: 'U' if a % 12 in UP else 'D'
period = min(k for k in range(1, w + 1) if all(gcd((a + k - 1) % w + 1, w) == 1 and lane(a + k) == lane(a) for a in lit))
check('colour pattern of the stripes repeats every', period, 60)
tracks = [a for a in lit if a % 12 in (5, 11) and gcd(a + 2, w) == 1]
check('twin tracks', len(tracks), 24)
check('twin tracks up / down', (sum(a % 12 == 11 for a in tracks), sum(a % 12 == 5 for a in tracks)), (12, 12))
check('tracks that wrap to the next row', [a for a in tracks if a + 2 > w], [239])
percol = {}
for p in range(5, 100001):
    if sv[p]: a = (p - 1) % w + 1; percol[a] = percol.get(a, 0) + 1
cnt = [percol.get(a, 0) for a in lit]
check('primes per stripe, 5..100,000: min, max', (min(cnt), max(cnt)), (137, 164))
check('mean per stripe x 10', round(10 * sum(cnt) / len(cnt)), 1498)
check('primes 5..100,000 off the stripes', sorted(set(percol) - set(lit)), [5])
tw = [p for p in range(5, 99999) if sv[p] and sv[p + 2]]
off = [(p, p + 2) for p in tw if ((p - 1) % w + 1) not in tracks]
check('twin pairs from (5, 7) to 100,000', len(tw), 1223)
check('twin pairs off the tracks', off, [(5, 7)])
# widths 100 and 277
check('columns of 100 that can hold primes', sum(gcd(a, 100) == 1 for a in range(1, 101)), 40)
check('column 1 of width 100: 1, 101, 201 mod 12', [1 % 12, 101 % 12, 201 % 12], [1, 5, 9])
check('277 is prime and = 1 mod 12', (sympy.isprime(277), 277 % 12), (True, 1))
# the race to 10^7
D = 0; nup = ndn = 0; behind = 0; runs = []; cur = None; xs = []; Ds = []
for p in range(5, M + 1):
    if not sv[p]: continue
    if p % 12 in UP: nup += 1; D -= 1
    else: ndn += 1; D += 1
    xs.append(p); Ds.append(D)
    if D < 0:
        behind += 1
        if cur is None: cur = [p, p, 0, 0]
        cur[1] = p; cur[2] += 1; cur[3] = min(cur[3], D)
    elif cur is not None:
        runs.append(cur); cur = None
if cur: runs.append(cur)
check('primes 5..10^7 up, down', (nup, ndn), (332225, 332352))
check('primes after which the up lane leads', behind, 1476)
check('stretches where the up lane leads', len(runs), 102)
check('first stretch starts at', runs[0][0], 61463)
lg = max(runs, key=lambda r: r[2])
check('longest stretch', (lg[0], lg[1], lg[2]), (1172807, 1179403, 481))
check('last stretch ends at', runs[-1][1], 2233787)
for p in (61463, 1172807, 1179403, 2233787): check(f'{p:,} is prime (sympy)', sympy.isprime(p), True)
half_li = lambda x: 0.5 * float(mpmath.li(math.sqrt(x)) - mpmath.li(2))
g = [math.exp(math.log(1e3) + (math.log(M) - math.log(1e3)) * k / 50000) for k in range(50001)]
r = [Ds[bisect.bisect_right(xs, x) - 1] / half_li(x) for x in g]
avg = sum(r) / len(r)
print(f'     D(x) / (Li(sqrt x)/2) averaged over ln x, 10^3..10^7 (50,001 points): {avg:.4f}')
check('that average, to 2 places', round(avg, 2), 0.99)
# the exponents
EXP = [5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,136279841]
check('exponents 5..136,279,841', len(EXP), 50)
check('exponents up / down', (sum(p % 12 in UP for p in EXP), sum(p % 12 in DOWN for p in EXP)), (16, 34))
check('up to no. 51', (sum(p % 12 in UP for p in EXP[:-1]), len(EXP) - 1), (16, 49))
print('all ok' if all(ok) else f'{ok.count(False)} FAILED')
