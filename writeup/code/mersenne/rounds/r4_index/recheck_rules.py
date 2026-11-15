"""Recompute the best real step rules in plain Python integers (no numpy, no grammar.py), from the derived
sequences rebuilt by check.py's method (sympy), and count their exact hits and out-of-sample hits."""
import sympy
from fractions import Fraction
import math
REAL = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423,
        9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433,
        1257787, 1398269, 2976221, 3021377, 6972593, 13466917, 20996011, 24036583, 25964951, 30402457,
        32582657, 37156667, 42643801, 43112609, 57885161, 74207281, 77232917, 82589933]
pi = [int(sympy.primepi(a)) for a in REAL]
pblen = [int(sympy.primepi(a)) - int(sympy.primepi(2 ** (a.bit_length() - 1) - 1)) for a in REAL]


def surv_index():  # survivor list written by check.py (plain pow loop, all j <= 1000)
    surv = [int(x) for x in open('check_surv.txt').read().split()]
    d = {v: i + 1 for i, v in enumerate(surv)}
    return [d[a] for a in REAL if a < 10 ** 6]


tf = surv_index()


def binpal(x):
    L = x.bit_length(); return (x << L) | int(bin(x)[2:][::-1], 2)


def rotl1(x):
    L = x.bit_length(); return ((x << 1) & ((1 << L) - 1)) | (x >> (L - 1))


def rnd(fr):  # round half to even like numpy rint
    return round(fr)


rules = {
    'pi: round(sqrt(binpal(x)))': (pi, lambda x: rnd(math.sqrt(binpal(x)))),
    'pblen: floor(e*rotl1(x))': (pblen, lambda x: math.floor(rotl1(x) * math.e)),
    'tf: ceil(x*3/e)': (tf, lambda x: math.ceil(x * 3 / math.e)),
    'tf: round((3x+2)/e)': (tf, lambda x: rnd((3 * x + 2) / math.e)),
}
for nm, (b, f) in rules.items():
    hits = [i for i in range(len(b) - 1) if f(b[i]) == b[i + 1]]
    print(f'{nm:30s} {len(hits)} exact transitions {[(b[i], b[i + 1]) for i in hits]}; '
          f'among targets b_26..: {sum(1 for i in hits if i + 1 >= 25)}')
