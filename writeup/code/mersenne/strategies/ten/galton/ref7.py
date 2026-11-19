"""Independent reference for prime7.html (the Galton board): the bins for every rule and source, numbers up to X.
usage: ref7.py X out.json        (the coin rule is random and left out; rows 10)"""
import sys, json, math
X = int(sys.argv[1]); R = 10
MEXP = [2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,136279841]
sv = bytearray([1]) * (X + 1); sv[0] = sv[1] = 0
for i in range(2, int(X ** 0.5) + 1):
    if sv[i]: sv[i * i::i] = bytearray(len(sv[i * i::i]))
OP = [p for p in range(3, 400) if all(p % d for d in range(2, int(p ** 0.5) + 1))][:R]
def bins(rule, nums):
    b = [0] * (R + 1)
    for n in nums:
        if rule == 'bin': c = sum((n >> i) & 1 for i in range(1, R + 1))
        elif rule == 'half': c = sum(1 for p in OP if n % p > p / 2)
        else: c = sum(1 for p in OP if n % p and pow(n % p, (p - 1) // 2, p) == 1)
        b[c] += 1
    return b
out = {}
for rule in ('bin', 'half', 'qr'):
    lo = 2 ** (R + 1) if rule == 'bin' else OP[-1] + 1
    srcs = {'primes': [n for n in range(lo, X + 1) if sv[n]], 'odd': range(lo + (1 - lo % 2), X + 1, 2), 'c6': [n for n in range(lo, X + 1) if n % 2 and n % 3],
            'all': range(lo, X + 1), 'mexp': [p for p in MEXP if lo <= p <= X]}
    for s, nums in srcs.items():
        out['%s|%s' % (rule, s)] = bins(rule, nums)
json.dump(out, open(sys.argv[2], 'w')); print(len(out), 'boards')
