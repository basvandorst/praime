"""Independent reference for prime6.html (the prime race): counts, time in the lead, changes of the leader and the
first lead of every horse, for the primes (or the Mersenne exponents) up to X.   usage: ref6.py X out.json"""
import sys, json
from math import gcd
X = int(sys.argv[1])
MEXP = [2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,136279841]
lim = min(X, 10 ** 8)
sv = bytearray([1]) * (lim + 1); sv[0] = sv[1] = 0
for i in range(2, int(lim ** 0.5) + 1):
    if sv[i]: sv[i * i::i] = bytearray(len(sv[i * i::i]))
P = [i for i in range(lim + 1) if sv[i]]
out = {}
for src in ('primes', 'mexp'):
    seq = P if src == 'primes' else [p for p in MEXP if p <= X]
    for q in (3, 4, 5, 7, 8, 10, 11, 12, 24, 30):
        cls = [r for r in range(1, q) if gcd(r, q) == 1]; ix = {r: i for i, r in enumerate(cls)}
        cnt = [0] * len(cls); lead = [0] * len(cls); first = [None] * len(cls); k = 0; leader = -1; ch = 0; chg = []
        for p in seq:
            if q % p == 0: continue
            cnt[ix[p % q]] += 1; k += 1
            m = max(cnt); L = cnt.index(m) if cnt.count(m) == 1 else -1
            if L >= 0:
                lead[L] += 1
                if first[L] is None: first[L] = [k, p]
                if L != leader:
                    if leader >= 0: ch += 1
                    chg.append([p, cls[L]]); leader = L
        out['%s|%d' % (src, q)] = dict(cls=cls, cnt=cnt, lead=lead, first=first, k=k, ch=ch, chg=chg[:300])
json.dump(out, open(sys.argv[2], 'w')); print(len(out), 'races')
