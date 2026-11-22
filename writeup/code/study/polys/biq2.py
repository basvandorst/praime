"""Quartic g(m) = f(s m + t), f(n) = N(n - alpha), alpha in O_K, K = Q(sqrt a, sqrt b). Filter: g has no root mod p for
every prime p <= 23 (necessary for a run > 46) ; then longest run of distinct prime |g|. usage: biq2.py NFIELDS K"""
import sys, json, itertools, time, math
import gmpy2
NF = int(sys.argv[1]); K = int(sys.argv[2])
fields = json.load(open('biq_fields.json'))[:NF]
SMALLP = [2, 3, 5, 7, 11, 13, 17, 19, 23]
best = []; t0 = time.time(); npass = 0
for least, a, b in fields:
    sab = math.sqrt(a * b)
    for c1, c2, c3 in itertools.product(range(-K, K + 1), repeat=3):
        A0, A1, A2, A3 = 2 * c1 + 2 * c2 + c3, 2 * c1 + c3, 2 * c2 + c3, c3
        if A1 == 0 or A2 == 0 or (c1, c2, c3) < (-c1, -c2, -c3): continue
        P2, Q2, R2, PQR = A1 * A1 * a, A2 * A2 * b, A3 * A3 * a * b, A1 * A2 * A3 * a * b
        c_2 = -2 * (P2 + Q2 + R2); c_1 = -8 * PQR; c_0 = P2 * P2 + Q2 * Q2 + R2 * R2 - 2 * (P2 * Q2 + P2 * R2 + Q2 * R2)
        f = lambda n: (lambda yy: (yy ** 4 + c_2 * yy * yy + c_1 * yy + c_0) // 256)(4 * n - A0)
        roots = [(A0 + s1 * A1 * math.sqrt(a) + s2 * A2 * math.sqrt(b) + s1 * s2 * A3 * sab) / 4 for s1 in (1, -1) for s2 in (1, -1)]
        lo, hi = min(roots), max(roots)
        if hi - lo > 300: continue
        for s in (1, 2, 3):
            for t in range(s):
                if all(f(s * m + t) % p for p in SMALLP for m in range(p)) is False: continue
                mlo = math.floor((lo - 60 * s - t) / s); mhi = math.ceil((hi + 60 * s - t) / s)
                npass += 1
                seen = {}; start = mlo; bl = 0; bst = None
                for m in range(mlo, mhi + 1):
                    v = abs(f(s * m + t))
                    if v > 1 and gmpy2.is_prime(v):
                        if v in seen and seen[v] >= start: start = seen[v] + 1
                        seen[v] = m
                        if m - start + 1 > bl: bl, bst = m - start + 1, start
                    else:
                        start = m + 1
                if bl >= 15: best.append((bl, a, b, c1, c2, c3, s, t, bst))
    best.sort(reverse=True)
print('candidates passing the mod-p filter', npass, 'time', round(time.time() - t0), 's')
print(best[:12])
json.dump(best, open('biq2_best.json', 'w'))
