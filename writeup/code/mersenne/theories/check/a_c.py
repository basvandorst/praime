# A: trial-factoring survival by p mod 12 (plain pow over all k, no mod-8 or primality filter)
# C: known Mersenne exponents by p mod 12
from sympy import primerange

cnt = {1: 0, 5: 0, 7: 0, 11: 0}
surv = {1: 0, 5: 0, 7: 0, 11: 0}
for p in primerange(100000, 800000):
    c = p % 12
    cnt[c] += 1
    ok = True
    for k in range(1, 101):
        q = 2 * k * p + 1
        if pow(2, p, q) == 1:
            ok = False
            break
    if ok:
        surv[c] += 1
print("A counts", cnt, "total", sum(cnt.values()))
print("A survivors", surv)
print("A shares", {c: round(surv[c] / cnt[c], 4) for c in cnt})

M = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423,
     9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433,
     1257787, 1398269, 2976221, 3021377, 6972593, 13466917, 20996011, 24036583, 25964951, 30402457,
     32582657, 37156667, 42643801, 43112609, 57885161, 74207281, 77232917, 82589933]
print("C len", len(M))
cc = {1: 0, 5: 0, 7: 0, 11: 0}
for p in M:
    if p > 3:
        cc[p % 12] += 1
print("C counts", cc, "sum", sum(cc.values()))
