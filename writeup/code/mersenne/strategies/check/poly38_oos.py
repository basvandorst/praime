"""Out-of-sample check of the near-miss ray 19m^2 + 15m + 1 (polygon spiral P = 38, corner 34): the 20 exponents above 500,000,
and the heuristic expected count of exponents on the ray between 500,000 and 136,279,841 (and up to 10^9)."""
import math, gmpy2
BIG = [756839,859433,1257787,1398269,2976221,3021377,6972593,13466917,20996011,24036583,25964951,30402457,32582657,37156667,42643801,43112609,57885161,74207281,77232917,82589933,136279841]
f = lambda m: 19 * m * m + 15 * m + 1
def on_ray(p):
    d = 225 + 76 * (p - 1); r = math.isqrt(d)
    return r * r == d and (r - 15) % 38 == 0
print('exponents above 500,000 on the ray:', [p for p in BIG if on_ray(p)])
print('check of the three hits:', [(m, f(m)) for m in (2, 32, 48)])
w = lambda p: min(1.0, math.exp(0.5772156649) * math.log2((6 if p % 4 == 1 else 2) * p) / p)
for hi in (136279841, 10**9):
    e = sum(w(f(m)) for m in range(1, 10**5) if 500000 < f(m) <= hi and gmpy2.is_prime(f(m)))
    print(f'expected exponents on the ray in (500,000, {hi:,}]: {e:.4f}')
