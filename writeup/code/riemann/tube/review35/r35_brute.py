import math
from sympy import primerange, isprime
d1 = 0.8276794767155048879910469896773736467539; B0 = -0.02202465783872691842463476000992357551451; C = 0.15655695399428649
chi = [0, 1, -1, -1, 1]
def Te(y): return math.fsum(y ** (-n) / (n * (n + 1)) for n in range(2, 400, 2))
def G(y):
    return math.fsum([math.log(p) * chi[p % 5] ** k * (y / p ** k - 1) for p in primerange(2, int(y) + 1) for k in range(1, 64) if p ** k <= y] + [y * d1])
def ratio(y, sgn=1, const=1 + B0, te=1):
    return abs(G(y) - sgn * math.log(y) - const + te * Te(y)) / math.sqrt(y)
for y in [1, 1.5, 2, 2.5, 3, 4.5, 7, 10, 30, 100, 3109, 1350457]:
    print(y, round(ratio(y), 9), ' wrong ln sign:', round(ratio(y, -1), 4), ' const B0 only:', round(ratio(y, 1, B0), 4), ' T_e sign flipped:', round(ratio(y, 1, 1 + B0, -1), 4))
print('primes to 1e9 count check skipped; pi(1e6)=', sum(1 for _ in primerange(2, 10**6)))
