"""Independent recheck of the key numbers: primes from sympy, angles from scratch, exact direct
sums, a rejection-sampling null (instead of the DP sampler), sympy Legendre symbols."""
import numpy as np, math, random
from sympy import primerange, legendre_symbol, totient
N = 500000
pr = np.array(list(primerange(2, N+1)), dtype=np.int64)
print('pi(N) =', len(pr))
# wheel: spokes are n mod 708 (as (n-1)%708 + 1)
sp = pr % 708
cnt = np.bincount(sp, minlength=708)
cop = [a for a in range(708) if math.gcd(a, 708) == 1]
print('wheel: spokes with >=1 prime', (cnt > 0).sum(), '; coprime spokes', len(cop), '; non-coprime spokes with primes', [a for a in range(708) if cnt[a] and math.gcd(a, 708) > 1])
cc = cnt[cop]; print('  primes per coprime spoke: min', cc.min(), 'max', cc.max(), 'mean', cc.mean().round(1), 'chi2', (((cc-cc.mean())**2)/cc.mean()).sum().round(1), 'dof', len(cc)-1)
# polar: arms mod 44 / 710 beyond 710
big = pr[pr > 710]
print('polar: residues mod 44 with primes', len(set((big % 44).tolist())), '(phi=%d)' % totient(44), '; mod 710:', len(set((big % 710).tolist())), '(phi=%d)' % totient(710))
# direct exact harmonic sums in blocks (no histogram): polar k=355, 71 in 25000-blocks n>=50000; sunflower k=1292 in 3000 blocks
def blockamp(theta_fn, k, n0, n1, D):
    num = 0.0; den = 0.0
    for a in range(n0, n1, D):
        b = min(a+D, n1)
        if b-a < D//2: continue
        p = pr[(pr >= a) & (pr < b)].astype(float)
        S = np.exp(1j*k*theta_fn(p)).sum()
        num += abs(S)**2 - len(p); den += len(p)**2
    return math.sqrt(max(num, 0)/den)
print('polar k=355 (pred 0.977): %.3f   k=71 (pred 0.250): %.3f' % (blockamp(lambda n: n, 355, 50000, N+1, 25000), blockamp(lambda n: n, 71, 50000, N+1, 25000)))
ga = math.pi*(3-math.sqrt(5))
print('sunflower k=1292 D=3000 (pred 0.893): %.3f   k=963 (pred 0.470): %.3f' % (blockamp(lambda n: n*ga, 1292, 100000, N+1, 3000), blockamp(lambda n: n*ga, 963, 100000, N+1, 3000)))
print('wheel k=354 (pred 1): %.3f  k=236 (pred .5): %.3f  k=101 (pred 0): %.3f' % tuple(blockamp(lambda n: -2*math.pi*(n-1)/708, k, 1, N+1, 7080) for k in (354, 236, 101)))
# Sacks: two curves, Bateman-Horn with sympy Legendre symbols, product over p < 1e5
ps = list(primerange(3, 100000)); isp = set(pr.tolist())
for label, f, disc, two in (('x^2+x+41', lambda x: x*x+x+41, 1-4*41, 2.0), ('x^2+1', lambda x: x*x+1, -4, 1.0), ('x^2+58', lambda x: x*x+58, -4*58, 1.0)):
    Cf = two
    for p in ps:
        om = 1 if disc % p == 0 else 1 + legendre_symbol(disc % p, p)
        Cf *= (1 - om/p)/(1 - 1/p)
    vals = [f(x) for x in range(10, 800) if 2 <= f(x) <= N]
    obs = sum(v in isp for v in vals); exp_ = Cf*sum(1/math.log(v) for v in vals)
    print(f'sacks {label}: C = {Cf:.3f}, primes {obs}, BH {exp_:.1f}')
# exponents: rejection-sampling null, 1000 accepted sets
E = [2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091]
g = math.exp(0.5772156649015329)
w = np.array([min(1.0, g*math.log2((6 if p % 4 == 1 else 2)*p)/p) for p in pr.tolist()])
rng = np.random.default_rng(7); acc = []
while len(acc) < 1000:
    U = rng.random((100, len(pr))) < w
    for row in U[U.sum(1) == 31]:
        acc.append(pr[row])
acc = np.array(acc[:1000])
def st(S):
    S = np.asarray(S, float)
    ray = np.abs(np.exp(1j*S).sum(-1))**2/S.shape[-1]          # polar: angle n rad
    h = 2.0**np.floor(np.log2(S)); t = (S-h+0.5)/h
    return ray, (t < 0.25).sum(-1), (np.mod(np.log2(S), 1) < 0.25).sum(-1)
o = st(np.array(E)[None, :]); nl = st(acc)
for nm, a, b in zip(('polar Rayleigh k=1', 'binring count t<0.25', 'doubling count frac(log2)<0.25'), o, nl):
    print(f'{nm}: obs {a[0]:.2f}, null mean {b.mean():.2f}, p = {(1+np.sum(b >= a[0]))/(1+len(b)):.3f}')
