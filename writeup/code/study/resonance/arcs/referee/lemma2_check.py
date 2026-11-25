# Referee's own brute-force check of Lemma 2 (arcs_theorem.md) and of the ray count / Corollary 1 (a) algebra.
from math import gcd, isqrt, cos, atan2, pi
import cmath
from sympy import isprime, factorint, totient, mobius

def prim_w(h):
    out = []
    for x in range(-h, h+1):
        y2 = h*h - x*x
        y = isqrt(y2)
        if y*y == y2:
            for yy in ({y, -y}):
                if gcd(abs(x), abs(yy)) == 1:
                    out.append((x, yy))
    return sorted(set(out))

def gauss_primes(N):  # elements a+bi with a^2+b^2 <= N, all associates
    P = []
    for a in range(-isqrt(N), isqrt(N)+1):
        for b in range(-isqrt(N), isqrt(N)+1):
            n = a*a + b*b
            if n < 2 or n > N: continue
            if a != 0 and b != 0:
                if isprime(n): P.append((a, b))
            else:
                q = abs(a) + abs(b)
                if isprime(q) and q % 4 == 3: P.append((a, b))
    return P

ok = True
GP = gauss_primes(20000)
for h in [5, 13, 25, 65, 85]:
    W = prim_w(h)
    om = len(factorint(h))
    cnt_ok = len(W) == 4 * 2**om
    # units of Z[i]/h: zeta = a+bi mod h with gcd(N zeta, h) = 1
    units = [(a, b) for a in range(h) for b in range(h) if gcd(a*a + b*b, h) == 1]
    Phi = len(units)
    res = []
    for (x, y) in W:
        fib = {}
        tot = 0
        for (a, b) in units:
            r = (a*x + b*y) % h          # Re((a+bi)(x-yi)) = ax + by
            fib[r] = fib.get(r, 0) + 1
            tot += cmath.exp(-2j*pi*r/h)
        image_ok = set(fib) == {r for r in range(h) if gcd(r, h) == 1}
        fib_ok = len(set(fib.values())) == 1 and list(fib.values())[0] == int(totient(h))
        mean = tot/Phi
        mean_ok = abs(mean - mobius(h)/totient(h)) < 1e-12
        # (iv): Gaussian primes z not dividing h: Re(z conj w) != 0 mod p for p | h
        bad = 0
        for (a, b) in GP:
            n = a*a + b*b
            if h % n == 0 or (b == 0 and h % abs(a) == 0) or (a == 0 and h % abs(b) == 0):
                continue
            r = a*x + b*y
            if any(r % p == 0 for p in factorint(h)): bad += 1
        res.append(image_ok and fib_ok and mean_ok and bad == 0)
    print(f"h={h}: {len(W)} primitive w (4*2^omega={4*2**om}); Phi(h)={Phi}=phi(h)^2={int(totient(h))**2}; "
          f"all w: image=(Z/h)*, fibres=phi(h), mean=mu/phi={float(mobius(h)/totient(h)):+.6f}, no prime with Re(z w)=0 mod p|h "
          f"(primes N<=2e4): {all(res)}")
    ok &= cnt_ok and all(res) and Phi == int(totient(h))**2
    if h == 65: print("   h=65 rays:", W)
# Corollary 1 (a): sum over rays of hypotenuse h of xi^m = 4 prod_{p|h} 2cos(2 m theta_p), m = 4, squarefree h
m = 4
def theta(p):
    for a in range(1, isqrt(p)+1):
        b = isqrt(p - a*a)
        if a*a + b*b == p: return atan2(b, a)
for h in [5, 13, 65, 85, 1105]:
    W = prim_w(h)
    s = sum(complex(x/h, y/h)**m for (x, y) in W)
    prod = 4
    for p in factorint(h): prod *= 2*cos(2*m*theta(p))
    print(f"h={h}: sum_rays xi^m = {s.real:+.12f}{s.imag:+.1e}i, 4*prod 2cos(2m theta_p) = {prod:+.12f}")
    ok &= abs(s - prod) < 1e-9
print("ALL OK" if ok else "SOME FAILED")
