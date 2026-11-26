"""Compare densities: C1 prime-level M_k (fact 2060), Tomczak's M_{inf,k} (arXiv 2606.08353),
BKLMY/Tomczak depth density M_{l,k} (arXiv 2603.25564 Thm 1.1 = 2606.08353 def.), and the C2 recipe
(fact 2061) applied to the prime set {l}.  Pointwise densities in y = p/N."""
import math
from sympy import primerange, factorint

# Artin's constant by direct product to 2e6 plus tail sum_{p>P} 1/p^2 approx
PMAX = 2_000_000
logA = 0.0
for p in primerange(2, PMAX):
    logA += math.log1p(-1.0 / (p * (p - 1)))
# tail: sum_{p>P} 1/p^2 ~ 1/(P log P)
logA += -1.0 / (PMAX * math.log(PMAX))
ARTIN = math.exp(logA)
Z2 = math.pi ** 2 / 6


def U(n, x):
    # Chebyshev second kind
    if n == 0:
        return 1.0
    u0, u1 = 1.0, 2 * x
    for _ in range(n - 1):
        u0, u1 = u1, 2 * x * u1 - u0
    return u1


def A(x):
    """prod_{p not| x} (p^2-p-1)/(p(p-1)); A(0)=1."""
    if x == 0:
        return 1.0
    val = ARTIN
    for p in factorint(abs(x)):
        val /= (1 - 1.0 / (p * (p - 1)))
    return val


def m_loc(l, r):
    """local mean m_l(r) of fact 2060"""
    if l == 2:
        return 4 / 3 if r % 2 == 0 else 2 / 3
    return l * l / (l * l - 1) if r % l == 0 else 1 - 1 / ((l - 1) * (l * l - 1))


def nu_prime(r):
    """nu'(r) = prod_l m_l(r): product over primes; equals ARTIN-type closed form"""
    # closed form via generic product * corrections at primes dividing r and at 2
    # generic odd factor g_l = 1 - 1/((l-1)(l^2-1)); product over odd l of g_l
    return NU1 / (2 / 3) * m_loc(2, r) * math.prod(
        (m_loc(p, r) / m_loc(p, 1)) for p in factorint(r) if p != 2)


# nu'(1) by direct product
lg = math.log(2 / 3)
for p in primerange(3, PMAX):
    lg += math.log1p(-1.0 / ((p - 1) * (p * p - 1)))
lg += -1.0 / (PMAX ** 2 * 2 * math.log(PMAX))
NU1 = math.exp(lg)


def M_C1(y, k):
    s = 0.0
    r = 1
    while r * r < 4 * y:
        s += nu_prime(r) * U(k - 2, r / (2 * math.sqrt(y))) * math.sqrt(4 * y - r * r)
        r += 1
    val = (2 * math.pi * math.sqrt(y) + (-1) ** (k // 2 - 1) * (12 / math.pi) * s) / (k - 1)
    if k == 2:
        val -= 12 * y
    return val


def M_tom_inf(v, k):
    s = 0.0
    smax = int(2 * math.sqrt(v)) + 1
    for t in range(-smax, smax + 1):
        if t * t < 4 * v:
            s += A(t) * math.sqrt(v - t * t / 4) * U(k - 2, t / (2 * math.sqrt(v)))
    val = (-1) ** (k // 2 + 1) * 2 * math.pi / (k - 1) * s
    if k == 2:
        val -= 12 * v
    return val


def M_depth(v, l, k):
    s = 0.0
    tmax = int(2 * l * math.sqrt(v)) + 1
    for t in range(-tmax, tmax + 1):
        if t * t < 4 * l * l * v:
            w = (1.0 if t % l == 0 else 0.0) - 1.0 / l
            s += w * A(l * t) * math.sqrt(v - t * t / (4 * l * l)) * U(k - 2, t / (2 * l * math.sqrt(v)))
    return (-1) ** (k // 2 + 1) * 2 * math.pi / ((k - 1) * (1 - 1 / l)) * s


def M_C2(y, S, k):
    """fact 2061 recipe for prime set S"""
    pi_fac = math.prod(1 + 1 / q for q in S)
    n_phi = math.prod(q / (q - 1) for q in S)
    s = 0.0
    r = 1
    while r * r < 4 * y:
        nuN = nu_prime(r) / math.prod(m_loc(q, r) for q in S)
        s += nuN * U(k - 2, r / (2 * math.sqrt(y))) * math.sqrt(4 * y - r * r)
        r += 1
    val = (2 * math.pi * pi_fac * math.sqrt(y) + (-1) ** (k // 2 - 1) * (12 / math.pi) * n_phi * s) / (k - 1)
    if k == 2:
        val -= 12 * n_phi * y
    return val


if __name__ == "__main__":
    print(f"Artin {ARTIN:.10f}  nu'(1) {NU1:.7f}  zeta(2)*Artin {Z2*ARTIN:.7f}")
    # identity nu'(r) = zeta(2) A(r)
    worst = max(abs(nu_prime(r) - Z2 * A(r)) for r in range(1, 200))
    print(f"max_r<200 |nu'(r) - zeta(2)A_r| = {worst:.2e}")
    ys = [0.2, 0.5, 1.0, 2.0, 3.3]
    for k in (2, 4, 8):
        print(f"\nk={k}")
        print("   y    C1 M_k    Tomczak M_inf   |diff|")
        for y in ys:
            a, b = M_C1(y, k), M_tom_inf(y, k)
            print(f"{y:5.2f} {a:10.6f} {b:12.6f}   {abs(a-b):.1e}")
    k = 2
    print("\nk=2: depth density M_{l,2} (2603.25564) vs C2 recipe for prime set {l} vs C1 (l->inf)")
    for l in (3, 5, 7, 101):
        print(f" l={l}")
        for y in (0.2, 0.5, 1.0):
            d, c2, c1 = M_depth(y, l, k), M_C2(y, [l], k), M_C1(y, k)
            print(f"   y={y:4.2f}  depth {d:9.5f}   C2({{l}}) {c2:9.5f}   C1 {c1:9.5f}   depth-C1 {d-c1:+.4f}  C2-C1 {c2-c1:+.4f}")
    # Zubrilina value at 0.2 for reference via C3 constants is not recomputed here
