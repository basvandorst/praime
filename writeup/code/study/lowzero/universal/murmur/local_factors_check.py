"""Exact local-factor checks linking the one-point law (with lines) to murmuration densities.

Notation (repo, family_proof.md section 1 and Proposition D'):
  per-form line amplitude, local factor at p:  S_p(v)/(1-1/p),
  S_p(v) = sum_{m-h=v} lambda(p^m) mu_f(p^h) p^{-(m+h)/2},
  good p: mu_f(1)=1, mu_f(p)=-lambda(p), mu_f(p^2)=1, mu_f(p^h)=0 (h>=3), lambda(p^j)=U_j(x/2);
  bad p (p||N): lambda(p^m)=lambda(p)^m, lambda(p)^2=1/p, mu_f(p)=-lambda(p), mu_f(p^h)=0 (h>=2).

Checks:
 (1) Sato-Tate (harmonic) average at good p, v in -6..6: 1 (v=0), 1/(p-1) (v=-2), 0 else
     -> line amplitudes mu^2(m)/phi(m) at r=1/m^2  (Theorem D, l=1; Miller-Montague chi(u)).
 (2) Plancherel (natural weights) average at good p and the bad-prime value, mixed with
     P(p|N) = (p-1)/(p^2+p-1) (squarefree N weighted by phi(N)), equals Zubrilina's local factor
     L_p(v) read off from her Theorem 1.3 Bessel series:
        L_p(0) = (p^4-2p^2-p+1)/(p^4-2p^2+p) * (1+Q(p)/p),  L_p(-2) = (p^4-2p^2-p+1)/(p^4-2p^2+p) * p Q(p),
        L_p(2e) = p^{-2e} L_p(0),  L_p(v)=0 otherwise,  Q(p) = p^2/(p^4-2p^2-p+1).
 (3) Numerical constants: alpha, beta, gamma of Zubrilina; sum Q(d) = beta/alpha; alpha*sum Q(d)/d = gamma/pi;
     prod_p L_p(0) = gamma/(2 pi^2).
 (4) Sawin-Sutherland prime-conductor prediction: coefficient of 2 pi sqrt(u) J_1(4 pi sqrt u) equals
     sum_{q sqfree} mu(q)/q^2 * prod_{p|q} l'_{p,2} with l'_{p,2} = -1/p  ->  zeta(3)/zeta(6);
     and the one-point-law average E[c_E] = prod_p (1 - E[lambda(p^2)]/p) with E[a_{p^2}] = l'_{p,2} = -1/p.
"""
import sympy as sp
import mpmath as mp

p, x = sp.symbols('p x', positive=True)


def U(j, xx):
    # Chebyshev U_j(x/2): lambda(p^j) for a good prime, x = lambda(p)
    if j < 0:
        return sp.Integer(0)
    a, b = sp.Integer(1), xx
    if j == 0:
        return a
    for _ in range(j - 1):
        a, b = b, sp.expand(xx * b - a)
    return b


def S_good(v):
    tot = 0
    mu = {0: sp.Integer(1), 1: -x, 2: sp.Integer(1)}
    for h in range(0, 3):
        m = v + h
        if m < 0:
            continue
        tot += U(m, x) * mu[h] * p ** sp.Rational(-(m + h), 2)
    return sp.expand(tot)


def S_bad(v, lam):
    # lam = lambda(p) = +-p^{-1/2}
    tot = 0
    mu = {0: sp.Integer(1), 1: -lam}
    for h in range(0, 2):
        m = v + h
        if m < 0:
            continue
        tot += lam ** m * mu[h] * p ** sp.Rational(-(m + h), 2)
    return sp.simplify(tot)


def poly_moments(expr, moment):
    """Replace x^j by moment(j) in a polynomial in x."""
    P = sp.Poly(sp.expand(expr), x)
    return sp.simplify(sum(c * moment(mon[0]) for mon, c in zip(P.monoms(), P.coeffs())))


def st_moment(j):
    # Sato-Tate moments of x = 2cos(theta): E[x^j] = Catalan(j/2) for even j
    return sp.catalan(j // 2) if j % 2 == 0 else 0


def pl_moment(j):
    # p-adic Plancherel measure: E[U_n(x/2)] = p^{-n/2} (n even), 0 (n odd); get raw moments x^j
    # x^j = sum_n c_{j,n} U_n(x/2); compute by inverting
    tot = 0
    # express x^j in Chebyshev basis via recursion: x*U_n = U_{n+1}+U_{n-1}
    coeffs = {0: sp.Integer(1)}
    for _ in range(j):
        new = {}
        for n, c in coeffs.items():
            new[n + 1] = new.get(n + 1, 0) + c
            if n >= 1:
                new[n - 1] = new.get(n - 1, 0) + c
        coeffs = new
    for n, c in coeffs.items():
        if n % 2 == 0:
            tot += c * p ** sp.Rational(-n, 2)
    return tot


print("(1) Sato-Tate (harmonic) average of S_p(v)/(1-1/p), good p:")
for v in range(-6, 7):
    val = sp.simplify(poly_moments(S_good(v), st_moment) / (1 - 1 / p))
    print(f"   v={v:+d}: {sp.factor(val)}")

Q = p ** 2 / (p ** 4 - 2 * p ** 2 - p + 1)
alpha_p = (p ** 4 - 2 * p ** 2 - p + 1) / (p ** 4 - 2 * p ** 2 + p)


def L_zub(v):
    if v == 0:
        return sp.simplify(alpha_p * (1 + Q / p))
    if v == -2:
        return sp.simplify(alpha_p * p * Q)
    if v > 0 and v % 2 == 0:
        return sp.simplify(p ** (-v) * alpha_p * (1 + Q / p))
    return sp.Integer(0)


Pbad = (p - 1) / (p ** 2 + p - 1)
Pgood = p ** 2 / (p ** 2 + p - 1)
print("\n(2) natural (Plancherel + bad primes, squarefree levels) vs Zubrilina local factor L_p(v):")
allok = True
for v in range(-6, 7):
    good = poly_moments(S_good(v), pl_moment) / (1 - 1 / p)
    lam = p ** sp.Rational(-1, 2)
    bad = (S_bad(v, lam) + S_bad(v, -lam)) / 2 / (1 - 1 / p)   # sign of lambda(p) equidistributed
    mix = sp.simplify(Pgood * good + Pbad * bad)
    z = L_zub(v)
    ok = sp.simplify(mix - z) == 0
    allok &= ok
    print(f"   v={v:+d}: one-point-law average = {sp.factor(mix)}   Zubrilina = {sp.factor(z)}   equal: {ok}")
print("   ALL EQUAL:", allok)
print("   prime-level natural local factors (p != N): v=0:", sp.factor(poly_moments(S_good(0), pl_moment) / (1 - 1 / p)),
      " v=-2:", sp.factor(poly_moments(S_good(-2), pl_moment) / (1 - 1 / p)),
      " v=2:", sp.factor(poly_moments(S_good(2), pl_moment) / (1 - 1 / p)))
print("   bad prime v=0:", sp.factor((S_bad(0, p ** sp.Rational(-1, 2)) + S_bad(0, -p ** sp.Rational(-1, 2))) / 2 / (1 - 1 / p)))

# (3) constants
mp.mp.dps = 30
primes = list(sp.primerange(2, 2_000_000))


def eprod(f):
    s = mp.mpf(0)
    for q in primes:
        s += mp.log(f(mp.mpf(q)))
    return mp.e ** s


alpha = 2 * mp.pi * eprod(lambda q: (1 - q - 2 * q ** 2 + q ** 4) / (q ** 4 - 2 * q ** 2 + q))
beta = 2 * mp.pi * eprod(lambda q: (-1 + q ** 2 + q ** 3) / (q * (-1 + q + q ** 2)))
gamma = 12 * eprod(lambda q: q * (1 + q) / (-1 + q + q ** 2))
sumQ = eprod(lambda q: 1 + q ** 2 / (q ** 4 - 2 * q ** 2 - q + 1))
sumQd = eprod(lambda q: 1 + q / (q ** 4 - 2 * q ** 2 - q + 1))
L0 = eprod(lambda q: (q - 1) * (q + 1) ** 2 / (q * (q ** 2 + q - 1)))
print("\n(3) constants (Euler products over p < 2e6):")
print("   alpha =", mp.nstr(alpha, 12), " beta =", mp.nstr(beta, 12), " gamma =", mp.nstr(gamma, 12))
print("   sum Q(d) =", mp.nstr(sumQ, 12), "  beta/alpha =", mp.nstr(beta / alpha, 12))
print("   alpha*sum Q(d)/d =", mp.nstr(alpha * sumQd, 12), "  gamma/pi =", mp.nstr(gamma / mp.pi, 12))
print("   prod_p L_p(0) =", mp.nstr(L0, 12), "  gamma/(2 pi^2) =", mp.nstr(gamma / (2 * mp.pi ** 2), 12),
      "  (alpha/2pi) sum Q(d)/d =", mp.nstr(alpha / (2 * mp.pi) * sumQd, 12))
print("   prime level natural: 6/pi^2 =", mp.nstr(6 / mp.pi ** 2, 12))
L2_0 = mp.mpf((2 - 1) * (2 + 1) ** 2) / (2 * (4 + 2 - 1))
L2_m2 = mp.mpf(4) / ((2 - 1) * (4 + 2 - 1))
print("   squarefree-natural line amplitudes: A(1)=%s  A(1/4)=%s  A(4)=%s  A(1/9)=%s" % (
    mp.nstr(L0, 6), mp.nstr(L0 * L2_m2 / L2_0, 6), mp.nstr(L0 / 4, 6),
    mp.nstr(L0 * (mp.mpf(9) / (2 * (9 + 3 - 1))) / (mp.mpf(2 * 16) / (3 * (9 + 3 - 1))), 6)))

# (4) Sawin-Sutherland prime-conductor prediction, resonant coefficient m = q
z3z6 = mp.zeta(3) / mp.zeta(6)
tot = mp.mpf(0)
for q in range(1, 200001):
    if any(e > 1 for e in sp.factorint(q).values()):
        continue
    mu = (-1) ** len(sp.factorint(q))
    prodl = mp.mpf(1)
    for r in sp.factorint(q):
        prodl *= mp.mpf(-1) / r
    tot += mu * prodl / mp.mpf(q) ** 2
print("\n(4) Sawin-Sutherland m=q coefficient (q<=2e5):", mp.nstr(tot, 10), "  zeta(3)/zeta(6) =", mp.nstr(z3z6, 10))
EcE = eprod(lambda q: 1 + 1 / q ** 3)
print("    one-point law E[c_E] = prod(1 - E[lambda(p^2)]/p) with E[lambda(p^2)] = -1/p^2:", mp.nstr(EcE, 10))
