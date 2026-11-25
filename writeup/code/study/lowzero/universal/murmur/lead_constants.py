"""Constants for the leads (all from Euler products; heuristic predictions, not theorems).

(1) Elliptic curves of prime conductor (ordered by conductor), Sawin-Sutherland local model (their Sec. 1.5:
    a random curve reduces to E_p with probability 1/(p |Aut E_p|)), gives E[a_{p^2}] = l'_{p,2} = -1/p (no weight-4 cusp forms of level 1),
    i.e. E[lambda(p^2)] = -1/p^2. One-point law with lines, per-curve local factor S_p(v)/(1-1/p):
      v = 0 : 1 - lambda(p^2)/p  -> mean 1 + p^{-3};   v = -2: 1/(p-1) (curve-independent).
    Family line amplitudes A(1/q^2) = prod_{p|q} 1/(p-1) * prod_{p not | q} (1 + p^{-3}),  A(1) = zeta(3)/zeta(6).
    Cross-check: the same numbers are the coefficients of 2 pi sqrt(u) J_1(4 pi sqrt(u) m/q) with m/q = 1/q (m = g, q = q g)
    in Sawin-Sutherland's prime-conductor prediction (local factors l'_{p,nu}).
(2) Sym^2 family (degree 3), one-point-law constant c = prod_p G_p(1/p), G_p(X) = 1 - (|lambda(p)|^2 - 1) X + lambda(p^3) X^2 (fact 2057),
    lambda = lambda_{Sym^2 f}: Sato-Tate mean of G_p(1/p) = 1 (harmonic family constant 1);
    Plancherel mean = 1 - p^{-2} + p^{-5} (natural weights, good p).
"""
import mpmath as mp
import sympy as sp

mp.mp.dps = 20
P = list(sp.primerange(2, 3_000_000))


def eprod(f, skip=()):
    s = mp.mpf(0)
    for q in P:
        if q in skip:
            continue
        s += mp.log(f(mp.mpf(q)))
    return mp.e ** s


z = mp.zeta(3) / mp.zeta(6)
print("(1) prime-conductor elliptic curves, family line amplitudes (S-S local model + one-point law):")
print("    A(1)    = zeta(3)/zeta(6) =", mp.nstr(z, 10))
for q in (2, 3, 5, 6, 7):
    f = sp.factorint(q)
    v = z
    for p in f:
        v = v / ((p - 1) * (1 + mp.mpf(p) ** -3))
    print(f"    A(1/{q}^2) =", mp.nstr(v, 10), f"   (harmonic GL(2) family: 1/phi({q}) = {mp.nstr(mp.mpf(1)/sp.totient(q), 6)})")
# cross-check A(1/4) directly from Sawin-Sutherland's double sum: pairs (m, q) = (g, 2g), g odd squarefree
tot = mp.mpf(0)
for g in range(1, 200001, 2):
    fg = sp.factorint(g)
    if any(e > 1 for e in fg.values()):
        continue
    mu = (-1) ** len(fg)
    prodl = mp.mpf(1)
    for p in fg:
        prodl *= mp.mpf(-1) / p
    m, q = g, 2 * g
    # term: mu(gcd(m,q)) / (q m phi(q/gcd)) * prod_{p|m} l'_{p,2} ; amplitude A = term * (q/m)  (A = c/(2 pi sqrt(rho)), c = 2 pi * term)
    term = mu * prodl / (q * m * 1)
    tot += term * mp.mpf(q) / m
print("    A(1/4) from the Sawin-Sutherland double sum (g <= 2e5):", mp.nstr(tot, 10))
c_sym2_nat = eprod(lambda p: 1 - p ** -2 + p ** -5)
print("\n(2) Sym^2 family: harmonic (Sato-Tate) constant 1; natural (Plancherel) constant prod(1 - p^-2 + p^-5) =", mp.nstr(c_sym2_nat, 10),
      "  (compare GL(2) natural 6/pi^2 =", mp.nstr(6 / mp.pi ** 2, 10), ")")
