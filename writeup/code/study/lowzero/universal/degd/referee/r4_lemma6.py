"""Referee check 4: Lemma 6 and Route A local algebra, symbolic (sympy) and numeric.
F_p(X) = sum_k conj(lambda(p^k)) mu_L(p^k) X^k, lambda(p^k) = h_k(alpha), mu_L(p^k) = (-1)^k e_k(alpha).
(i) F_p(1) = 0 under unitarity: conj multiset {conj alpha_j} = {1/alpha_j} (no GRC needed?);
(iii) d = 2, 3 closed forms of G_p = F_p/(1-X); general top coefficient;
Route A, q = p: resonant weight via the Hecke recurrence lambda(p^j) = sum_i (-1)^{i+1} e_i lambda(p^{j-i})."""
import itertools, cmath, math, random
import numpy as np
import sympy as sp

X = sp.Symbol('X')
def hk(vars_, k):
    return sp.Add(*[sp.Mul(*c) for c in itertools.combinations_with_replacement(vars_, k)]) if k > 0 else sp.Integer(1)
def ek(vars_, k):
    return sp.Add(*[sp.Mul(*c) for c in itertools.combinations(vars_, k)]) if k > 0 else sp.Integer(1)

for d in (1, 2, 3, 4):
    a = sp.symbols(f'a1:{d+1}', nonzero=True)
    ainv = [1 / x for x in a]
    # unitarity: conj(h_k(alpha)) = h_k(conj alpha) = h_k(1/alpha)
    F = sum(hk(ainv, k) * (-1) ** k * ek(a, k) * X ** k for k in range(d + 1))
    F1 = sp.simplify(sp.together(F.subs(X, 1)))
    G, rem = sp.div(sp.together(F * sp.Mul(*a) ** d).as_numer_denom()[0], (1 - X), X)
    print(f'd={d}: F_p(1) = {F1};  remainder of F/(1-X): {sp.simplify(rem)}')
    if d in (2, 3):
        Gp = sp.expand(sp.cancel(F / (1 - X)))
        e1, ed = ek(a, 1), ek(a, d)
        lam1c = hk(ainv, 1)             # conj lambda(p)
        absl2 = sp.expand(e1 * lam1c)    # |lambda(p)|^2
        if d == 2:
            claim = 1 - (absl2 - 1) * X
        else:
            claim = 1 - (absl2 - 1) * X + ed * hk(ainv, 3) * X ** 2
        print(f'   d={d}: G_p - claim = {sp.simplify(sp.together(Gp - claim))}')
    # top coefficient of G_p: (-1)^{d+1} omega conj lambda(p^d)
    Gp = sp.cancel(F / (1 - X))
    top = sp.Poly(sp.together(Gp).as_numer_denom()[0], X)
    den = sp.together(Gp).as_numer_denom()[1]
    deg = top.degree()
    lead = sp.simplify(top.LC() / den - (-1) ** (d + 1) * ek(a, d) * hk(ainv, d))
    print(f'   degree of G_p = {deg} (claim d-1 = {d-1}); top coeff - (-1)^(d+1) omega conj lambda(p^d) = {lead}')

# (i) without GRC: non-tempered unitary parameters (complementary-series type pairs alpha, 1/conj(alpha))
rng = np.random.default_rng(7)
worst = 0
worst_noGRC_d2 = None
for d in (2, 3, 4, 5):
    for trial in range(300):
        # build a multiset closed under alpha -> 1/conj(alpha): pairs (r u, u/r) plus unimodular singletons
        al = []
        while len(al) < d:
            if d - len(al) >= 2 and rng.random() < 0.6:
                rr = math.exp(rng.uniform(-0.4, 0.4)); u = cmath.exp(2j * math.pi * rng.random())
                al += [rr * u, u / rr]
            else:
                al.append(cmath.exp(2j * math.pi * rng.random()))
        al = np.array(al)
        # check closure: {conj alpha} == {1/alpha}
        Pp = np.poly(al)  # [1, -e1, e2, ...]  = coefficients of prod(1 - a X) in increasing powers
        # lambda(p^k) = h_k(alpha) via 1/P
        K = d + 2
        lam = [1 + 0j]
        for k in range(1, K + 1):
            lam.append(-sum(Pp[i] * lam[k - i] for i in range(1, min(k, d) + 1)))
        F1 = sum(np.conj(lam[k]) * Pp[k] for k in range(d + 1))
        worst = max(worst, abs(F1))
print(f'(i) without GRC (non-tempered unitary pairs r u, u/r with r up to e^0.4): max |F_p(1)| = {worst:.1e}')
# counter-check: a non-unitary multiset (not closed under alpha -> 1/conj alpha) gives F_p(1) != 0
al = np.array([1.3, 0.5j, -0.9])
Pp = np.poly(al); lam = [1 + 0j]
for k in range(1, 6): lam.append(-sum(Pp[i] * lam[k - i] for i in range(1, min(k, 3) + 1)))
print(f'    non-unitary control: F_p(1) = {sum(np.conj(lam[k]) * Pp[k] for k in range(4)):.4f}')

# Route A, q = p, direct via the Hecke recurrence (d = 2, 3, symbolic):
# T_p = sum_n c_p(n) lam(n) n^{-1/2} I_n = -S(1) + p * sum_k lam(pk) (pk)^{-1/2} I_{pk},
# lam(pk) = sum_{i=1}^d (-1)^{i+1} e_i lam(pk/p^i)  =>  sum_k lam(pk)(pk)^{-1/2} I_{pk} = sum_i (-1)^{i+1} e_i p^{-i/2} S(p^i),
# S(r) := sum_j lam(j) j^{-1/2} I_j(r) = eps sum_m conj lam(m) m^{-1/2} Hhat(ln(r/m))  (G3); resonant m = r: eps conj lam(r) r^{-1/2} int H.
p = sp.Symbol('p', positive=True)
for d in (2, 3):
    a = sp.symbols(f'a1:{d+1}', nonzero=True); ainv = [1 / x for x in a]
    reso = -1 + p * sum((-1) ** (i + 1) * ek(a, i) * p ** sp.Rational(-i, 2) * hk(ainv, i) * p ** sp.Rational(-i, 2) for i in range(1, d + 1))
    F = sum(hk(ainv, k) * (-1) ** k * ek(a, k) * X ** k for k in range(d + 1))
    cLp = F.subs(X, 1 / p) / (1 - 1 / p)
    claim = (p - 1) * (1 - cLp)
    print(f'Route A d={d}, q=p: direct resonant weight - (p-1)(1 - c_Lp) = {sp.simplify(sp.together(reso - claim))}')
    if d == 2:
        # level-1, trivial central character form: weight (p-1)/p * lambda(p^2) (family note 9a)
        b = sp.Symbol('b', nonzero=True)
        sub = {a[0]: b, a[1]: 1 / b}
        lam_p2 = hk([b, 1 / b], 2)
        print(f'   d=2, omega=1: (p-1)(1-c_Lp) - (p-1)/p*lambda(p^2) = {sp.simplify(claim.subs(sub) - (p - 1) / p * lam_p2)}')
