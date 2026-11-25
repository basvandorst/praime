# Symbolic check: B_p(X) = F_p(X) * prod_{i,j} (1 - a_i b_j X)^{-1} = 1 + O(X^2) for d = 2, 3 (b_j = conj a_j = 1/a_j under Ramanujan)
import sympy as sp
X = sp.symbols('X')
for d in (2, 3):
    a = sp.symbols('a1:%d' % (d + 1)); b = sp.symbols('b1:%d' % (d + 1))
    K = 4
    # lam(p^k) = h_k(a), conj lam(p^k) = h_k(b); mu(p^k) = (-1)^k e_k(a)
    def h(k, v):
        return sp.Poly(sp.series(1 / sp.prod([1 - x * X for x in v]), X, 0, k + 1).removeO(), X).coeff_monomial(X ** k)
    def e(k, v):
        return sp.Poly(sp.expand(sp.prod([1 + x * X for x in v])), X).coeff_monomial(X ** k) if k <= len(v) else 0
    F = sum(h(k, b) * (-1) ** k * e(k, a) * X ** k for k in range(0, d + 1))
    RS = sp.prod([1 / (1 - ai * bj * X) for ai in a for bj in b])
    Bser = sp.series(F * RS, X, 0, 3).removeO()
    c1 = sp.simplify(sp.expand(Bser).coeff(X, 1)); c2 = sp.factor(sp.expand(Bser).coeff(X, 2))
    print('d=%d: coefficient of X in B_p: %s ; coefficient of X^2: %s' % (d, c1, c2))
    # with b = 1/a and det = 1 (d = 2: a2 = 1/a1)
    if d == 2:
        t = sp.symbols('t')
        sub = {a[0]: t, a[1]: 1 / t, b[0]: 1 / t, b[1]: t}
        print('   unitary, det 1: X^2 coefficient = %s ; F_p = %s' % (sp.factor(sp.simplify(c2.subs(sub))), sp.factor(sp.simplify(F.subs(sub)))))
