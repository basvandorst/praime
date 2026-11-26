# Exact local algebra (sympy) for the averaged one-point law of S4 quartic fields (degree-3 Artin L-function zeta_K/zeta).
# Per field, line r with v = v_l(r): F(v) = (l/(l-1)) sum_{h>=max(0,-v)} mu_L(l^h) lam(l^(v+h)) l^(-h)   (method of the cubic case)
# Family average at l with Bhargava masses (sum over etale quartic algebras of 1/|Aut| 1/Disc), total 1 + 1/l + 2/l^2 + 1/l^3.
import sympy as sp
l, X = sp.symbols('l X', positive=True)
J = 24
# type: (mass, 1/L_l(X) as polynomial in X)   1/L = prod_P (1 - X^f_P) / (1 - X)
T = {'1111': (sp.Rational(1, 24), (1 - X)**3), '211': (sp.Rational(1, 4), (1 - X)**2 * (1 + X)), '22': (sp.Rational(1, 8), (1 + X)**2 * (1 - X)),
     '31': (sp.Rational(1, 3), 1 - X**3), '4': (sp.Rational(1, 4), 1 + X + X**2 + X**3),
     '1^2 11': (1 / (2 * l), (1 - X)**2), '1^2 2': (1 / (2 * l), 1 - X**2),
     '1^2 1^2': (1 / (2 * l**2), 1 - X), '2^2': (1 / (2 * l**2), 1 + X), '1^3 1': (1 / l**2, 1 - X), '1^4': (1 / l**3, sp.Integer(1))}
Z = 1 + 1 / l + 2 / l**2 + 1 / l**3
assert sp.simplify(sum(m for m, _ in T.values()) - Z) == 0
def coeffs(name):
    p = sp.Poly(sp.expand(T[name][1]), X); mu = [p.coeff_monomial(X**k) for k in range(4)]
    lam = sp.Poly(sp.series(1 / T[name][1], X, 0, J + 4).removeO(), X); lam = [lam.coeff_monomial(X**k) for k in range(J + 4)]
    return mu, lam
def F(name, v):
    mu, lam = coeffs(name)
    return sp.cancel(l / (l - 1) * sum(mu[h] * lam[v + h] * l**(-h) for h in range(max(0, -v), 4) if v + h >= 0))
if __name__ == '__main__':
    # a_K(l) for each type and Chebotarev / Bhargava check of E[a] at fixed l
    for n in T: print('%-8s mass %-10s lam(l), lam(l^2), lam(l^3): %s   F(v=-3..3): %s' % (n, T[n][0], coeffs(n)[1][1:4], [F(n, v) for v in range(-3, 4)]))
    b = {v: sp.factor(sp.cancel(sum(T[n][0] * F(n, v) for n in T) / Z)) for v in range(-5, J)}
    print('family average b_l(v):'); [print('  v=%2d:' % v, b[v]) for v in range(-5, 13)]
    E = sum(b[v] * X**v for v in range(-3, J))
    # closed form: guess E_l(X) * (1-X^2)(1-X^3)(1-X^4) is a Laurent polynomial
    for name, P in (('(1-X^2)(1-X^3)(1-X^4)', (1 - X**2) * (1 - X**3) * (1 - X**4)), ('(1-X^2)(1-X^3)', (1 - X**2) * (1 - X**3)), ('(1-X)(1-X^2)(1-X^3)', (1-X)*(1-X**2)*(1-X**3))):
        H = sp.expand(sp.series(sp.expand(E * P * X**3), X, 0, J - 6).removeO())
        print(name, '* E * X^3 =', sp.factor(H) if sp.Poly(H, X).degree() < 15 else 'degree %d (no closure)' % sp.Poly(H, X).degree())
