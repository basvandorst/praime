# Exact local algebra (sympy) for the averaged one-point law of S3 cubic fields.
# Per field, line r with v = v_l(r): local factor (LL normalisation, = l^{v/2} x the dictionary normalisation A)
#   F(v) = (l/(l-1)) * sum_{h >= max(0,-v)} mu_L(l^h) * lam(l^(v+h)) * l^(-h)      (HL sieve at l, exact dual by fact 2056(ii))
# Family average at l with Bhargava masses: (111) 1/6, (21) 1/2, (3) 1/3, (1^2 1) 1/l, (1^3) 1/l^2, total 1 + 1/l + 1/l^2.
import sympy as sp, mpmath as mp
l, X, t = sp.symbols('l X t', positive=True)
J = 16
types = {  # name: (mass, coefficients of 1/L_l(X) = mu_L, chi(l))
    '111': (sp.Rational(1, 6), [1, -2, 1], 1), '21': (sp.Rational(1, 2), [1, 0, -1], -1),
    '3': (sp.Rational(1, 3), [1, 1, 1], 1), '1^21': (1 / l, [1, -1], 0), '1^3': (1 / l**2, [1], 0)}
Z = 1 + 1 / l + 1 / l**2

def lam_seq(mu, n):                      # coefficients of 1/(sum mu_k X^k)
    lam = [sp.Integer(1)]
    for j in range(1, n):
        lam.append(-sum(mu[k] * lam[j - k] for k in range(1, len(mu)) if j - k >= 0))
    return lam

def F(name, v):
    mass, mu, chi = types[name]; lam = lam_seq(mu, J + 4)
    tot = sum(mu[h] * lam[v + h] * l**(-h) for h in range(max(0, -v), len(mu)) if v + h >= 0)
    return sp.cancel(l / (l - 1) * tot)

ok = True
for name in ('111', '21', '3'):
    mass, mu, chi = types[name]; lam = lam_seq(mu, J + 4)
    for v in range(-4, 10):
        vor = (lam[v] - chi * lam[v + 2] / l) if v >= 0 else (-chi * lam[1] / l if v == -1 else (sp.Integer(chi) / (l * (l - 1)) if v == -2 else 0))
        ok &= sp.cancel(F(name, v) - vor) == 0
print('good types: HL-sieve factor == Voronoi line weights (lines_lib), v = -4..9:', ok)
for name in ('1^21', '1^3'):
    print('ramified %-5s F(v), v=-2..4:' % name, [F(name, v) for v in range(-2, 5)])
b = {v: sp.factor(sp.cancel(sum(types[n][0] * F(n, v) for n in types) / Z)) for v in range(-4, J)}
print('family average b_l(v):'); [print('  v=%2d:' % v, b[v]) for v in range(-4, 9)]
print('b_l(0) == l^3/(l^3-1):', sp.cancel(b[0] - l**3 / (l**3 - 1)) == 0)
print('beta_v = b(v)/b(0):', {v: sp.factor(sp.cancel(b[v] / b[0])) for v in range(-2, 9)})
P = 1 / ((1 - X**2) * (1 - X**3))
E = (P * (1 - X / l) + 1 / (l * (1 - X)) + 1 / (l * (l - 1)) - 1 / (l * (l - 1) * X)) / Z
num = sp.series(sp.cancel(E * X), X, 0, J - 1).removeO()
print('closed form of E_l(X) = sum_v b_l(v) X^v matches, v = -1..%d:' % (J - 3), all(sp.cancel(sp.expand(num).coeff(X, v + 1) - b[v]) == 0 for v in range(-1, J - 2)))
H = sp.factor(sp.cancel(E * (1 - X**2) * (1 - X**3)))
print('H_l(X) = E_l(X)(1-X^2)(1-X^3) =', H)
print('H_l(l^-1/2) =', sp.factor(sp.cancel(H.subs(X, t).subs(l, t**-2))), '  (t = l^-1/2; claim 1 - t^3)')
print('H_l(l^-1/3) =', sp.factor(sp.cancel(H.subs(X, t).subs(l, t**-3))), '  (t = l^-1/3; claim (1+t^3)(1-t^5)/(1+t^3+t^6))')
# consequences: D(s) = zeta(2s) zeta(3s) H(s);  H(1/2) = 1/zeta(3/2);  H(1/3) = zeta(3)/(zeta(2) zeta(5/3))
mp.mp.dps = 30
A_half = mp.zeta(1.5) / mp.zeta(1.5) / 2
z3 = mp.zeta(3); H13 = z3 / (mp.zeta(2) * mp.zeta(mp.mpf(5) / 3))
for sgn, kap in (('-', (2 * mp.pi)**(-mp.mpf(2) / 3) * mp.gamma(mp.mpf(1) / 3) / mp.gamma(mp.mpf(2) / 3)),
                 ('+', (2 * mp.pi)**(-mp.mpf(2) / 3) * mp.gamma(mp.mpf(1) / 3)**2 * (1 + mp.cos(mp.pi / 3)) / mp.pi)):
    A13 = 2 * mp.pi * mp.zeta(mp.mpf(2) / 3) / 3 * H13 * kap
    KC = 1 / mp.sqrt(3) if sgn == '-' else 1
    cTT = KC * 8 * z3 * mp.zeta(mp.mpf(1) / 3) / (mp.gamma(mp.mpf(2) / 3)**3 * mp.zeta(mp.mpf(5) / 3))
    print('D%s0: y^(-1/3) coefficient of the averaged law %s; Taniguchi-Thorne secondary term (5/6)(K/C)12 zeta(3) kappa_TT: %s' % (sgn, mp.nstr(A13, 15), mp.nstr(cTT, 15)))
print('y^(-1/2) coefficient of the line sum: zeta(3/2) H(1/2) / 2 = 1/2 (prime squares, E[lam(p^2)] = 1)')
