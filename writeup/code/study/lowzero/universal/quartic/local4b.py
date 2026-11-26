# closed form of the family-averaged local series E_l(X) = sum_v b_l(v) X^v and identities used by the prediction
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
from local4 import T, F, Z, l, X, J
b = {v: sp.cancel(sum(T[n][0] * F(n, v) for n in T) / Z) for v in range(-3, J)}
b0 = l**4 / ((l - 1) * (l**3 + l**2 + 2 * l + 1))
h = 1 - 1 / (X * l**2) + X**2 / l - X**4 / l - X**5 / l + X**5 / l**2 - X**6 / l + X**7 / l**2
Eclosed = b0 * (1 - X**2 / l) * h / ((1 - X**2) * (1 - X**3) * (1 - X**4))
ser = sp.expand(sp.series(Eclosed * X, X, 0, J - 2).removeO())
print('E_l(X) = b_l(0)(1 - X^2/l) h_l(X)/((1-X^2)(1-X^3)(1-X^4)) matches b_l(v), v=-1..%d:' % (J - 4), all(sp.cancel(ser.coeff(X, v + 1) - b[v]) == 0 for v in range(-1, J - 3)))
print('b_l(0) == 1/(1 + l^-2 - l^-3 - l^-4):', sp.cancel(b0 - 1 / (1 + l**-2 - l**-3 - l**-4)) == 0)
t = sp.symbols('t', positive=True)
print('h_l(l^-1/2) / (1 - t^3)(1 + t^4 - t^6 - t^8)  (t = l^-1/2):', sp.cancel(h.subs(X, t).subs(l, t**-2) / ((1 - t**3) * (1 + t**4 - t**6 - t**8))))
P = sp.expand((1 - X**2 / l) * h); print('P_l(X) = (1 - X^2/l) h_l(X) =', P)
# per-type identity: sum_v F(type, v) l^(-v/2) = l/(l-1)
from local4 import coeffs
ok = True
for n in T:
    mu, lam = coeffs(n); Lx = sum(lam[k] * X**k for k in range(J + 4)); Linv = sum(mu[k] * X**k for k in range(4))
    ok &= all(sp.cancel(F(n, v) - sp.expand(sp.Rational(1) * l / (l - 1) * sp.expand(Linv.subs(X, 1 / (l * X)) * Lx)).coeff(X, v)) == 0 for v in range(-3, 8))
print('per type: sum_v F(v) X^v = (l/(l-1)) L^-1(1/(lX)) L(X) (checked v=-3..7):', ok)
