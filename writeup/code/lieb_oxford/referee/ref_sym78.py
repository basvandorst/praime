"""Point 3(f): one configuration of (78) re-derived symbolically by the referee's own route
(shells of ball 2 about e1 + Archimedes + Newton potential of ball 1), partial overlap with R1 <= 1, R2 <= 1,
R1 + R2 > 1.  Also the region-A closed form and the sphere formula (43) in the overlap configuration."""
import sympy as sp
R1, R2, rho, t, y = sp.symbols('R1 R2 rho t y', positive=True)
# antiderivative in t of V1(sqrt t), V1 = ball potential
P_in = (3 * R1**2 * t - t**2 / 2) / (2 * R1**3)
P_out = 2 * sp.sqrt(t) - sp.Rational(3, 4) * R1
# shell radius rho about e1: |y|^2 uniform on [(1-rho)^2, (1+rho)^2]
# piece 1: rho in [0, 1-R1]: whole shell outside ball 1 -> avg = (P_out((1+rho)^2) - P_out((1-rho)^2))/(4 rho) = 1
avg1 = sp.simplify((P_out.subs(t, (1 + rho)**2) - P_out.subs(t, (1 - rho)**2)) / (4 * rho))
# piece 2: rho in [1-R1, R2]: (1-rho)^2 <= R1^2 <= (1+rho)^2 -> split at R1^2
avg2 = ((P_out.subs(t, (1 + rho)**2) - P_out.subs(t, R1**2)) + (P_in.subs(t, R1**2) - P_in.subs(t, (1 - rho)**2))) / (4 * rho)
avg2 = sp.simplify(avg2.subs(sp.sqrt((1 + rho)**2), 1 + rho))
W = sp.integrate(3 * rho**2 / R2**3 * avg1, (rho, 0, 1 - R1)) + sp.integrate(3 * rho**2 / R2**3 * avg2, (rho, 1 - R1, R2))
a, b = 1 / R1, 1 / R2
Psi = sp.simplify(a**3 * b**3 * (1 - W))
T1 = (a + b - a * b)**4 / 160 * (a**2 * b**2 - 5 * a**2 - 5 * b**2 + 4 * a * b**2 + 4 * b * a**2 + 20 * a * b)
print("avg on piece 1 =", avg1)
print("partial overlap, R1<=1, R2<=1: Psi - T1 of (78) simplifies to", sp.simplify(sp.expand(Psi - T1)))
# region A (ball 1 inside ball 2, R1 + 1 <= R2): all shells of ball 2 with rho < 1 - R1 see avg 1 ... use direct Newton:
# average over ball 1 (radius R1 about 0) of the potential of ball 2 (radius R2 about e1), |x - e1| <= 1 + R1 <= R2
xs, rr = sp.symbols('xs rr', positive=True)
V2 = lambda s2: (3 * R2**2 - s2) / (2 * R2**3)        # s2 = |x-e1|^2 inside ball 2
# average of |x - e1|^2 over uniform ball radius R1 about 0 = 1 + 3 R1^2/5
WA = V2(1 + sp.Rational(3, 5) * R1**2)
PA_mine = sp.expand(a**3 * b**3 * (1 - WA))
aa, bb = sp.symbols('a b', positive=True)
PA_note = aa**3 * bb**3 * (1 - sp.Rational(3, 2) * bb + bb**3 / 2) + sp.Rational(3, 10) * aa * bb**6
print("region A: mine - note PA =", sp.simplify(PA_mine.subs({R1: 1 / aa, R2: 1 / bb}) - PA_note))
