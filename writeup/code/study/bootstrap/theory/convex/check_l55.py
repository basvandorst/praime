import mpmath as mp
mp.mp.dps = 25
for Z in [0.4, 1.3, 3.0]:
    lhs = mp.quad(lambda t: mp.besselk(0, Z*t)*mp.besselk(0, Z/t), [0, 1, mp.inf])
    print("Z=%g  int K0(Zt)K0(Z/t)dt = %s   (pi/Z)K0(2Z) = %s" % (Z, mp.nstr(lhs, 18), mp.nstr(mp.pi/Z*mp.besselk(0, 2*Z), 18)))
# second identity: int_0^inf v^-2 K0(2 pi n / v) K0(2 pi x v) dv = K0(4 pi sqrt(n x))/(2n)
for n, x in [(1, 1.5), (2, 1.2), (3, 2.7)]:
    lhs = mp.quad(lambda v: v**-2*mp.besselk(0, 2*mp.pi*n/v)*mp.besselk(0, 2*mp.pi*x*v), [0, 1, mp.inf])
    print("n=%d x=%g  lhs = %s  rhs = %s" % (n, x, mp.nstr(lhs, 15), mp.nstr(mp.besselk(0, 4*mp.pi*mp.sqrt(n*x))/(2*n), 15)))
