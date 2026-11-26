# Consistency check of Part II Thm 5.4(c),(d) at x = 1: the c=1 part of a_n K_0(4 pi sqrt(n x)) sqrt x.
# n(-I2dot(z))K0(z) * 16 pi^2 = z^2 (-I2dot(z)) K0(z); compare with 1 - 3/(2z) + c3/z^3 + ...
import mpmath as mp
mp.mp.dps = 40
def mI2d(z): return mp.besselk(2,z) + 2/z*mp.besseli(1,z) - 2/z**2*mp.besseli(0,z)
for z in [10, 20, 40, 80, 160, 320]:
    z = mp.mpf(z)
    v = z**2*mI2d(z)*mp.besselk(0,z)
    r = v - (1 - mp.mpf(3)/(2*z))
    print("z=%4d  z^2(-I2dot)K0 = %s   residual*z^2 = %s  residual*z^3 = %s" % (int(z), mp.nstr(v,20), mp.nstr(r*z**2,10), mp.nstr(r*z**3,10)))
