# Continuous Lambda-relaxation test: dLambda = e^u du on [u0, inf) (psi(x) = x - x0).
# Feasible iff  Re H(1/2+it) = b(t) - F(t) >= 0 for all t, where
# b(t) = (1/2)(Re psi(1/4+it/2) - log pi) + (1/2)/(1/4+t^2),  F(t) = Re L_1(1/2+it),
# L_1(s) = -int_0^{u0} e^{-(s-1)u} du.
import mpmath as mp
mp.mp.dps = 30
def b(t):
    return 0.5*(mp.re(mp.digamma(mp.mpf(1)/4 + 1j*t/2)) - mp.log(mp.pi)) + 0.5/(0.25 + t*t)
def F(t, u0):
    s = mp.mpf(1)/2 + 1j*t
    return mp.re(-(1 - mp.exp(-(s-1)*u0))/(s-1))
for x0 in [1.5, 1.8, 2.0, 2.2, 2.5, 3.0, 4.0]:
    u0 = mp.log(x0)
    worst = min((b(t) - F(t, u0), t) for t in [k*0.01 for k in range(0, 4000)])
    print("x0=%.2f  min_t [b-F] = %.5f at t=%.2f" % (x0, worst[0], worst[1]))
print("b(t) at t=0..14:", [round(float(b(t)),3) for t in range(15)])
