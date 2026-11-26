# Numerical sanity checks for global.md
import mpmath as mp
mp.mp.dps = 25
pi = mp.pi
out = []
# (1) Poisson identity: (1/2pi) int 2w/(w^2+t^2) a(t) dt = psi(1/4+w/2) - log pi,  a(t) = Re psi(1/4+it/2) - log pi
def a(t): return mp.re(mp.digamma(mp.mpf(1)/4 + 1j*t/2)) - mp.log(pi)
for w in [mp.mpf('0.6'), mp.mpf(1), mp.mpf('2.5')]:
    lhs = mp.quad(lambda t: 2*w/(w*w+t*t)*a(t), [-mp.inf, -50, -5, 0, 5, 50, mp.inf])/(2*pi)
    rhs = mp.digamma(mp.mpf(1)/4 + w/2) - mp.log(pi)
    out.append("(1) w=%s  lhs=%s  rhs=%s  diff=%.2e" % (w, mp.nstr(lhs,15), mp.nstr(rhs,15), float(abs(lhs-rhs))))
# (2) corrected Lemma H' identity for zeta: sum_{n!=0} c sinc^2(pi c n) = 1 - c  (0<c<=1)
for c in [mp.mpf('0.3'), mp.mpf('0.77'), mp.mpf(1)]:
    S = 2*mp.nsum(lambda n: c*(mp.sin(pi*c*n)/(pi*c*n))**2, [1, mp.inf])
    out.append("(2) c=%s  sum=%s  1-c=%s" % (c, mp.nstr(S,12), mp.nstr(1-c,12)))
# (3) Hermite +1 eigenfunction F = h4 - 12.1 h0, h_k(x) = H_k(sqrt(2pi) x) e^{-pi x^2}
def F(x): return (mp.hermite(4, mp.sqrt(2*pi)*x) - mp.mpf('12.1'))*mp.exp(-pi*x*x)
def FT(xi): return 2*mp.quad(lambda x: F(x)*mp.cos(2*pi*x*xi), [0, 2, 6, mp.inf])
for xi in [0, mp.mpf('0.4'), 1, mp.mpf('1.7')]:
    out.append("(3) xi=%s  F=%s  FT(F)=%s" % (xi, mp.nstr(F(xi),12), mp.nstr(FT(xi),12)))
mn = min(F(1 + k*mp.mpf('0.01')) for k in range(0, 500))
out.append("(3) F(0)=%s ; min of F on [1,6] = %s (positive => F>=0 on |x|>=1)" % (mp.nstr(F(0),6), mp.nstr(mn,6)))
# (4) Herglotz form for zeta at real s: xi'/xi(s) = int w/(w^2+t^2) dmu_0 = 2 sum_{gamma>0} w/(w^2+gamma^2), w = s-1/2
def xilogder(s):
    return 1/s + 1/(s-1) - mp.log(pi)/2 + mp.digamma(s/2)/2 + mp.diff(lambda z: mp.log(mp.zeta(z)), s)
N = 400
gam = [mp.im(mp.zetazero(k)) for k in range(1, N+1)]
for s in [mp.mpf(2), mp.mpf('1.2')]:
    w = s - mp.mpf(1)/2
    part = 2*sum(w/(w*w+g*g) for g in gam)
    T = (gam[-1] + mp.im(mp.zetazero(N+1)))/2
    tail = 2*mp.quad(lambda t: w/(w*w+t*t)*mp.log(t/(2*pi))/(2*pi), [T, mp.inf])
    out.append("(4) s=%s  xi'/xi=%s  zero sum (400 zeros + smooth tail)=%s" % (s, mp.nstr(xilogder(s),10), mp.nstr(part+tail,10)))
print("\n".join(out))
