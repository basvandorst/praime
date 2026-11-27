# Referee, independent float/mpmath evaluation of the certificate (own code; does not import the repo).
# G is built by direct Gauss-Legendre quadrature of the autocorrelations f*f~ (not from a closed form).
# F(kappa) by two routes: (i) u-form derived from Gauss' integral for digamma, (ii) Fourier side with mpmath digamma.
# Usage: ref_float.py <json>
import sys, json, numpy as np
from fractions import Fraction as Fr
import mpmath as mp

js = json.load(open(sys.argv[1]))
T = float(Fr(js['T']))
S = np.array([[float(Fr(x)) for x in v] for v in js['s']])
R = np.array([[float(Fr(x)) for x in v] for v in js['r']])
n = S.shape[1]; K = (n - 1)//2
ks = np.arange(-K, K + 1)
a = np.log(32.0)/2
om = ks*np.pi/a

xg, wg = np.polynomial.legendre.leggauss(400)

def trig(c, x, deriv=False):
    """sum_k c_k e^{i om_k x} (or its derivative)"""
    E = np.exp(1j*np.outer(x, om))
    return E @ (c*(1j*om if deriv else 1))

def G_quad(u):
    """G(u) = Re[sum_i f_i*f_i~(u) + sum_j (psi_j'*psi_j'~ - T^2 psi_j*psi_j~)(u)], u in [0,2a], by quadrature."""
    u = abs(u)
    if u >= 2*a: return 0.0
    lo, hi = -a, a - u
    x = (hi - lo)/2*xg + (hi + lo)/2; w = (hi - lo)/2*wg
    tot = 0j
    for s in S:
        tot += np.sum(w*trig(s, x + u)*np.conj(trig(s, x)))
    for r in R:
        tot += np.sum(w*trig(r, x + u, True)*np.conj(trig(r, x, True)))
        tot -= T*T*np.sum(w*trig(r, x + u)*np.conj(trig(r, x)))
    return tot.real

# boundary values of psi_j (side condition)
print('psi_j(+-a):', [abs(trig(r, np.array([a]))[0]) for r in R])
c0 = G_quad(0.0)
print('c0 = G(0) = %.15e' % c0)
for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
    print('G(m log %d):' % p, ' '.join('%.3e' % G_quad(m*np.log(p)) for m in range(1, 6) if m*np.log(p) < 2*a))

# Fourier transform, closed form derived independently: fhat(t) = int_{-a}^{a} f(u) e^{itu} du
def fhat(c, t):
    t = np.asarray(t, float)
    # sum_k c_k * 2 sin((om_k + t) a)/(om_k + t); handle t = -om_k
    z = (om[None, :] + t[:, None])*a
    sinc = np.where(np.abs(z) < 1e-12, 1.0, np.sin(z)/np.where(z == 0, 1, z))
    return (2*a*sinc) @ c
def Ghat(t):
    t = np.atleast_1d(np.asarray(t, float))
    if len(t) > 50000:
        return np.concatenate([Ghat(t[i:i + 50000]) for i in range(0, len(t), 50000)])
    out = np.zeros(len(t))
    for sgn in (1, -1):
        tt = sgn*t
        gh = sum(np.abs(fhat(s, tt))**2 for s in S) + (tt**2 - T*T)*sum(np.abs(fhat(r, tt))**2 for r in R)
        out += gh/2
    return out
tt = np.linspace(T, 400, 400001)
gv = Ghat(tt)
print('min Ghat on [T,400] (sampled 400001): %.3e ; Ghat(0) = %.6e' % (gv.min(), Ghat(0.0)[0]))
# Fourier inversion check of G(0): (1/2pi) int Ghat
ug, uw = np.polynomial.legendre.leggauss(200)

# u-form of F
U = (ug + 1)*a; UW = uw*a          # nodes on [0, 2a]
GU = np.array([G_quad(u) for u in U])
EG = 0.5772156649015329
def F_u(kap):
    integ = 2*(np.exp(-2*U)*c0 - np.exp(-(0.5 + kap)*U)*GU)/(-np.expm1(-2*U))
    return -(EG + np.log(np.pi))*c0 - c0*np.log1p(-2.0**-10) + np.sum(UW*integ)
def F_u_quad(kap):
    # second independent quadrature: mpmath quad on G_quad (slow) -- only a few points
    f = lambda u: 2*(mp.e**(-2*u)*c0 - mp.e**(-(0.5 + kap)*u)*G_quad(float(u)))/(1 - mp.e**(-2*u))
    return float(-(EG + mp.log(mp.pi))*c0 - c0*mp.log(1 - mp.mpf(2)**-10) + mp.quad(f, [0, 0.5, 1, 2, 2*a]))
# Fourier route: F(kappa) = (1/pi) int_0^inf Ghat(t) [-log pi + Re psi((1/2+kappa+it)/2)] dt
def F_four(kap, tmax=3000.0, nn=600001):
    t = np.linspace(0, tmax, nn)
    gh = Ghat(t)
    ps = np.array([float(mp.re(mp.digamma((0.5 + kap + 1j*x)/2))) for x in t[::50]])
    # digamma is smooth: interpolate on the coarse grid
    psf = np.interp(t, t[::50], ps)
    val = np.trapezoid(gh*(-np.log(np.pi) + psf), t)/np.pi
    # tail: Ghat ~ 4 sin^2(at) Q(t), average 2 Q(t); Q ~ q2/t^2; Re psi ~ log(t/2)
    q2 = (Ghat(np.array([tmax*1.0]))[0])  # not used directly
    # estimate tail by averaging Ghat over periods: use mean of t^2 Ghat on last stretch
    tl = t[t > tmax - 200]; m2 = np.mean(tl**2*Ghat(tl))
    tail = m2*(np.log(tmax/2) + 1 - np.log(np.pi))/tmax/np.pi
    return val + tail
print('check (1/2pi) int Ghat vs c0:', np.trapezoid(Ghat(np.linspace(0, 3000, 600001)), np.linspace(0, 3000, 600001))/np.pi, c0)
for kap in [0.0, 1.0, 5.5, 5.8738, 6.5, 20.0]:
    print('F(%g): u-form GL %.12e   mpmath quad %.12e' % (kap, F_u(kap), F_u_quad(kap)), flush=True)
for kap in [0.0, 5.8738]:
    print('F(%g): Fourier route %.9e' % (kap, F_four(kap)), flush=True)

# sup of F over kappa >= 0 (float scan; rigorous version in ref_arb.py)
kg = np.arange(0, 60, 0.001)
Fv = np.array([F_u(k) for k in kg])
i = int(np.argmax(Fv)); print('float sup F on [0,60]: %.12e at %.4f ; F(60) = %.4e, F(1e3) = %.4e, F(-0.4999) = %.4e' % (Fv[i], kg[i], F_u(60.0), F_u(1e3), F_u(-0.4999)))
km = np.arange(-0.4999, 0, 0.0005); print('max F on (-1/2, 0]: %.6e' % max(F_u(k) for k in km))
supF = Fv[i]

# primes: w_{p,m} = -2 log p p^{-m/2} G(m log p); M_p(R) = max_theta sum_m w R^m cos m theta
th = np.linspace(0, np.pi, 200001)
def Mp(p, R):
    ws = [(m, -2*np.log(p)*p**(-m/2)*G_quad(m*np.log(p))) for m in range(1, 8) if m*np.log(p) < 2*a - 1e-12]
    P = sum(w*R**m*np.cos(m*th) for m, w in ws)
    j = int(np.argmax(P))
    # refine
    f = lambda x: float(sum(w*R**m*mp.cos(m*x) for m, w in ws))
    xm = float(mp.findroot(lambda x: sum(-m*w*R**m*mp.sin(m*x) for m, w in ws), th[j])) if 0 < j < len(th) - 1 else th[j]
    return max(P[j], f(xm)), th[j], ws
primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
for label, Rs in [('Ramanujan all', {p: 1.0 for p in primes}), ('Ram at 2, R=p else', {p: (1.0 if p == 2 else float(p)) for p in primes}),
                  ('trivial R=p everywhere', {p: float(p) for p in primes})]:
    tot = 0.0; parts = []
    for p in primes:
        m_, t_, ws = Mp(p, Rs[p]); tot += m_; parts.append('%d:%.4e' % (p, m_))
    print('%s: sum M_p = %.10e ; B = supF + sum = %.10e' % (label, tot, supF + tot))
    print('   ', ' '.join(parts))
print('w_{2,m}:', Mp(2, 1.0)[2])
