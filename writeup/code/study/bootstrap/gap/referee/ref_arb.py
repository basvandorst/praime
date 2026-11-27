# Referee, rigorous (python-flint arb) evaluation of B, by a method different from the proof's code:
#   * F(kappa) = C0 + int_0^{2a} phi_kappa(u) du by N-point Gauss-Legendre with certified nodes (arb.legendre_p_root)
#     and the Bernstein-ellipse error bound (Trefethen, SIAM Rev. 2008, Thm 4.5; used with the weaker rho^(-2N+2)),
#     M = max |phi| on the ellipse boundary from acb balls covering the boundary (maximum modulus principle).
#   * sup over kappa-intervals by the mean value form F(K) in F(k_c) + F'(K)[-r, r], F'(K) by the same quadrature.
#   * tail kappa >= K1: F' < 0 from h(u) = 2u e^{-u/2} G(u)/(1-e^{-2u}) <= -eta on [0,u1], |h| <= Mh (interval bounds).
#   * M_p: polynomial in x = cos(theta), grid + second-order remainder, all in arb.
# Every float that enters a bound is an arb ball; no float margins.
# Usage: ref_arb.py <json> [kmax_tail]
import sys, json
from fractions import Fraction as Fr
import flint
from flint import arb, acb, ctx
ctx.prec = 128

js = json.load(open(sys.argv[1]))
T = arb(flint.fmpq(*Fr(js['T']).as_integer_ratio()))
def q(x): f = Fr(x); return arb(flint.fmpq(f.numerator, f.denominator))
S = [[q(x) for x in v] for v in js['s']]
R = [[q(x) for x in v] for v in js['r']]
n = len(S[0]); K = (n - 1)//2
ks = list(range(-K, K + 1))
pi = arb.pi(); a = arb(5)/2*arb(2).log(); twoa = 2*a
om = [k*pi/a for k in ks]
# D = sum s s^T + (pi/a)^2 K R K - T^2 R
D = [[arb(0)]*n for _ in range(n)]
for i in range(n):
    for j in range(n):
        ss = sum((v[i]*v[j] for v in S), arb(0)); rr = sum((v[i]*v[j] for v in R), arb(0))
        D[i][j] = ss + (pi/a)**2*ks[i]*ks[j]*rr - T*T*rr
diag = [D[i][i] for i in range(n)]
e = [sum((D[l][k]*(-1)**((ks[l] - ks[k]) % 2)/(ks[l] - ks[k]) for l in range(n) if l != k), arb(0)) for k in range(n)]
c0 = twoa*sum(diag, arb(0))
print('c0 =', c0)

def G(z):
    """right branch (0 <= u <= 2a) of G, analytic continuation; z arb or acb"""
    if isinstance(z, acb):
        return (twoa - z)*sum((diag[i]*(om[i]*z).cos() for i in range(n)), acb(0)) + twoa/pi*sum((e[i]*(om[i]*z).sin() for i in range(n)), acb(0))
    return (twoa - z)*sum((diag[i]*(om[i]*z).cos() for i in range(n)), arb(0)) + twoa/pi*sum((e[i]*(om[i]*z).sin() for i in range(n)), arb(0))

# cross-check G with the float quadrature values of ref_float.py (printed there)
for u in ['0.6931471805599453', '1.0986122886681098', '1.3862943611198906', '2.0794415416798357', '2.5']:
    print('G(%s) = %s' % (u, G(arb(u)).str(12)))
print('G(2a) =', G(twoa).str(5))

EG = arb.const_euler()
C0 = -(EG + pi.log())*c0 - c0*(1 - arb(2)**-10).log()
def phi(z, kap):
    one = acb(1) if isinstance(z, acb) else arb(1)
    return 2*((-2*z).exp()*c0 - (-(arb('0.5') + kap)*z).exp()*G(z))/(one - (-2*z).exp())
def dphi(z, kap):      # d/dkappa of -(...)  -> F'(kappa) integrand
    one = acb(1) if isinstance(z, acb) else arb(1)
    return 2*z*(-(arb('0.5') + kap)*z).exp()*G(z)/(one - (-2*z).exp())

# Gauss-Legendre nodes on [0, 2a]
N = int(sys.argv[3]) if len(sys.argv) > 3 else 120
nodes = []
for k in range(N):
    x, w = arb.legendre_p_root(N, k, weight=True)
    nodes.append((a + a*x, a*w))
Gn = [G(x) for x, w in nodes]
E2 = [(-2*x).exp() for x, w in nodes]
den = [1 - v for v in E2]

# Bernstein ellipse: foci 0, 2a; semi-minor axis 1.0  ->  rho
b = arb(1)
rho = (b + (b*b + a*a).sqrt())/a
print('rho =', rho.str(8), ' semi-major extends left of 0 by', (a*(rho + 1/rho)/2 - a).str(6))
def ellipse_M(f, kap_ball, nseg=3000):
    Mx = arb(0)
    twopi = 2*pi
    for j in range(nseg):
        th = (arb(j) + arb('0.5') + arb(0, '0.5'))*twopi/nseg       # ball covering [j, j+1]*2pi/nseg
        ei = acb(th.cos(), th.sin())
        t = (rho*ei + ei.conjugate()/rho)/2
        z = a + a*t
        v = abs(f(z, kap_ball))
        Mx = Mx.union(v) if j else v
    return arb(Mx.upper())
def quad_err(Mb):
    return a*arb(64)/15*Mb*rho**(-2*N + 2)/(rho*rho - 1)

def F_point(kap):
    s = sum((w*2*(E2[i]*c0 - (-(arb('0.5') + kap)*x).exp()*Gn[i])/den[i] for i, (x, w) in enumerate(nodes)), arb(0))
    return C0 + s
def dF_ball(kb):
    return sum((w*2*x*(-(arb('0.5') + kb)*x).exp()*Gn[i]/den[i] for i, (x, w) in enumerate(nodes)), arb(0))

# tail: h(u) = 2u e^{-u/2} G(u)/(1-e^{-2u}); fac(u) = 2u/(1-e^{-2u}) increasing, >= 1.
# F'(kappa) = int_0^{2a} e^{-kappa u} h(u) du <= [-eta (1 - e^{-kappa u1}) + Mh e^{-kappa u1}]/kappa  if h <= -eta on [0,u1], |h| <= Mh on [u1,2a]
def fac_hi(u):
    return (2*u/(1 - (-2*u).exp())).upper() if u > 0 else arb(1)
def G_ball(ul, ur):
    c = (ul + ur)/2; r = (ur - ul)/2
    return G(c + arb(0, r.upper()))
npc = 4096; step = twoa/npc
gb = [G_ball(step*i, step*(i + 1)) for i in range(npc)]
habs = [arb(fac_hi(step*(i + 1)))*abs(gb[i]) for i in range(npc)]
# suffix max of |h|
suf = [None]*(npc + 1); suf[npc] = arb(0)
for i in range(npc - 1, -1, -1): suf[i] = arb(arb(suf[i + 1].upper()).union(arb(habs[i].upper())).upper())
best = None; cur_eta = None
for i in range(npc):
    if not gb[i].upper() < 0: break
    et = arb(((-(step*(i + 1)/2)).exp()*(-gb[i].upper())).lower())
    cur_eta = et if cur_eta is None else arb(min(cur_eta.lower(), et.lower()))
    u1 = step*(i + 1); Mh = suf[i + 1]
    Kt = (1 + Mh/cur_eta).log()/u1
    if best is None or Kt.upper() < best[0].upper(): best = (Kt, u1, cur_eta, Mh)
Kt, u1, eta, Mh = best
print('tail: h <= -eta on [0,u1], u1 = %s, eta >= %s; |h| <= %s on [u1,2a]; F decreasing for kappa > %s' % (u1.str(6), eta.str(6), Mh.str(6), Kt.str(6)), flush=True)
KTAIL = int(Kt.upper().floor().unique_fmpz()) + 1 if hasattr(Kt.upper().floor(), 'unique_fmpz') else int(float(Kt.upper())) + 1

KLO, KHI = arb('-0.5'), arb(KTAIL)
kall = KLO.union(KHI)
M_phi = ellipse_M(phi, kall); M_dphi = ellipse_M(dphi, kall)
errF = quad_err(M_phi); errdF = quad_err(M_dphi)
print('Bernstein: M_phi <= %s, M_dphi <= %s ; quadrature errors <= %s, %s' % (M_phi.str(5), M_dphi.str(5), errF.str(3), errdF.str(3)))
errF = arb(0, errF.upper()); errdF = arb(0, errdF.upper())

for kap in ['0', '1', '5.5', '5.8738', '6.5', '-0.4']:
    print('F(%s) = %s' % (kap, (F_point(arb(kap)) + errF).str(15)), flush=True)

# sup over [-1/2, 12] by adaptive intervals
def upper_on(kl, kr):
    kc = (kl + kr)/2; r = (kr - kl)/2
    Fc = F_point(kc) + errF
    dF = dF_ball(kc + arb(0, r.upper())) + errdF
    return (Fc + dF*arb(0, r.upper())).upper()
stack = [(arb(i)/16, arb(i + 1)/16) for i in range(-8, KTAIL*16)]
best_lower = max(F_point(arb(x)/1000).lower() for x in range(5850, 5900))
target = best_lower + arb('1e-8')
supU = arb(-1); cnt = 0; worst = None
while stack:
    kl, kr = stack.pop(); cnt += 1
    u = upper_on(kl, kr)
    if u > target and (kr - kl) > arb('1e-7'):
        m = (kl + kr)/2; stack += [(kl, m), (m, kr)]; continue
    if u > supU: supU = u; worst = (kl, kr)
print('sup F on [-1/2,%d] <=' % KTAIL + ' %s (intervals evaluated %d; attained near [%s,%s]); best point value >= %s' % (supU.str(12), cnt, worst[0].str(8), worst[1].str(8), best_lower.str(12)))

# primes
primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
def wpm(p, m):
    lp = arb(p).log()
    return -2*lp*arb(p)**(-arb(m)/2)*G(m*lp)
def Mp(p, Rr):
    ws = [(m, wpm(p, m)*Rr**m) for m in range(1, 6) if m*arb(p).log() < twoa]
    # Chebyshev: cos(m th) = T_m(x)
    Tm = {1: [0, 1], 2: [-1, 0, 2], 3: [0, -3, 0, 4], 4: [1, 0, -8, 0, 8], 5: [0, 5, 0, -20, 0, 16]}
    co = [arb(0)]*6
    for m, w in ws:
        for j, c in enumerate(Tm[m]): co[j] = co[j] + w*c
    P = lambda x: sum((co[j]*x**j for j in range(6)), arb(0))
    dP = lambda x: sum((j*co[j]*x**(j - 1) for j in range(1, 6)), arb(0))
    d2 = sum((abs(co[j])*j*(j - 1) for j in range(2, 6)), arb(0))*arb('1.001')   # |P''| on [-1-h,1+h]
    ng = 20000; h = arb(2)/ng; best = None
    for i in range(ng + 1):
        x = -1 + h*i
        v = P(x) + abs(dP(x))*h/2 + d2*h*h/8
        best = v if best is None else best.union(v)
    return arb(best.upper())
for label, Rs in [('Ramanujan at all p<=31', {p: arb(1) for p in primes}), ('Ramanujan at 2 only, R=p else', {p: (arb(1) if p == 2 else arb(p)) for p in primes}),
                  ('R2=1.0365, others 1', {p: (arb('1.0365') if p == 2 else arb(1)) for p in primes}), ('R=p everywhere', {p: arb(p) for p in primes}),
                  ('Ram at 2 and 5, R=p else', {p: (arb(1) if p in (2, 5) else arb(p)) for p in primes}), ('Ram at 2, R=sqrt(p) else', {p: (arb(1) if p == 2 else arb(p).sqrt()) for p in primes})]:
    tot = arb(0); parts = []
    for p in primes:
        m_ = Mp(p, Rs[p]); tot += m_; parts.append('%d:%s' % (p, m_.str(4)))
    print('%s: sum M_p <= %s ; B <= %s' % (label, tot.str(10), (supU + tot).str(10)))
    print('    ' + ' '.join(parts), flush=True)
