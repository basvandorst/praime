# Referee: end-to-end check of the explicit formula normalisation and of Lemma 3 on true L-functions (PARI zeros).
# Own code: G by quadrature of the autocorrelations, Ghat by the closed form 2 sin(at) sum (-1)^k s_k/(t+om_k),
# F(kappa) by the u-form (Gauss-Legendre), c(n) from PARI's Dirichlet coefficients by the log-derivative recursion.
# For each L:  W(G) = c0 log N + sum_j F(kappa_j) - 2 sum_{n<32} Re c(n) n^{-1/2} G(log n)
#   (a) zero side sum_rho Ghat(gamma) (zeros |gamma| <= H from lfunzeros, H a zero of sin(aH), plus smooth tail and the
#       resonance of n = 32) must equal W(G);
#   (b) Lemma 3: W(G)/d <= B (needs no GRH); (c) the zeros with Ghat(gamma) < 0 lie in (-T, T).
# Usage: ref_zeros.py <json> <B>
import sys, json, numpy as np, mpmath as mp, cypari2
from scipy.special import psi as cpsi
from fractions import Fraction as Fr
pari = cypari2.Pari(); pari.allocatemem(500*10**6); pari.set_real_precision(38)

js = json.load(open(sys.argv[1])); Bcert = float(sys.argv[2])
T = float(Fr(js['T']))
S = np.array([[float(Fr(x)) for x in v] for v in js['s']]); R = np.array([[float(Fr(x)) for x in v] for v in js['r']])
n = S.shape[1]; K = (n - 1)//2; ks = np.arange(-K, K + 1); a = np.log(32.0)/2; om = ks*np.pi/a
xg, wg = np.polynomial.legendre.leggauss(400)
def trig(c, x, d=False): return np.exp(1j*np.outer(x, om)) @ (c*(1j*om if d else 1))
def G_quad(u):
    u = abs(u)
    if u >= 2*a: return 0.0
    lo, hi = -a, a - u; x = (hi - lo)/2*xg + (hi + lo)/2; w = (hi - lo)/2*wg; tot = 0j
    for s in S: tot += np.sum(w*trig(s, x + u)*np.conj(trig(s, x)))
    for r in R: tot += np.sum(w*trig(r, x + u, True)*np.conj(trig(r, x, True))) - T*T*np.sum(w*trig(r, x + u)*np.conj(trig(r, x)))
    return tot.real
def fh(c, t):
    t = np.atleast_1d(np.asarray(t, float)); z = om[None, :] + t[:, None]
    near = np.abs(z) < 1e-9
    val = np.where(near, a, np.sin(a*z)/np.where(near, 1, z))      # int_{-a}^{a} e^{izu} du / 2
    return 2*(val @ c)
def Ghat(t):
    t = np.atleast_1d(np.asarray(t, float))
    if len(t) > 20000: return np.concatenate([Ghat(t[i:i + 20000]) for i in range(0, len(t), 20000)])
    out = 0
    for sg in (1, -1):
        tt = sg*t; z = om[None, :] + tt[:, None]; near = np.abs(z) < 1e-9
        M = 2*np.where(near, a, np.sin(a*z)/np.where(near, 1, z))
        fs = M @ S.T; fr = M @ R.T
        out = out + (np.sum(fs**2, axis=1) + (tt**2 - T*T)*np.sum(fr**2, axis=1))/2
    return out
c0 = G_quad(0.0)
U, UW = np.polynomial.legendre.leggauss(200); U = (U + 1)*a; UW = UW*a
GU = np.array([G_quad(u) for u in U]); EG = float(mp.euler)
def F(kap):
    integ = 2*(np.exp(-2*U)*c0 - np.exp(-(0.5 + kap)*U)*GU)/(-np.expm1(-2*U))
    return -(EG + np.log(np.pi))*c0 - c0*np.log1p(-2.0**-10) + np.sum(UW*integ)
Glog = {m: G_quad(np.log(m)) for m in range(2, 33)}
print('c0 = %.13e, G(log 32) = %.2e' % (c0, Glog[32]))

def smooth_density(t, N, kaps):
    return (np.log(N) + sum(-np.log(np.pi) + np.real(cpsi((0.5 + k + 1j*np.asarray(t))/2)) for k in kaps))/(2*np.pi)
_cache = {}
def grid(H, tmax=20000.0):
    if H not in _cache:
        t = np.linspace(H, tmax, int((tmax - H)*40) + 1); _cache[H] = (t, Ghat(t))
    return _cache[H]
def tail_smooth(H, N, kaps, tmax=20000.0):
    # (both signs) 2 * int_H^inf Ghat(t) * density(t) dt ; density interpolated on a coarse grid (smooth), Ghat exact
    t, gt = grid(H, tmax)
    d = smooth_density(t, N, kaps)
    val = np.trapezoid(gt*d, t)
    # beyond tmax: Ghat = 4 sin^2(at) Q(t), mean 2Q(t), Q ~ q/t^2 ; density ~ (log N + d log(t/(2 pi)))/(2 pi)
    tl = np.linspace(tmax - 2*np.pi/a*50, tmax, 20001); q = np.trapezoid(tl**2*Ghat(tl), tl)/(tl[-1] - tl[0])  # mean of t^2 Ghat = 2q
    dd = len(kaps)
    rest = q*(np.log(N) + dd*(np.log(tmax/(2*np.pi)) + 1))/(2*np.pi)/tmax
    return 2*(val + rest), q
def tail_res32(H, c32, tmax=20000.0):
    # n = 32 term of the zero density: -(1/pi) Re(c(32) 32^{-1/2} e^{-it log 32}) per unit t (each sign)
    t, gt = grid(H, tmax)
    v = -(1/np.pi)*np.real(c32*32**-0.5*np.exp(-1j*t*np.log(32)))
    return 2*np.trapezoid(gt*v, t)

def run(name, expr, H=200.0, selfdual=True):
    L = pari(expr); N, k, Vga = pari.lfunparams(L); N = float(N); k = float(k)
    kaps = [float(v) + (k - 1)/2 for v in Vga]; d = len(kaps)
    an = [complex(x) for x in pari.lfunan(L, 64)]
    an = [an[i]/(i + 1)**((k - 1)/2) for i in range(64)]           # analytic normalisation
    c = [0j]*65
    for m in range(2, 65):
        c[m] = an[m - 1]*np.log(m) - sum(c[dv]*an[m//dv - 1] for dv in range(2, m) if m % dv == 0)
    W = c0*np.log(N) + sum(F(kp) for kp in kaps) - 2*sum((c[m]*m**-0.5).real*Glog[m] for m in range(2, 32))
    # H: a zero of sin(aH) (Ghat and Ghat' vanish there)
    H = np.pi/a*np.ceil(H*a/np.pi)
    if selfdual:
        z = np.array([float(x) for x in pari.lfunzeros(L, H)])
        zs = np.concatenate([-z[z > 1e-9], z])
    else:
        zs = np.array([float(x) for x in pari.lfunzeros(L, [-H, H])])
    Zs = float(np.sum(Ghat(zs)))
    ts, q = tail_smooth(H, N, kaps); tr = tail_res32(H, c[32])
    neg = zs[Ghat(zs) < 0]
    print('%-14s N=%g kappa=%s d=%d | W = %.10e  zeros(|g|<=%.1f, %d) = %.10e  +smooth tail %.3e +res32 %.3e -> diff %.2e | W/d = %.4e <= B=%.4e: %s | zeros with Ghat<0: %s'
          % (name, N, [round(x, 3) for x in kaps], d, W, H, len(zs), Zs, ts, tr, Zs + ts + tr - W, W/d, Bcert, W/d <= Bcert,
             np.round(neg, 4).tolist()), flush=True)
    # local roots at 2 (|alpha| <= 1 ?) from the power sums c(2^m)/log 2
    return W

print('Ghat min on [T, 300] sampled: %.3e' % Ghat(np.linspace(T, 300, 300001)).min())
run('Delta (k=12)', 'lfunmf(mfinit([1,12],1), mfbasis(mfinit([1,12],1))[1])')
run('k=16 level 1', 'lfunmf(mfinit([1,16],1), mfbasis(mfinit([1,16],1))[1])')
run('11a (k=2)', 'lfuncreate(ellinit([0,-1,1,-10,-20]))')
run('37a (rank 1)', 'lfuncreate(ellinit([0,0,1,-1,0]))')
run('chi_-3', 'lfuncreate(-3)')
run('chi_-4', 'lfuncreate(-4)')
run('chi_5 (real)', 'lfuncreate(5)')
run('chi mod 5 ord4', 'lfuncreate([znstar(5,1),[1]])', selfdual=False)
run('chi mod 7 ord3', 'lfuncreate([znstar(7,1),[2]])', selfdual=False)
run('Sym2 11a', 'lfunsympow(ellinit([0,-1,1,-10,-20]),2)')
