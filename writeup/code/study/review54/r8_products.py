# Products zeta_K = zeta * L(rho): finite-height terms from the lines of L(rho) swept by the zeta chirp.
# Full stationary-phase sum over all rationals r = a/b (not only r = 1/k), exact integrals, own line factors.
import os, sys, math, json, numpy as np, cypari2
from scipy.special import loggamma
from fractions import Fraction
from math import gcd
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6)
D = os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/universal/'
T = 300.0

def th(t, Q, mu):
    t = np.asarray(t, float); d = len(mu)
    return 0.5 * t * math.log(Q / math.pi ** d) + sum(np.imag(loggamma((0.5 + m + 1j * t) / 2)) for m in mu)

def fac(n):
    f = {}; p = 2
    while p * p <= n:
        while n % p == 0: f[p] = f.get(p, 0) + 1; n //= p
        p += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

for key, pol, dd, c in (('R23', 'x^3-x-1', 23, 1.853352), ('R31', 'x^3+x-1', 31, 1.879242)):
    f = pari(pol)
    cache = {}
    def local(p):
        if p in cache: return cache[p]
        if p == dd: mu = [1, -1]; lam = lambda k: 1.0
        else:
            nr = len(pari.polrootsmod(f, p))
            if nr == 3: mu = [1, -2, 1]; lam = lambda k: k + 1.0
            elif nr == 0: mu = [1, 1, 1]; lam = lambda k: [1.0, -1.0, 0.0][k % 3]
            else: mu = [1, 0, -1]; lam = lambda k: 1.0 if k % 2 == 0 else 0.0
        M = lambda h: mu[h] if h < len(mu) else 0.0
        Fp = sum(lam(k) * M(k) * p ** -k for k in range(4))
        def rho(v):
            if v >= 0: return sum(lam(v + h) * M(h) * p ** -h for h in range(60)) / Fp
            return sum(lam(m) * M(-v + m) * p ** -m for m in range(60)) / Fp
        cache[p] = (rho, lam, nr if p != dd else -1)
        return cache[p]
    def rho_r(a, b):
        v = 1.0
        for p, e in fac(a).items(): v *= local(p)[0](e)
        for p, e in fac(b).items(): v *= local(p)[0](-e)
        return v
    # psi_rho(x)
    x = T / (2 * math.pi)
    psi = 0.0
    for n in range(2, int(x) + 1):
        fa = fac(n)
        if len(fa) != 1: continue
        p, k = list(fa.items())[0]
        rh, lam, cls = local(p)
        if cls == -1: s = 1.0
        elif cls == 3: s = 2.0
        elif cls == 0: s = -1.0 if k % 3 else 2.0
        else: s = 0.0 if k % 2 else 2.0
        psi += s * math.log(p)
    # zeros
    gr = np.load(D + 'zeros/zeros_%s.npz' % key)['g']; gr = np.sort(gr[(gr > 0) & (gr <= T)])
    gz = np.load(D + 'zeros/zeros_Zeta.npz')['g']; gz = np.sort(gz[(gz > 0) & (gz <= T)])
    thz = lambda t: th(t, 1.0, [0.0]); thr = lambda t: th(t, float(dd), [0.0, 1.0])
    S_rho_zeta = np.exp(2j * thz(gr)).sum(); S_rho_K = np.exp(2j * (thz(gr) + thr(gr))).sum()
    S_zeta_zeta = np.exp(2j * thz(gz)).sum(); S_zeta_K = np.exp(2j * (thz(gz) + thr(gz))).sum()
    S_rho_rho = np.exp(2j * thr(gr)).sum()
    # line sum: (1/2pi) int_0^T e^{2i theta_zeta(t)} A_r r^{it} dt, A_r = -c rho(r)/sqrt(ab)
    tt = np.linspace(1e-6, T, 300001); dt = tt[1] - tt[0]
    E = np.exp(2j * thz(tt))
    def integ(r):
        fv = E * np.exp(1j * tt * math.log(r))
        return ((fv.sum() - 0.5 * (fv[0] + fv[-1])) * dt) / (2 * math.pi)
    tot_a1 = 0j; tot_all = 0j; by_a = {}
    for a in range(1, 13):
        for b in range(1, int(a * x * 1.25) + 3):
            if gcd(a, b) != 1: continue
            r = a / b
            A = -c * rho_r(a, b) / math.sqrt(a * b)
            v = A * integ(r)
            by_a[a] = by_a.get(a, 0j) + v
            tot_all += v
            if a == 1: tot_a1 += v
    # the simple a = 1 stationary formula
    simple = -c * sum(rho_r(1, k) for k in range(1, int(x) + 1))
    print('%s: x = T/2pi = %.2f ; zeros: rho %d, zeta %d' % (key, x, len(gr), len(gz)))
    print('   rho zeros, theta_zeta: measured %.2f%+.2fi ; -psi_rho(x) = %.2f (theorem: Gonek lemma with L_rho\'/L_rho)' % (S_rho_zeta.real, S_rho_zeta.imag, -psi))
    print('   zeta zeros, theta_zeta: measured %.2f%+.2fi ; -psi(x) = %.2f' % (S_zeta_zeta.real, S_zeta_zeta.imag, -sum(math.log(list(fac(n).keys())[0]) for n in range(2, int(x) + 1) if len(fac(n)) == 1)))
    print('   rho zeros, theta_rho: measured %.2f%+.2fi ; law -c T/2pi = %.2f' % (S_rho_rho.real, S_rho_rho.imag, -c * x))
    print('   rho zeros, theta_K: measured %.2f%+.2fi ; simple a=1 formula -c sum rho(1/k) = %.2f ; exact integrals a=1: %.2f%+.2fi ; all a<=12: %.2f%+.2fi'
          % (S_rho_K.real, S_rho_K.imag, simple, tot_a1.real, tot_a1.imag, tot_all.real, tot_all.imag))
    print('      by a: ' + ' '.join('a=%d:%.1f%+.1fi' % (a, v.real, v.imag) for a, v in sorted(by_a.items())))
    print('   zeta zeros, theta_K: measured %.2f%+.2fi' % (S_zeta_K.real, S_zeta_K.imag))
    tot_meas = S_rho_K + S_zeta_K
    print('   all zeros theta_K: measured %.2f%+.2fi ; prediction (rho lines, all a) %.2f%+.2fi ; residual %.2f%+.2fi'
          % (tot_meas.real, tot_meas.imag, tot_all.real, tot_all.imag, (tot_meas - tot_all).real, (tot_meas - tot_all).imag))
    sys.stdout.flush()
