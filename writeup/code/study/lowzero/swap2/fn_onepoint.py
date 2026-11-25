"""One-point Gram sum on the USED set Fn: the mean law -kappa (T2 - T1)/2pi against a per-function resonance predictor
(stationary phase of int e^{2 i theta(t)} (L'/L)(1/2 + it) w(t) dt: the orbit n resonates at t_n = 2 pi n / q with amplitude
1/sqrt(q) and phase i^a e(-n/q)):
   R_w(chi) = -(i^a / sqrt q) sum_n chi(n) Lambda(n) e(-n/q) w(2 pi n / q),
   G_w(chi) = sum_gamma e^{2 i theta(gamma)} w(gamma),   w = sharp window (100, 600] or a Hann taper on [100, 600].
Mean of R_w: -eps kappa int w dt / 2pi (PNT in progressions: sum_b chi(b) e(-b/q) = chi(-1) tau(chi))."""
import sys, math, os, time
import numpy as np
from scipy.special import loggamma
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hlib import *

T1, T2 = 100.0, 600.0
t0 = time.time()


def theta(t, q, a):
    t = np.asarray(t, float)
    return t / 2 * math.log(q / math.pi) + np.imag(loggamma(0.25 + a / 2 + 0.5j * t))


def hann(t):
    x = (np.asarray(t) - T1) / (T2 - T1)
    return np.where((x > 0) & (x <= 1), np.sin(math.pi * x) ** 2, 0.0)


def sharp(t):
    t = np.asarray(t); return ((t > T1) & (t <= T2)).astype(float)


z = np.load(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'data', 'riemann', 'lzeros')) + '/r138_Fn.npz')
Ds = sorted(int(k[1:]) for k in z.files if int(k[1:]) != 1153)
pari = pari_instance()
nmax = int(1300 * T2 / (2 * math.pi)) + 10
pr = primes_upto(nmax)
# prime powers <= nmax
pp, lp, base = [], [], []
for p in pr:
    pk = p; k = 1
    while pk <= nmax:
        pp.append(pk); lp.append(math.log(p)); base.append((p, k)); pk *= p; k += 1
pp = np.array(pp); lp = np.array(lp); kk = np.array([b[1] for b in base]); pidx = np.searchsorted(pr, [b[0] for b in base])
rows = []
for D in Ds:
    q = abs(D); a = 0 if D > 0 else 1
    chi_p = np.array([int(x) for x in pari('apply(p->kronecker(%d,p), primes([2,%d]))' % (D, nmax))], float)
    chi_n = chi_p[pidx] ** kk
    g = np.sort(z['D%d' % D]); ph = np.exp(2j * theta(g, q, a))
    out = [D, phi_ratio(q)]
    for w in (sharp, hann):
        G = (ph * w(g)).sum()
        R = -(1j ** a / math.sqrt(q)) * np.sum(chi_n * lp * np.exp(-2j * math.pi * (pp % q) / q) * w(2 * math.pi * pp / q))
        out += [G, R]
    rows.append(out)
rows = np.array(rows, dtype=object)
kap = rows[:, 1].astype(float)
for k, (nm, wint) in enumerate((('sharp', T2 - T1), ('hann', (T2 - T1) / 2))):
    G = rows[:, 2 + 2 * k].astype(complex); R = rows[:, 3 + 2 * k].astype(complex); M = -kap * wint / (2 * math.pi)
    print('%s window: mean law slope Re G/M %.4f; resonance R: slope Re R/M %.4f (Im R mean %.2f);' % (
        nm, (G.real * M).sum() / (M * M).sum(), (R.real * M).sum() / (M * M).sum(), R.imag.mean()))
    print('   sd over L-functions of (G - M): Re %.2f Im %.2f;  of (G - R): Re %.2f Im %.2f;  corr(Re G - M, Re R - M) %.3f, corr Im %.3f' % (
        (G - M).real.std(), (G - M).imag.std(), (G - R).real.std(), (G - R).imag.std(),
        np.corrcoef((G - M).real, (R - M).real)[0, 1], np.corrcoef(G.imag, R.imag)[0, 1]))
    b = np.polyfit((R - M).real, (G - M).real, 1); print('   regression of (G - M) on (R - M), Re: slope %.3f icpt %.2f' % (b[0], b[1]))
print('%.0fs' % (time.time() - t0))
