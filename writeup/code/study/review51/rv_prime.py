# The Voronoi resonance over primes, without zeros: S_k(x) = sum_{n<=x} Lambda_E(n) (1-n/x)^k e(-2 sqrt(n/N))
# for x up to 1e7 (a_p by PARI ellap), against the law (2 sqrt2/3) e^{-3 pi i/4} w c (x^3/N)^{1/4} * (3/4) B(3/4,k+1);
# and the band-limited stationary-phase sum R (heights 20-400, sharp and Hann) against the zero sum G.
import os, sys, time
import numpy as np
import cypari2
from scipy.special import loggamma, beta
SC = os.path.dirname(os.path.abspath(__file__))+'/'
CU = np.load(SC + 'rv_curves.npz'); LW = np.load(SC + 'rv_law.npz')
ids, Nall, wall, call = CU['ids'], CU['N'], CU['w'], LW['c']
sel = [int(s) for s in sys.argv[1].split(',')]
XMAX = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10**7
pari = cypari2.Pari(); pari.allocatemem(600 * 10**6)
primes = np.array([int(x) for x in pari('primes(primepi(%d))' % XMAX)], dtype=np.int64)
Z = np.load(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/data/study/r144/fresh_zeros.npz')
out = {}
for k in sel:
    t0 = time.time()
    s = ids[k]; N = int(Nall[k]); w = int(wall[k]); c = float(call[k])
    a = [int(x.replace('m', '-')) for x in s.split('_')[1:]]
    pari('global(EE); EE = ellinit(%s)' % str(a))
    ap = np.array([int(x) for x in pari('my(v=vectorsmall(%d),i=0);forprime(p=2,%d,i++;v[i]=ellap(EE,p));v' % (len(primes), XMAX))], float)
    tap = time.time() - t0
    pf = primes.astype(float); lam = ap / np.sqrt(pf); bad = (N % primes == 0)
    th = np.arccos(np.clip(lam / 2, -1, 1))
    ns, Ls = [], []
    for kk in range(1, 25):
        m = pf ** kk <= XMAX
        if not m.any():
            break
        if kk == 1:
            coef = lam[m]
        else:
            coef = np.where(bad[m], lam[m] ** kk, 2 * np.cos(kk * th[m]))
        ns.append(pf[m] ** kk); Ls.append(coef * np.log(pf[m]))
    n = np.concatenate(ns); L = np.concatenate(Ls); o = np.argsort(n); n = n[o]; L = L[o]
    ph = np.exp(-4j * np.pi * np.sqrt(n / N))
    # raw and Riesz sums on an x grid
    xs = np.unique(np.round(np.geomspace(1e4, XMAX, 40)))
    res = []
    cs = np.cumsum(L * ph); csn = np.cumsum(L * ph * n); csn2 = np.cumsum(L * ph * n * n)
    for x in xs:
        i = np.searchsorted(n, x, side='right') - 1
        S0 = cs[i]; S1 = S0 - csn[i] / x; S2 = S0 - 2 * csn[i] / x + csn2[i] / x ** 2
        law0 = (2 * np.sqrt(2) / 3) * np.exp(-3j * np.pi / 4) * w * c * (x ** 3 / N) ** 0.25
        laws = [law0 * 0.75 * beta(0.75, kk + 1) for kk in range(3)]
        res.append([x, 2 * np.pi * np.sqrt(x / N)] + [ (S * np.conj(Lw)).real / abs(Lw) ** 2 for S, Lw in zip((S0, S1, S2), laws)] +
                   [abs(S0 - laws[0]) / abs(laws[0])])
    res = np.array(res)
    # band sum with the exact stationary point t_n: ln(N t^2/4pi^2) + 1/(6t^2) + 1/(60t^4) = ln n
    t = 2 * np.pi * np.sqrt(n / N)
    for it in range(6):
        g = np.log(N * t * t / (4 * np.pi ** 2)) + 1 / (6 * t * t) + 1 / (60 * t ** 4) - np.log(n)
        dg = 2 / t - 1 / (3 * t ** 3) - 1 / (15 * t ** 5)
        t = t - g / dg
    th2 = t * np.log(N) - 2 * t * np.log(2 * np.pi) + 2 * loggamma(1 + 1j * t).imag      # 2 theta(t_n)
    f = th2 - t * np.log(n)
    fpp = 2 / t - 1 / (3 * t ** 3)
    term = -(1 / (2 * np.pi)) * L / np.sqrt(n) * np.sqrt(2 * np.pi / fpp) * np.exp(1j * (f + np.pi / 4))
    inb = (t > 20) & (t <= 400)
    Rs = term[inb].sum()
    hann = np.sin(np.pi * (t - 20) / 380) ** 2
    Rh = (term * hann)[inb].sum()
    z = Z[s]; g = z[(z > 20) & (z <= 400)]
    e2 = np.exp(1j * (g * np.log(N) - 2 * g * np.log(2 * np.pi) + 2 * loggamma(1 + 1j * g).imag))
    Gs = e2.sum(); Gh = (e2 * np.sin(np.pi * (g - 20) / 380) ** 2).sum()
    lawG = -w * c * 380 / (2 * np.pi)
    nmax_needed = N * 400 ** 2 / (4 * np.pi ** 2)
    print('%s N=%d w=%+d c=%.3f | ellap to %.0e: %.1fs | n needed for t=400: %.2e' % (s, N, w, c, XMAX, tap, nmax_needed))
    print('   zeros G sharp %8.2f%+8.2fi  primes R sharp %8.2f%+8.2fi | Hann: zeros %8.2f%+8.2fi primes %8.2f%+8.2fi | law sharp %.2f, Hann %.2f' % (
        Gs.real, Gs.imag, Rs.real, Rs.imag, Gh.real, Gh.imag, Rh.real, Rh.imag, lawG, lawG / 2))
    for row in res[::4].tolist() + [res[-1].tolist()]:
        print('   x %.3e (t_x %6.1f): Re(S/law) k=0 %.3f k=1 %.3f k=2 %.3f |S0-law|/|law| %.3f' % tuple(row))
    out[s] = dict(res=res, Rs=Rs, Rh=Rh, Gs=Gs, Gh=Gh, lawG=lawG, N=N, w=w, c=c)
    np.save(SC + 'rv_prime_%s.npy' % s, out[s], allow_pickle=True)
