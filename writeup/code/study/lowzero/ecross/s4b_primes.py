"""Step 4b: the main term from primes alone, no zero used: S(x) = sum_{n <= x} Lambda_E(n) e(-2 sqrt(n/N)) against the
predicted resonance (2 sqrt 2/3) e^{-3 pi i/4} w c_E (x^3/N)^{1/4}, and a smoothed version with the weight
n^{-1/4} sin^2(pi (t_n - 20)/380), t_n = 2 pi sqrt(n/N) in (20, 400) (= the Hann one-point sum up to the factor
-(e^{3 pi i/4}/sqrt 2) N^{-1/4}).  Control: the same sums at a non-resonant frequency e(-2 sqrt(n/N) * 1.13).
Output: s4b_out.txt."""
import math, sys, time, json
import numpy as np, cypari2
C = np.load('curves.npz'); ids = [str(s) for s in C['ids']]
pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)
out = open('s4b_out.txt', 'w')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
P('prime sums S(x) = sum_{n<=x} Lambda_E(n) e(-a sqrt(n/N)); law at a = 2: (2 sqrt2/3) e^{-3 pi i/4} w c_E (x^3/N)^{1/4}')
for cid in ('N431_1_m1_1_m9_m8', 'N563_1_1_1_m15_16', 'N590_1_m1_0_1_5', 'N782_1_m1_1_0_1'):
    i = ids.index(cid); N = int(C['N'][i]); w = int(C['w'][i]); c = float(C['c'][i]); bad = json.loads(str(C['bad'][i]))
    M = int(N * 402 ** 2 / (4 * math.pi ** 2)); a = [int(s.replace('m', '-')) for s in cid.split('_')[1:]]
    pari('E = ellinit(%s); P = primes([2, %d]);' % (str(a), M))
    ap = np.array([int(x) for x in pari('Vecsmall(apply(p -> ellap(E, p), P))')], float); pr = np.array([int(x) for x in pari('Vecsmall(P)')], float)
    lam = ap / np.sqrt(pr); isb = np.isin(pr.astype(np.int64), bad); th = np.arccos(np.clip(lam / 2, -1, 1))
    ns, Ls = [], []; k = 1
    while True:
        m = pr ** k <= M
        if not m.any(): break
        ns.append(pr[m] ** k); Ls.append(np.where(isb[m], lam[m] ** k, 2 * np.cos(k * th[m])) * np.log(pr[m])); k += 1
    n = np.concatenate(ns); Lam = np.concatenate(Ls); o = np.argsort(n); n, Lam = n[o], Lam[o]
    for aa in (2.0, 2.0 * 1.13):
        e = Lam * np.exp(-2j * math.pi * aa * np.sqrt(n / N)); S = np.cumsum(e)
        xs = [M // 8, M // 4, M // 2, M]
        vals = [S[np.searchsorted(n, x, 'right') - 1] for x in xs]
        law = [(2 * math.sqrt(2) / 3) * np.exp(-3j * math.pi / 4) * w * c * (x ** 3 / N) ** 0.25 if aa == 2.0 else 0 for x in xs]
        tn = 2 * math.pi * np.sqrt(n / N); wt = np.where((tn > 20) & (tn < 400), np.sin(math.pi * (tn - 20) / 380) ** 2, 0) * n ** -0.25
        sm = np.sum(e * wt)
        # law of the smoothed sum: sum mu(n) n^{-1/4} wt, mu(n) = (w c/sqrt2) e^{-3 pi i/4} (nN)^{-1/4}  ->  integral
        xx = np.linspace(1, M, 400001); wtx = np.where((2 * math.pi * np.sqrt(xx / N) > 20) & (2 * math.pi * np.sqrt(xx / N) < 400), np.sin(math.pi * (2 * math.pi * np.sqrt(xx / N) - 20) / 380) ** 2, 0)
        lsm = (w * c / math.sqrt(2)) * np.exp(-3j * math.pi / 4) * N ** -0.25 * np.trapezoid(xx ** -0.5 * wtx, xx) if aa == 2.0 else 0
        P('  %s N=%d w=%+d c=%.3f a=%.2f: S(x) at x = M/8..M: %s | law %s | smoothed %.1f%+.1fi, law %.1f%+.1fi; random-sign scale sqrt(sum Lam^2 wt^2) %.1f' % (
            cid, N, w, c, aa, ' '.join('%.0f%+.0fi' % (v.real, v.imag) for v in vals), ' '.join('%.0f%+.0fi' % (complex(v).real, complex(v).imag) for v in law),
            sm.real, sm.imag, complex(lsm).real, complex(lsm).imag, math.sqrt(np.sum((Lam * wt) ** 2))))
out.close()
