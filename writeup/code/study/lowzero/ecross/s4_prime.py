"""Step 4 (development curves, smallest conductors): the one-point sum is a prime sum.  By the explicit formula and
stationary phase, G_w = sum_gamma e^{2 i theta(gamma)} w(gamma) = -(1/2pi) sum_n Lambda_E(n) n^{-1/2} hhat(ln n) + small,
hhat(ln n) at the stationary point 2 theta'(t_n) = ln n, t_n ~ 2 pi sqrt(n/N), phase -2 t_n + pi/2 + ...: the GL(2)
Voronoi/Bessel kernel e(-2 sqrt(n/N)).  Hann taper on [20, 400].  Lambda_E(p^k) = (alpha^k + beta^k) ln p.
Output: s4_out.txt."""
import os, math, sys, time
import numpy as np, cypari2
sys.path.insert(0, '.')
import elib as L
T1, T2 = 20.0, 400.0; T = T2 - T1
C = np.load('curves.npz'); ids = [str(s) for s in C['ids']]
DEVZ = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'r148', 'dev_zeros.npz'))
pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)
taper = lambda t: np.where((t > T1) & (t < T2), np.sin(math.pi * (t - T1) / T) ** 2, 0.0)
out = open('s4_out.txt', 'w')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')
P('one-point sum (Hann) from the zeros against the prime sum (stationary phase, exact theta), development curves')
for cid in ('N431_1_m1_1_m9_m8', 'N563_1_1_1_m15_16', 'N590_1_m1_0_1_5'):
    i = ids.index(cid); N = int(C['N'][i]); w = int(C['w'][i]); c = float(C['c'][i])
    import json; bad = json.loads(str(C['bad'][i]))
    M = int(N * (T2 + 2) ** 2 / (4 * math.pi ** 2)) + 10
    t0 = time.time()
    ainv = [int(x) for x in pari('ellinit(ellfromeqn(0))') ] if False else None
    # curve from its id: a-invariants in the id string
    a = [int(s.replace('m', '-')) for s in cid.split('_')[1:]]
    pari('E = ellinit(%s); P = primes([2, %d]);' % (str(a), M))
    ap = np.array([int(x) for x in pari('Vecsmall(apply(p -> ellap(E, p), P))')], float)
    pr = np.array([int(x) for x in pari('Vecsmall(P)')], float)
    lam = ap / np.sqrt(pr); isbad = np.isin(pr.astype(np.int64), bad)
    th = np.arccos(np.clip(lam / 2, -1, 1))
    ns, Ls, ks = [], [], []
    k = 1
    while True:
        m = pr ** k <= M
        if not m.any(): break
        Lam = np.where(isbad[m], lam[m] ** k, 2 * np.cos(k * th[m])) * np.log(pr[m])
        ns.append(pr[m] ** k); Ls.append(Lam); ks.append(np.full(m.sum(), k)); k += 1
    n = np.concatenate(ns); Lam = np.concatenate(Ls); kk = np.concatenate(ks)
    R = L.stationary_R(N, taper, n, Lam, T1, T2)
    R1 = L.stationary_R(N, taper, n[kk == 1], Lam[kk == 1], T1, T2)
    g = np.sort(DEVZ[cid]); b = g[(g > T1) & (g < T2)]
    Gz = np.sum(np.exp(2j * L.theta_E(N, b)) * taper(b))
    law = -w * c * T / (2 * math.pi) / 2
    P('  %s N=%d w=%+d c_E=%.3f: n <= %d (%d primes, %.0fs)' % (cid, N, w, c, M, len(pr), time.time() - t0))
    P('     zeros G = %.2f %+.2fi | primes R = %.2f %+.2fi (primes only %.2f %+.2fi) | law (Hann) %.2f | |G - R| = %.2f' % (
        Gz.real, Gz.imag, R.real, R.imag, R1.real, R1.imag, law, abs(Gz - R)))
    # the main term is carried by the primes: cumulative prime sum over t-bands
    edges = np.linspace(T1, T2, 6)
    tn = 2 * math.pi * np.sqrt(n / N)
    parts = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (tn >= lo) & (tn < hi)
        parts.append(L.stationary_R(N, taper, n[m], Lam[m], T1, T2).real)
    zp = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        bb = b[(b >= lo) & (b < hi)]; zp.append(np.sum(np.exp(2j * L.theta_E(N, bb)) * taper(bb)).real)
    lawp = [-w * c / (2 * math.pi) * (lambda a_, b_: (b_ - a_) / 2 - T / (4 * math.pi) * (math.sin(2 * math.pi * (b_ - T1) / T) - math.sin(2 * math.pi * (a_ - T1) / T)))(lo, hi) for lo, hi in zip(edges[:-1], edges[1:])]
    P('     by height fifth: primes ' + ' '.join('%7.2f' % x for x in parts) + ' | zeros ' + ' '.join('%7.2f' % x for x in zp) + ' | law ' + ' '.join('%7.2f' % x for x in lawp))
out.close()
