"""Numerical tables for the proof note (5 curves): zeros vs primes (sharp, Hann, Gaussian-edged), twists, Voronoi over
integers, sharp Voronoi partial sums, prime-side law.  Output: num_tables_out.txt and num_tables.json."""
import os, sys, math, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, cypari2
import gl2lib as G
from exact_check import In_all

OUT = os.path.dirname(os.path.abspath(__file__)) + '/'
CURVES_ALL = ['N307_0_0_1_m8_m9', 'N335_0_0_1_m2_2', 'N373_0_1_1_m2_m2', 'N611_0_0_1_m1_1', 'N870_1_0_0_m5_m3']
CURVES = sys.argv[1:] or CURVES_ALL
WINS = [(20, 400), (20, 140), (140, 260), (260, 400), (60, 230), (230, 380)]
pred = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'pred_fresh.npz'))
cmap = dict(zip([str(s) for s in pred['ids']], pred['c']))
f = open(OUT + 'num_tables_out.txt', 'a')
def out(s=''):
    print(s); f.write(s + '\n'); f.flush()
def cs(z):
    return '%8.3f%+8.3fi' % (z.real, z.imag)
res = {}
pari = cypari2.Pari(); pari.allocatemem(450 * 10 ** 6)
for cid in CURVES:
    C = G.load_curves()[cid]; N, w = C['N'], C['w']
    gam = G.load_zeros(cid)
    t0 = time.time()
    X2 = int(2 * N * 402.0 ** 2 / (4 * math.pi ** 2)) + 10        # r = 1/2 needs 2x
    n, Lam, ap50 = G.prime_power_data(pari, C['ainv'], C['bad'], X2)
    X = X2 // 2
    # independent c_E: direct partial product prod_{good p <= P} (1 - lambda(p^2)/p) prod_{p|N} (1 + 1/p), P = X2
    pr_mask = np.array([True] * len(n))
    pari('P__ = primes([2, %d]);' % X2)
    pr = np.array(pari('Vecsmall(P__)'), float)
    ap = np.array(pari('Vecsmall(apply(p -> ellap(E__, p), P__))'), float)
    good = ~np.isin(pr.astype(np.int64), np.array(C['bad']))
    lam2 = ap[good] ** 2 / pr[good] - 1
    c_direct = math.exp(np.sum(np.log1p(-lam2 / pr[good])) + sum(math.log1p(1 / q) for q in C['bad']))
    c = cmap[cid]
    out('=== %s  N=%d  w=%+d  c_E (frozen B/L(1,Sym2)) = %.4f, direct partial product to %d: %.4f  (%d prime powers, %.0fs)' % (
        cid, N, w, c, X2, c_direct, len(n), time.time() - t0))
    R = res[cid] = dict(N=N, w=w, c=c, c_direct=c_direct, rows=[])
    # ---- Hann and sharp windows
    out('  window      | weight | zeros Z            | primes P1 (exact th) | P1+2nd order       | P2 (e(-2sqrt(n/N)))  | |Z-P1| |Z-P1b| |Z-P2| | law -w c intH/2pi')
    for kind in ('hann', 'box'):
        for (a, b) in WINS:
            H = G.Weight(kind, a, b, 0.0)
            Z = G.zero_side(N, gam, H)
            P1 = G.prime_side(N, n, Lam, H, mode='exact')
            P1b = G.prime_side(N, n, Lam, H, mode='exact', second_order=True) if kind == 'hann' else P1
            P2 = G.prime_side(N, n, Lam, H, mode='simple')
            law = -w * c * H.integral() / G.TWO_PI
            out('  [%3d,%3d]   | %-6s | %s | %s | %s | %s | %6.3f %6.3f %6.3f | %8.3f' % (
                a, b, 'sharp' if kind == 'box' else 'Hann', cs(Z), cs(P1), cs(P1b), cs(P2), abs(Z - P1), abs(Z - P1b), abs(Z - P2), law))
            R['rows'].append(dict(win=[a, b], weight='sharp' if kind == 'box' else 'hann', Z=[Z.real, Z.imag], P1=[P1.real, P1.imag],
                                  P1b=[P1b.real, P1b.imag], P2=[P2.real, P2.imag], law=law))
    # ---- Gaussian-edged boxes: zeros, exact quadrature (Prop 1), stationary phase
    for (a, b, D) in [(70, 180, 10.0), (180, 330, 10.0), (125, 275, 20.0)]:
        H = G.Weight('box', a, b, D)
        lo, hi = H.support()
        Z = G.zero_side(N, gam, H)
        xlo = N * lo ** 2 / (4 * math.pi ** 2) * math.exp(-1.0); xhi = N * hi ** 2 / (4 * math.pi ** 2) * math.exp(1.0)
        m = (n >= xlo) & (n <= xhi)
        Pex = -(1 / G.TWO_PI) * np.sum(Lam[m] * n[m] ** -0.5 * In_all(N, H, np.log(n[m]), h=0.05))
        P1 = G.prime_side(N, n, Lam, H, mode='exact'); P1b = G.prime_side(N, n, Lam, H, mode='exact', second_order=True)
        law = -w * c * H.integral() / G.TWO_PI
        out('  box [%3d,%3d] Delta=%2d: zeros %s | primes exact I_n %s (|Z-P| %.1e) | stationary %s (%.3f), 2nd order %s (%.3f) | law %.3f' % (
            a, b, D, cs(Z), cs(Pex), abs(Z - Pex), cs(P1), abs(Z - P1), cs(P1b), abs(Z - P1b), law))
        R['rows'].append(dict(win=[a, b], weight='box%g' % D, Z=[Z.real, Z.imag], Pex=[Pex.real, Pex.imag], P1=[P1.real, P1.imag], P1b=[P1b.real, P1b.imag], law=law))
    # ---- twists r^{-i gamma} (Hann [20,400])
    H = G.Weight('hann', 20, 400)
    for r in (2.0, 0.5, 3.0, 1.5):
        Z = G.zero_side(N, gam, H, r=r)
        P1 = G.prime_side(N, n, Lam, H, r=r, mode='exact'); P2 = G.prime_side(N, n, Lam, H, r=r, mode='simple')
        out('  twist r=%.2f Hann[20,400]: zeros %s | primes P1 %s | P2 %s | |Z-P1| %.3f' % (r, cs(Z), cs(P1), cs(P2), abs(Z - P1)))
        R['rows'].append(dict(win=[20, 400], weight='hann', r=r, Z=[Z.real, Z.imag], P1=[P1.real, P1.imag], P2=[P2.real, P2.imag]))
    # ---- Voronoi over all integers
    t0 = time.time()
    nn, lam = G.lambda_all(pari, C['ainv'], X)
    for (a, b) in [(20, 400), (20, 140), (140, 260)]:
        H = G.Weight('hann', a, b)
        V1 = G.prime_side(N, nn, lam, H, mode='exact')             # -(1/2pi) sum lambda(n) n^-1/2 [stat. phase I_n]
        V2 = G.prime_side(N, nn, lam, H, mode='simple')
        V1b = G.prime_side(N, nn, lam, H, mode='exact', second_order=True)
        # exact dual side: -(w/2pi) sum_m lambda(m) m^-1/2 Hhat(-ln m), Hann: Hhat(u) = int_a^b sin^2(pi(t-a)/L) e^{itu} dt
        L_ = b - a; k = 2 * math.pi / L_
        mm = np.arange(1, 200, dtype=float); u = np.log(mm)
        def Fbox(u):
            return np.where(np.abs(u) < 1e-14, L_, (np.exp(1j * b * u) - np.exp(1j * a * u)) / (1j * np.where(np.abs(u) < 1e-14, 1, u)))
        Hh = 0.5 * Fbox(u) - 0.25 * (np.exp(-1j * k * a) * Fbox(u + k) + np.exp(1j * k * a) * Fbox(u - k))
        dual = -(w / G.TWO_PI) * np.sum(lam[:199] * mm ** -0.5 * Hh)
        out('  Voronoi Hann[%3d,%3d]: -(1/2pi) sum lambda(n) n^-1/2 I_n by stat. phase: %s (2nd order %s), simple %s | exact dual -(w/2pi) sum_m lambda(m) m^-1/2 Hhat = %s (m=1: %.3f) | diff %.3f %.3f' % (
            a, b, cs(V1), cs(V1b), cs(V2), cs(dual), -w * H.integral() / G.TWO_PI, abs(V1 - dual), abs(V1b - dual)))
        R['rows'].append(dict(win=[a, b], weight='hann', voronoi=True, V1=[V1.real, V1.imag], V1b=[V1b.real, V1b.imag], V2=[V2.real, V2.imag], dual=[dual.real, dual.imag]))
    # sharp Voronoi partial sums, and the prime analogue
    e = lam * np.exp(-4j * math.pi * np.sqrt(nn / N)); S = np.cumsum(e)
    ep = Lam * np.exp(-4j * math.pi * np.sqrt(n / N)); Sp = np.cumsum(ep)
    for xx in (X // 64, X // 16, X // 4, X):
        sv = S[xx - 1]; main = (2 * math.sqrt(2) / 3) * np.exp(-3j * math.pi / 4) * w * (xx ** 3 / N) ** 0.25
        j = np.searchsorted(n, xx, 'right') - 1
        out('  sharp: x=%9d  sum_{n<=x} lambda(n) e(-2sqrt(n/N)) = %s  vs (2sqrt2/3)e^{-3pi i/4} w (x^3/N)^{1/4} = %s  (ratio %.3f, |diff|/x^{1/2} %.2f) | primes: sum Lambda_E e(..) = %s vs c_E x main %s' % (
            xx, cs(sv), cs(main), (sv / main).real, abs(sv - main) / math.sqrt(xx), cs(Sp[j]), cs(c * main)))
        R['rows'].append(dict(sharp_x=int(xx), S=[sv.real, sv.imag], main=[main.real, main.imag], Sp=[Sp[j].real, Sp[j].imag]))
    out('  (Voronoi part %.0fs)' % (time.time() - t0))
    del nn, lam, e, S
json.dump(res, open(OUT + 'num_tables_%s.json' % CURVES[0], 'w'), indent=1)
