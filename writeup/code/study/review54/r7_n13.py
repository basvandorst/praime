# Dropped class: the weight-2 newform of level 13 with a character of order 6 (non-self-dual, complex root number).
# Own route: explicit complex coefficient vector -> lfuncreate (dual = conjugate), root number by lfunrootres,
# Z(t) = eps^{-1/2} e^{i theta(t)} L(1/2+it) on a grid, sign changes, refinement; then the one-point sum and its phase.
import os, sys, math, json, time, numpy as np, cypari2
from scipy.special import loggamma
pari = cypari2.Pari(); pari.allocatemem(1200 * 10 ** 6)
OUT = os.path.dirname(os.path.abspath(__file__))+'/'
T = float(sys.argv[1]) if len(sys.argv) > 1 else 120.0
charn = int(sys.argv[2]) if len(sys.argv) > 2 else 4
step = float(sys.argv[3]) if len(sys.argv) > 3 else 0.04
NC = 30000
t0 = time.time()
mf = pari('mfinit([13,2,Mod(%d,13)],0)' % charn)
Bs = pari.mfeigenbasis(mf)
print('char Mod(%d,13): order %s, newforms %d' % (charn, pari('znorder(Mod(%d,13))' % charn), len(Bs)))
F = Bs[0]
co = pari.mfcoefs(F, NC)
emb = pari.mfembed(F, co)
print('mfembed type', pari.type(emb), 'len', len(emb), 'type of first', pari.type(emb[0]))
embs = list(emb) if pari.type(emb[0]) == 't_VEC' else [emb]

def theta(t):
    t = np.asarray(t, float)
    return 0.5 * t * math.log(13 / math.pi ** 2) + np.imag(loggamma((1.0 + 1j * t) / 2)) + np.imag(loggamma((2.0 + 1j * t) / 2))

out = []
for ei, E in enumerate(embs):
    a = np.array([complex(E[n]) for n in range(0, NC + 1)])
    # mfcoefs returns [a0, a1, ..., aNC]
    a = a[1:]
    a = a / a[0]
    chi2 = (a[1] ** 2 - a[3]) / 2
    print('embedding %d: a(2)=%s a(3)=%s a(13)=%s | chi(2) from a(4): %s (arg/60deg = %.4f) | |a(13)|^2/13 = %.4f'
          % (ei, np.round(a[1], 6), np.round(a[2], 6), np.round(a[12], 6), np.round(chi2, 6), np.angle(chi2) / (math.pi / 3), abs(a[12]) ** 2 / 13))
    # constant c = prod_{p != 13} (1 - (|lam_p|^2 - 1)/p), direct, log-Cesaro over [P/10, P], bias nu = 1
    s = np.ones(NC + 1, bool); s[:2] = False
    for i in range(2, int(NC ** .5) + 1):
        if s[i]: s[i * i::i] = False
    pr = np.nonzero(s)[0]
    lam2 = np.array([abs(a[p - 1]) ** 2 / p for p in pr])
    lt = np.where(pr == 13, 0.0, np.log(1 - (lam2 - 1) / pr))
    cum = np.cumsum(lt); m = pr >= NC / 10; w = 1.0 / pr[m]
    c_dir = math.exp(((cum[m] + 1 / (np.sqrt(pr[m]) * np.log(pr[m]))) * w).sum() / w.sum())
    print('   c direct (P = %d, log-Cesaro + bias) = %.4f ; mean |lam_p|^2 = %.4f' % (NC, c_dir, lam2.mean()))
    vec = pari([pari.complex(float(x.real), float(x.imag)) for x in a])
    pari.set_real_precision(38)
    Ld = pari.lfuncreate([vec, 1, [0, 1], 2, 13, 1])
    rr = pari.lfunrootres(Ld)
    eps = complex(rr[2])
    print('   lfunrootres eps = %s (|eps| = %.10f)' % (np.round(eps, 10), abs(eps)))
    Ld = pari.lfuncreate([vec, 1, [0, 1], 2, 13, pari.complex(eps.real, eps.imag)])
    print('   lfuncheckfeq = %s' % pari.lfuncheckfeq(Ld))
    Li = pari.lfuninit(Ld, [1, T + 3])          # arithmetic s = 1 + it  <-> analytic 1/2 + it
    print('   lfuncheckfeq (init) = %s ; (%.0fs)' % (pari.lfuncheckfeq(Li), time.time() - t0)); sys.stdout.flush()
    seps = np.sqrt(eps)
    def Z(t):
        v = complex(pari.lfun(Li, pari.complex(1.0, float(t))))
        return np.exp(1j * theta(t)) * v / seps
    ts = np.arange(0.02, T + 1.0, step)
    zs = np.array([Z(t) for t in ts])
    relim = np.abs(zs.imag) / np.maximum(np.abs(zs), 1e-300)
    print('   grid %d points: max |Im Z|/|Z| = %.2e (median %.2e) ; (%.0fs)' % (len(ts), relim.max(), np.median(relim), time.time() - t0)); sys.stdout.flush()
    zr = zs.real
    # check also the opposite sign convention for the square root (Z real either way; sign irrelevant)
    idx = np.nonzero(np.sign(zr[:-1]) * np.sign(zr[1:]) < 0)[0]
    roots = []
    for i in idx:
        lo, hi = ts[i], ts[i + 1]; flo = zr[i]
        for _ in range(30):
            mid = 0.5 * (lo + hi); fm = Z(mid).real
            if np.sign(fm) == np.sign(flo): lo, flo = mid, fm
            else: hi = mid
        roots.append(0.5 * (lo + hi))
    g = np.array(roots); g = g[g <= T]
    # local minima of |Z| without sign change (possible missed close pairs)
    az = np.abs(zr)
    sus = [ts[i] for i in range(1, len(ts) - 1) if az[i] < az[i - 1] and az[i] < az[i + 1] and az[i] < 0.05 * np.median(az) and not (i in idx or i - 1 in idx)]
    n = len(g)
    off = np.arange(1, n + 1) - theta(g) / math.pi
    print('   zeros found %d in (0, %.0f]; n - theta(T)/pi = %.3f ; mean(n - theta(g_n)/pi) = %.3f (expected 1/2 - arg(eps^1/2)/pi = %.3f); drift %.3f ; suspicious minima %s'
          % (n, T, n - theta(T) / math.pi, off.mean(), 0.5 - np.angle(seps) / math.pi, off[-n // 5:].mean() - off[:n // 5].mean(), np.round(sus, 2)))
    gaps = np.diff(g); dens = (np.log(13 * (g[:-1] / (2 * math.pi)) ** 2)) / (2 * math.pi)
    print('   min normalised gap %.3f' % (gaps * dens).min())
    c_froz = 1.0215617
    G = np.exp(2j * theta(g)).sum()
    pred = -eps * c_froz * T / (2 * math.pi)
    print('   G = %.2f%+.2fi ; law -eps c T/2pi = %.2f%+.2fi ; arg G = %.1f deg ; law %.1f ; no root number 180 ; conj(eps) %.1f ; Re(G/pred) = %.3f , Im(G/pred) = %.3f'
          % (G.real, G.imag, pred.real, pred.imag, np.degrees(np.angle(G)), np.degrees(np.angle(-eps)), np.degrees(np.angle(-np.conj(eps))), (G / pred).real, (G / pred).imag))
    # Hann
    wgt = np.sin(math.pi * g / T) ** 2
    Gh = (wgt * np.exp(2j * theta(g))).sum()
    print('   Hann: G_h = %.2f%+.2fi ; law/2 = %.2f%+.2fi ; arg %.1f deg' % (Gh.real, Gh.imag, pred.real / 2, pred.imag / 2, np.degrees(np.angle(Gh))))
    # blocks
    K = 6; edges = np.linspace(0, T, K + 1)
    Gk = np.array([np.exp(2j * theta(g[(g > edges[k]) & (g <= edges[k + 1])])).sum() for k in range(K)])
    bk = Gk / (pred / K)
    print('   block ratios G_k/pred_k: %s ; mean %.3f%+.3fi, se %.3f' % (np.round(bk, 2), bk.mean().real, bk.mean().imag, bk.std(ddof=1) / math.sqrt(K)))
    out.append(dict(emb=ei, eps=[eps.real, eps.imag], c_direct=c_dir, zeros=g.tolist(), G=[G.real, G.imag], Gh=[Gh.real, Gh.imag], T=T))
    sys.stdout.flush()
json.dump(out, open(OUT + 'r7_n13_char%d_T%d.json' % (charn, int(T)), 'w'))
print('done (%.0fs)' % (time.time() - t0))
