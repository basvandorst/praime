# Order-6 nebentypus newform f of level 13, weight 2 (and its dual f-bar): zeros from sign changes of Z on a grid.
# Coefficients embedded EXACTLY at high precision (t -> exp(i pi/3) in the exact q-expansion); PARI's mfembed/lfunmf
# coefficients are accurate to ~1e-19 only, which ruins lfun for t > ~28 (cancellation e^{pi t/2}).
import os, sys, math, json, time, numpy as np, cypari2
from scipy.special import loggamma
pari = cypari2.Pari(); pari.allocatemem(1400 * 10 ** 6)
OUT = os.path.dirname(os.path.abspath(__file__))+'/'
T = float(sys.argv[1]); step = float(sys.argv[2])
N = 13
PD = int(40 + math.pi * T / 2 / math.log(10)) + 15
NC = int((math.pi * T / 2 + PD * 2.31) / (2 * math.pi / math.sqrt(N)) * 1.3) + 50
pari.set_real_precision(PD)
t0 = time.time()
pari('mf = mfinit([13,2,Mod(4,13)],0)'); pari('F = mfeigenbasis(mf)[1]'); pari('v = mfcoefs(F, %d)' % NC); pari('vt = variable(lift(v[3]))')
pari('z = exp(I*Pi/3)'); pari('W = vector(%d, n, subst(lift(v[n+1]), vt, z))' % NC); pari('Wc = conj(W)')
pari('th(y, A) = sum(n=1, #A, A[n]*exp(-2*Pi*n*y))'); pari('y0 = 0.9/sqrt(13)'); pari('EPS = th(1/(13*y0), W)/(13*y0^2*th(y0, Wc))')
eps = complex(pari('EPS'))
pari('Ld = lfuncreate([W, Wc, [0,1], 2, 13, EPS])')
print('PD %d, NC %d, eps %s, lfuncheckfeq %s (%.0fs)' % (PD, NC, np.round(eps, 12), pari('lfuncheckfeq(Ld)'), time.time() - t0)); sys.stdout.flush()
pari('Li = lfuninit(Ld, [1, 0, %f])' % (T + 2))
print('lfuninit (%.0fs)' % (time.time() - t0)); sys.stdout.flush()
def theta(t):
    t = np.asarray(t, float)
    return 0.5 * t * math.log(N / math.pi ** 2) + np.imag(loggamma((1.0 + 1j * t) / 2)) + np.imag(loggamma((2.0 + 1j * t) / 2))
def Lf(t): return complex(pari('lfun(Li, 1 + (%.15f)*I)' % t))
se = np.sqrt(eps)
def Zf(t): return np.exp(1j * theta(t)) * Lf(t) / se                      # zeros of f at height t
def Zd(t): return np.exp(1j * theta(t)) * np.conj(Lf(-t)) / np.conj(se)   # zeros of f-bar at height t
ts = np.arange(0.01, T + 0.5, step)
res = {}
for name, Zfun, e in (('N13', Zf, eps), ('N13d', Zd, np.conj(eps))):
    zv = np.array([Zfun(t) for t in ts])
    rel = np.abs(zv.imag) / np.abs(zv)
    zr = zv.real
    idx = np.nonzero(np.sign(zr[:-1]) * np.sign(zr[1:]) < 0)[0]
    roots = []
    for i in idx:
        lo, hi, flo = ts[i], ts[i + 1], zr[i]
        for _ in range(22):
            mid = 0.5 * (lo + hi); fm = Zfun(mid).real
            if np.sign(fm) == np.sign(flo): lo, flo = mid, fm
            else: hi = mid
        roots.append(0.5 * (lo + hi))
    g = np.array(roots); g = g[(g > 0) & (g <= T)]
    az = np.abs(zr)
    sus = [round(float(ts[i]), 2) for i in range(1, len(ts) - 1) if az[i] < az[i - 1] and az[i] < az[i + 1] and (i not in idx and i - 1 not in idx) and az[i] < 0.15 * np.median(az)]
    n = len(g); off = np.arange(1, n + 1) - theta(g) / math.pi
    sq = np.sqrt(e)
    print('%s: grid %d pts, max |Im Z|/|Z| %.1e ; zeros %d ; n - theta(T)/pi = %.3f ; mean(n - theta/pi) = %.3f (expected %.3f) ; drift %.3f ; min norm gap %.3f ; |Z| dips w/o sign change: %s (%.0fs)'
          % (name, len(ts), rel.max(), n, n - theta(T) / math.pi, off.mean(), 0.5 - np.angle(sq) / math.pi, off[-n // 6:].mean() - off[:n // 6].mean(),
             (np.diff(g) * np.log(N * (g[:-1] / (2 * math.pi)) ** 2) / (2 * math.pi)).min(), sus, time.time() - t0)); sys.stdout.flush()
    res[name] = dict(eps=[e.real, e.imag], zeros=g.tolist())
json.dump(dict(T=T, PD=PD, res=res), open(OUT + 'r7f_n13_T%d.json' % int(T), 'w'))
c = 1.0215617
for name in ('N13', 'N13d'):
    g = np.array(res[name]['zeros']); e = complex(*res[name]['eps'])
    for TT in (T / 2, T):
        gg = g[g <= TT]; G = np.exp(2j * theta(gg)).sum(); pred = -e * c * TT / (2 * math.pi)
        print('%s T=%.0f: G = %.2f%+.2fi ; law -eps c T/2pi = %.2f%+.2fi ; arg G = %.1f deg ; law %.1f ; no root number 180 ; -conj(eps) c %.1f ; G/pred = %.3f%+.3fi'
              % (name, TT, G.real, G.imag, pred.real, pred.imag, np.degrees(np.angle(G)), np.degrees(np.angle(-e)), np.degrees(np.angle(-np.conj(e))), (G / pred).real, (G / pred).imag))
    pred = -e * c * T / (2 * math.pi)
    for K in (4, 6, 8, 12, 20):
        edges = np.linspace(0, T, K + 1)
        bk = np.array([np.exp(2j * theta(g[(g > edges[k]) & (g <= edges[k + 1])])).sum() for k in range(K)]) / (pred / K)
        print('   K=%2d blocks: G/pred = %.3f%+.3fi  se(Re) %.3f se(Im) %.3f' % (K, bk.mean().real, bk.mean().imag, bk.real.std(ddof=1) / math.sqrt(K), bk.imag.std(ddof=1) / math.sqrt(K)))
    wgt = np.sin(math.pi * g / T) ** 2; Gh = (wgt * np.exp(2j * theta(g))).sum()
    print('   Hann: G_h/(pred/2) = %.3f%+.3fi ; arg G_h %.1f deg' % ((Gh / (pred / 2)).real, (Gh / (pred / 2)).imag, np.degrees(np.angle(Gh))))
print('done (%.0fs)' % (time.time() - t0))
