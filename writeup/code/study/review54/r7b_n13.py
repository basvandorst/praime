# Order-6 nebentypus newform of level 13, weight 2: zeros by sign changes of Z(t) from PARI lfun on lfuninit (lfunmf Ldata),
# values checked against my own incomplete-gamma approximate functional equation (mpmath), root number from the theta relation.
import os, sys, math, json, time, numpy as np, cypari2, mpmath
from scipy.special import loggamma
pari = cypari2.Pari(); pari.allocatemem(1400 * 10 ** 6)
OUT = os.path.dirname(os.path.abspath(__file__))+'/'
T = float(sys.argv[1]); charn = int(sys.argv[2]); step = float(sys.argv[3])
PD = int(40 + math.pi * T / 2 / math.log(10)) + 15
pari.set_real_precision(PD)
print('PARI realprecision', PD)
t0 = time.time()
mf = pari('mfinit([13,2,Mod(%d,13)],0)' % charn); F = pari.mfeigenbasis(mf)[0]
co = pari.mfcoefs(F, 2000); E = pari.mfembed(F, co)
a = np.array([complex(E[n]) for n in range(len(E))])[1:]
A_mp = [mpmath.mpc(*[mpmath.mpf(str(pari.real(E[n])).replace(' E','e')), mpmath.mpf(str(pari.imag(E[n])).replace(' E','e'))]) for n in range(1, min(len(E), 900))]
EPS_mp = None
N = 13
def fth(y, conj=False):
    n = np.arange(1, len(a) + 1); aa = np.conj(a) if conj else a
    return (aa * np.exp(-2 * np.pi * n * y)).sum()
y = 0.9 / math.sqrt(N)
eps = complex(fth(1 / (N * y)) / (N * y * y * fth(y, True)))     # eps = -eta
mpmath.mp.dps = PD
yy = mpmath.mpf(9) / 10 / mpmath.sqrt(N)
th1 = mpmath.fsum(A_mp[k] * mpmath.exp(-2 * mpmath.pi * (k + 1) / (N * yy)) for k in range(len(A_mp)))
th2 = mpmath.fsum(mpmath.conj(A_mp[k]) * mpmath.exp(-2 * mpmath.pi * (k + 1) * yy) for k in range(len(A_mp)))
EPS_mp = th1 / (N * yy * yy * th2)
print('char %d: a(2) = %s, chi(2) = %s, own eps = %s' % (charn, np.round(a[1], 6), np.round((a[1] ** 2 - a[3]) / 2, 6), np.round(eps, 10)))
def theta(t):
    t = np.asarray(t, float)
    return 0.5 * t * math.log(N / math.pi ** 2) + np.imag(loggamma((1.0 + 1j * t) / 2)) + np.imag(loggamma((2.0 + 1j * t) / 2))
def L_own(t):
    dps = int(30 + math.pi * t / 2 / math.log(10)) + 10
    with mpmath.workdps(dps):
        s = mpmath.mpc(1, t); x0 = 2 * mpmath.pi / mpmath.sqrt(N)
        tot = mpmath.mpc(0)
        nmax = int((math.pi * t / 2 + dps * 2.31) / float(x0)) + 5
        for n in range(1, min(nmax, len(a)) + 1):
            an = A_mp[n - 1]
            tot += an * mpmath.power(N, s / 2) * mpmath.power(2 * mpmath.pi * n, -s) * mpmath.gammainc(s, n * x0)
            tot += EPS_mp * mpmath.conj(an) * mpmath.power(N, (2 - s) / 2) * mpmath.power(2 * mpmath.pi * n, -(2 - s)) * mpmath.gammainc(2 - s, n * x0)
        return complex(tot / (mpmath.power(N, s / 2) * mpmath.power(2 * mpmath.pi, -s) * mpmath.gamma(s)))
Ld = pari.lfunmf(mf, F)
w = complex(pari.lfunrootres(Ld)[2])
print('PARI root number %s ; lfuncheckfeq %s' % (np.round(w, 10), pari.lfuncheckfeq(Ld)))
Li = pari.lfuninit(Ld, [1, T + 2])
print('lfuninit done (%.0fs)' % (time.time() - t0)); sys.stdout.flush()
def Lp(t): return complex(pari.lfun(Li, pari.complex(1.0, float(t))))
for tc in (7.3, 51.7, 0.75 * T, T):
    lo = L_own(tc); lp = Lp(tc)
    print('   t=%.1f  own %s  PARI %s  |rel diff| %.1e' % (tc, np.round(lo, 8), np.round(lp, 8), abs(lo - lp) / abs(lo)))
sys.stdout.flush()
seps = np.sqrt(eps)
ts = np.arange(0.01, T + 0.5, step)
zv = np.array([Lp(t) for t in ts]) * np.exp(1j * theta(ts)) / seps
print('grid %d points (%.0fs); max |Im Z|/|Z| = %.1e' % (len(ts), time.time() - t0, (np.abs(zv.imag) / np.abs(zv)).max())); sys.stdout.flush()
zr = zv.real
idx = np.nonzero(np.sign(zr[:-1]) * np.sign(zr[1:]) < 0)[0]
def Zr(t): return (Lp(t) * np.exp(1j * theta(t)) / seps).real
roots = []
for i in idx:
    lo, hi, flo = ts[i], ts[i + 1], zr[i]
    for _ in range(25):
        mid = 0.5 * (lo + hi); fm = Zr(mid)
        if np.sign(fm) == np.sign(flo): lo, flo = mid, fm
        else: hi = mid
    roots.append(0.5 * (lo + hi))
g = np.array(roots); g = g[(g > 0) & (g <= T)]
n = len(g)
off = np.arange(1, n + 1) - theta(g) / math.pi
exp_off = 0.5 - np.angle(seps) / math.pi
az = np.abs(zr); sus = [round(ts[i], 2) for i in range(1, len(ts) - 1) if az[i] < az[i - 1] and az[i] < az[i + 1] and (i not in idx and i - 1 not in idx) and az[i] < 0.1 * np.median(az)]
print('zeros %d in (0,%.0f]; n - theta(T)/pi = %.3f ; mean(n - theta(g_n)/pi) = %.3f (expected %.3f) ; drift %.3f ; |Z| minima without sign change below 0.1 median: %s'
      % (n, T, n - theta(T) / math.pi, off.mean(), exp_off, off[-n // 5:].mean() - off[:n // 5].mean(), sus))
dens = np.log(N * (g[:-1] / (2 * math.pi)) ** 2) / (2 * math.pi)
print('min normalised gap %.3f ; first zeros %s' % ((np.diff(g) * dens).min(), np.round(g[:5], 5)))
c = 1.0215617
for TT in (T / 2, T):
    gg = g[g <= TT]
    G = np.exp(2j * theta(gg)).sum(); pred = -eps * c * TT / (2 * math.pi)
    print('T=%.0f: G = %.2f%+.2fi ; law %.2f%+.2fi ; arg G %.1f deg ; law %.1f ; no root number 180.0 ; conj(eps) %.1f ; G/pred = %.3f%+.3fi'
          % (TT, G.real, G.imag, pred.real, pred.imag, np.degrees(np.angle(G)), np.degrees(np.angle(-eps)), np.degrees(np.angle(-np.conj(eps))), (G / pred).real, (G / pred).imag))
wgt = np.sin(math.pi * g / T) ** 2; Gh = (wgt * np.exp(2j * theta(g))).sum(); pred = -eps * c * T / (2 * math.pi)
print('Hann: G_h = %.2f%+.2fi ; law/2 = %.2f%+.2fi ; arg %.1f ; G_h/(pred/2) = %.3f%+.3fi' % (Gh.real, Gh.imag, pred.real / 2, pred.imag / 2, np.degrees(np.angle(Gh)), (Gh / (pred / 2)).real, (Gh / (pred / 2)).imag))
for K in (4, 6, 8, 12):
    edges = np.linspace(0, T, K + 1)
    bk = np.array([np.exp(2j * theta(g[(g > edges[k]) & (g <= edges[k + 1])])).sum() for k in range(K)]) / (pred / K)
    print('  K=%d blocks: G/pred mean %.3f%+.3fi, se(Re) %.3f se(Im) %.3f' % (K, bk.mean().real, bk.mean().imag, bk.real.std(ddof=1) / math.sqrt(K), bk.imag.std(ddof=1) / math.sqrt(K)))
json.dump(dict(char=charn, T=T, eps=[eps.real, eps.imag], zeros=g.tolist()), open(OUT + 'r7b_char%d_T%d.json' % (charn, int(T)), 'w'))
print('done (%.0fs)' % (time.time() - t0))
