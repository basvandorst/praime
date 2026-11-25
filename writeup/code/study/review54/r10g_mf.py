# Weight-24 level-1 eigenforms (Hecke field Q(sqrt 144169)), exact embedding at high precision; zeros from sign changes of Z.
# The law predicts c = -0.0858 for the embedding with lambda_2 = 1.7596 (sum POSITIVE for eps = +1) and c = 0.527 for the other.
import os, sys, math, json, time, numpy as np, cypari2
from scipy.special import loggamma
pari = cypari2.Pari(); pari.allocatemem(1400 * 10 ** 6)
OUT = os.path.dirname(os.path.abspath(__file__))+'/'
T = float(sys.argv[1]); step = float(sys.argv[2]); which = sys.argv[3]
NL = int(sys.argv[4]); k = int(sys.argv[5]); cfro = float(sys.argv[6]); tag = sys.argv[7]
PD = int(40 + math.pi * T / 2 / math.log(10)) + 15
NC = int((math.pi * T / 2 + PD * 2.31 + 60) / (2 * math.pi / math.sqrt(NL)) * 1.3) + 80
pari.set_real_precision(PD)
t0 = time.time()
pari('mf = mfinit([%d,%d],0)' % (NL, k)); pari('B = mfeigenbasis(mf)'); pari('F = B[1]'); pari('v = mfcoefs(F, %d)' % NC)
pari('pol = mffields(mf)[1]'); pari('vt = variable(pol)')
print('Hecke field polynomial', pari('pol'), ' a(2) =', pari('v[3]'))
pari('rts = polroots(pol)')
rts = [complex(pari('rts[%d]' % (i + 1))) for i in range(2)]
# choose the root giving lambda_2 = 1.7596 ('+') or -1.3867 ('-')
cand = []
for i in range(2):
    pari('zz = rts[%d]' % (i + 1))
    a2 = complex(pari('subst(lift(v[3]), vt, zz)')).real
    cand.append((a2 / 2 ** ((k - 1) / 2), i))
cand.sort()
idx = cand[-1][1] if which == '+' else cand[0][1]
pari('zz = real(rts[%d])' % (idx + 1))
pari('W = vector(%d, n, subst(lift(v[n+1]), vt, zz))' % NC)
lam2 = float(pari('W[2]')) / 2 ** ((k - 1) / 2)
print('chosen embedding: lambda_2 = %.4f ; a(1) = %s ; check a(6) - a(2)a(3) = %s' % (lam2, pari('W[1]'), pari('W[6]-W[2]*W[3]')))
pari('Ld0 = lfuncreate([W, 0, [0,1], %d, %d, 1])' % (k, NL)); pari('EPSr = lfunrootres(Ld0)[3]'); print('root number', pari('EPSr')); pari('Ld = lfuncreate([W, 0, [0,1], %d, %d, round(real(EPSr))])' % (k, NL))
print('PD %d NC %d lfuncheckfeq %s (%.0fs)' % (PD, NC, pari('lfuncheckfeq(Ld)'), time.time() - t0)); sys.stdout.flush()
pari('Li = lfuninit(Ld, [%f, 0, %f])' % (k / 2, T + 2))
print('lfuninit (%.0fs)' % (time.time() - t0)); sys.stdout.flush()
mu = [(k - 1) / 2, (k + 1) / 2]
def theta(t):
    t = np.asarray(t, float)
    return 0.5 * t * math.log(NL / math.pi ** 2) + sum(np.imag(loggamma((0.5 + m + 1j * t) / 2)) for m in mu)
EPSv = float(pari('round(real(EPSr))'))
def Z(t): return np.exp(1j * theta(t)) * complex(pari('lfun(Li, %f + (%.15f)*I)' % (k / 2, t))) / (1 if EPSv > 0 else 1j)
ts = np.arange(0.01, T + 0.5, step)
zv = np.array([Z(t) for t in ts]); zr = zv.real
print('grid %d (%.0fs), max |Im Z|/|Z| = %.1e' % (len(ts), time.time() - t0, (np.abs(zv.imag) / np.abs(zv)).max())); sys.stdout.flush()
idx = np.nonzero(np.sign(zr[:-1]) * np.sign(zr[1:]) < 0)[0]
roots = []
for i in idx:
    lo, hi, flo = ts[i], ts[i + 1], zr[i]
    for _ in range(22):
        mid = 0.5 * (lo + hi); fm = Z(mid).real
        if np.sign(fm) == np.sign(flo): lo, flo = mid, fm
        else: hi = mid
    roots.append(0.5 * (lo + hi))
g = np.array(roots); g = g[(g > 0) & (g <= T)]
n = len(g); off = np.arange(1, n + 1) - theta(g) / math.pi
az = np.abs(zr)
sus = [round(float(ts[i]), 2) for i in range(1, len(ts) - 1) if az[i] < az[i - 1] and az[i] < az[i + 1] and (i not in idx and i - 1 not in idx) and az[i] < 0.15 * np.median(az[max(0, i - 50):i + 50])]
print('zeros %d ; n - theta(T)/pi = %.3f ; mean(n - theta/pi) = %.3f (expected 0.5 for eps=+1, 0 for -1) ; drift %.3f ; dips %s' % (n, n - theta(T) / math.pi, off.mean(), off[-n // 6:].mean() - off[:n // 6].mean(), sus))
c = cfro
G = np.exp(2j * theta(g)).sum(); pred = -EPSv * c * T / (2 * math.pi)
print('%s lambda_2 %.4f c %.6f eps %+d : G = %.2f%+.2fi ; law -eps c T/2pi = %+.2f ; c=1 -> %.1f ; slope Re(G)/pred = %.3f' % (tag, lam2, c, EPSv, G.real, G.imag, pred, -EPSv * T / (2 * math.pi), G.real / pred))
for K in (4, 6, 8, 12, 20):
    edges = np.linspace(0, T, K + 1)
    Gk = np.array([np.exp(2j * theta(g[(g > edges[j]) & (g <= edges[j + 1])])).sum() for j in range(K)])
    print('   K=%2d blocks: se of Re(sum) %.2f' % (K, Gk.real.std(ddof=1) * math.sqrt(K)))
json.dump(dict(T=T, lam2=lam2, c=c, zeros=g.tolist(), G=[G.real, G.imag]), open(OUT + 'r10g_%s_T%d.json' % (tag, int(T)), 'w'))
print('done (%.0fs)' % (time.time() - t0))
