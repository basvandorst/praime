import os, sys, math, json, numpy as np, cypari2
from scipy.special import loggamma
pari = cypari2.Pari(); pari.allocatemem(1400 * 10 ** 6)
OUT = os.path.dirname(os.path.abspath(__file__))+'/'
T = 600.0; N = 13
PD = int(40 + math.pi * T / 2 / math.log(10)) + 15
NC = int((math.pi * T / 2 + PD * 2.31) / (2 * math.pi / math.sqrt(N)) * 1.3) + 50
pari.set_real_precision(PD)
pari('mf = mfinit([13,2,Mod(4,13)],0)'); pari('F = mfeigenbasis(mf)[1]'); pari('v = mfcoefs(F, %d)' % NC); pari('vt = variable(lift(v[3]))')
pari('z = exp(I*Pi/3)'); pari('W = vector(%d, n, subst(lift(v[n+1]), vt, z))' % NC); pari('Wc = conj(W)')
pari('th(y, A) = sum(n=1, #A, A[n]*exp(-2*Pi*n*y))'); pari('y0 = 0.9/sqrt(13)'); pari('EPS = th(1/(13*y0), W)/(13*y0^2*th(y0, Wc))')
eps = complex(pari('EPS'))
pari('Ld = lfuncreate([W, Wc, [0,1], 2, 13, EPS])'); pari('Li = lfuninit(Ld, [1, 0, %f])' % (T + 2))
def theta(t):
    t = np.asarray(t, float)
    return 0.5 * t * math.log(N / math.pi ** 2) + np.imag(loggamma((1.0 + 1j * t) / 2)) + np.imag(loggamma((2.0 + 1j * t) / 2))
def Zd(t): return (np.exp(1j * theta(t)) * np.conj(complex(pari('lfun(Li, 1 + (%.15f)*I)' % (-t)))) / np.conj(np.sqrt(eps))).real
ts = np.arange(431.9, 432.7, 0.002)
z = np.array([Zd(t) for t in ts])
idx = np.nonzero(np.sign(z[:-1]) * np.sign(z[1:]) < 0)[0]
new = [float(ts[i] - z[i] * (ts[i + 1] - ts[i]) / (z[i + 1] - z[i])) for i in idx]
print('fine grid near 432.3: sign changes at', np.round(new, 4), ' min |Z| %.3e at %.3f' % (np.abs(z).min(), ts[np.argmin(np.abs(z))]))
r = json.load(open(OUT + 'r7g_n13_400_600.json'))
old = np.array(r['res']['N13d']['zeros'])
print('zeros of the coarse run within [431.9, 432.7]:', old[(old > 431.9) & (old < 432.7)])
allz = np.sort(np.concatenate([old, [x for x in new if np.min(np.abs(old - x)) > 1e-3]]))
r['res']['N13d']['zeros'] = allz.tolist(); r['fixed'] = new
json.dump(r, open(OUT + 'r7g_n13_400_600_fixed.json', 'w'))
print('N13d zeros in (400,600] now', len(allz))
