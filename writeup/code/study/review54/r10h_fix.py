import math, json, numpy as np, cypari2
from scipy.special import loggamma
pari = cypari2.Pari(); pari.allocatemem(1400 * 10 ** 6)
T = 400.0; k = 24
PD = int(40 + math.pi * T / 2 / math.log(10)) + 15
NC = int((math.pi * T / 2 + PD * 2.31 + 60) / (2 * math.pi) * 1.3) + 80
pari.set_real_precision(PD)
pari('mf = mfinit([1,24],1)'); pari('F = mfeigenbasis(mf)[1]'); pari('v = mfcoefs(F, %d)' % NC); pari('pol = mffields(mf)[1]'); pari('vt = variable(pol)'); pari('rts = polroots(pol)')
best = None
for i in (1, 2):
    pari('zz = real(rts[%d])' % i); l2 = float(pari('subst(lift(v[3]), vt, zz)')) / 2 ** 11.5
    if best is None or l2 > best[0]: best = (l2, i)
pari('zz = real(rts[%d])' % best[1]); pari('W = vector(%d, n, subst(lift(v[n+1]), vt, zz))' % NC)
pari('Ld = lfuncreate([W, 0, [0,1], 24, 1, 1])'); pari('Li = lfuninit(Ld, [12, 0, %f])' % (T + 3))
def theta(t): return 0.5 * t * math.log(1 / math.pi ** 2) + np.imag(loggamma((12 + 1j * t) / 2)) + np.imag(loggamma((13 + 1j * t) / 2))
def Z(t): return (np.exp(1j * theta(t)) * complex(pari('lfun(Li, 12 + (%.15f)*I)' % t))).real
ts = np.arange(398.3, 401.5, 0.004)
z = np.array([Z(t) for t in ts])
idx = np.nonzero(np.sign(z[:-1]) * np.sign(z[1:]) < 0)[0]
print('lambda_2 %.4f ; sign changes in [398.3, 401.5]:' % best[0], np.round([ts[i] for i in idx], 3))
