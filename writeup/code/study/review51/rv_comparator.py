# Comparator: X12 predicted from the MEASURED one-point lines of the two curves (all r = a/b, ab <= M),
# the law only for the lines not measured; the remainder is the continuum. Self pairs X_ii for the scale
# of a fully correlated (isogenous) alternative. Zero counts vs theta/pi (completeness).
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/')
from rv_law import lam, bad, pr, primes, local_rho, w, cE
from rv_analysis import fit, show, Z, PF, PI, PJ, Xmat, ids, nc, c, RAT, N
import mpmath as mp

primes_small = [int(p) for p in primes[:40]]


def rho_gen(k, a, b):
    val = 1.0
    for (n, sgn) in ((a, 1), (b, -1)):
        x = n
        for p in primes_small:
            v = 0
            while x % p == 0:
                x //= p; v += 1
            if v:
                i = int(np.searchsorted(primes, p))
                R, F = local_rho(np.array([lam[k, i]]), np.array([bad[k, i]]), float(p), 8)
                val *= R[sgn * v + 2][0]
        assert x == 1
    return val


ampR = np.array([[-w[k] * c[k] * rho_gen(k, a, b) / np.sqrt(a * b) for (a, b) in RAT] for k in range(nc)])
abR = np.array([a * b for a, b in RAT])
Gw = Z['G_w20_396']; Tw = 376.0
for sig in (0.5, 1.0):
    M = Xmat('X_%s' % sig); ReX = M[PI, PJ].real
    pred = PF['pred_%s' % sig]
    kern = np.array([sig * np.sqrt(2 * np.pi) * np.exp(-sig ** 2 * np.log(a / b) ** 2 / 2) for a, b in RAT])
    print('=== sigma %.1f ===' % sig)
    print('  M    nlines  corr(pred,lawR)  mean(pred-lawR)  rms(pred-lawR) | slope[measR + (pred-lawR)] on pred | ReX on [measR + pred - lawR]: slope | continuum mean | continuum slope on pred')
    units2 = np.c_[PI, PJ]
    for Mx in (1, 2, 6, 12, 20, 30, 60):
        sel = np.where(abR <= Mx)[0]
        lawR = (ampR[PI][:, sel] * ampR[PJ][:, sel] * kern[sel]).sum(1) / (4 * np.pi ** 2)
        measR = ((Gw[PI][:, sel] * np.conj(Gw[PJ][:, sel])).real * kern[sel]).sum(1) / Tw ** 2
        comp = measR + (pred - lawR)
        th1, jk1, sw1, _ = fit(comp, np.c_[np.ones(len(pred)), pred], units2)
        th2, jk2, sw2, _ = fit(ReX, np.c_[np.ones(len(pred)), comp], units2)
        cont = ReX - comp
        th3, jk3, sw3, _ = fit(cont, np.ones((len(pred), 1)), units2)
        th4, jk4, sw4, _ = fit(cont, np.c_[np.ones(len(pred)), pred], units2)
        print('  %3d  %4d   %.4f   %+.4f   %.4f | %.4f +- %.4f | %.4f +- %.4f | %+.4f +- %.4f | %+.4f +- %.4f' % (
            Mx, len(sel), np.corrcoef(pred, lawR)[0, 1], np.mean(pred - lawR), np.std(pred - lawR), th1[1], jk1[1], th2[1], jk2[1], th3[0], jk3[0], th4[1], jk4[1]))
    # data-only comparator (no law at all), M = 30
    sel = np.where(abR <= 30)[0]
    measR = ((Gw[PI][:, sel] * np.conj(Gw[PJ][:, sel])).real * kern[sel]).sum(1) / Tw ** 2
    th, jk, sw, _ = fit(ReX, np.c_[np.ones(len(pred)), measR], units2)
    show('  data only: Re X on measured line products ab<=30 (no law)', th, jk, sw)
    print('  rms(ReX - measR30) %.4f vs rms(ReX - pred) %.4f' % (np.std(ReX - measR), np.std(ReX - pred)))

# self pairs (fully correlated alternative): X_ii by direct double sum
from rv_zeros import Xpair
Xii = np.array([Xpair(k, k, 0.5) for k in range(nc)])
print('self pairs X_ii (sigma 0.5): mean Re %.3f, sd %.3f, min %.3f, max %.3f; density L/2pi at t=200 for N=1000: %.3f' % (
    Xii.real.mean(), Xii.real.std(), Xii.real.min(), Xii.real.max(), np.log(1000 * 200 ** 2 / (4 * np.pi ** 2)) / (2 * np.pi)))

# zero counts vs theta/pi
zz = np.load(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/data/study/r144/fresh_zeros.npz')
dev = []
for k in range(nc):
    z = zz[ids[k]]
    th = lambda t: float(t / 2 * mp.log(N[k]) - t * mp.log(2 * mp.pi) + mp.im(mp.loggamma(1 + 1j * mp.mpf(t))))
    n = int(((z > 20) & (z <= 400)).sum())
    dev.append(n - (th(400) - th(20)) / mp.pi)
dev = np.array(dev, float)
print('zero count in (20,400] minus (theta(400)-theta(20))/pi: mean %.3f, sd %.3f, min %.3f, max %.3f' % (dev.mean(), dev.std(), dev.min(), dev.max()))
