# Development curves (20) with the reviewer's own code: one-point lines (13 r), pooled slope; combined with fresh.
import os, numpy as np, sys
import cypari2, mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/')
SC = os.path.dirname(os.path.abspath(__file__))+'/'
pari = cypari2.Pari(); pari.allocatemem(300 * 10**6)
D = np.load(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/data/study/r148/dev_zeros.npz')
dids = sorted(D.keys())
P = 200000
primes = np.array([int(x) for x in pari('primes(primepi(%d))' % P)])
pr = primes.astype(float)
LINES = [(1, 1), (2, 1), (1, 2), (3, 1), (1, 3), (4, 1), (1, 4), (5, 1), (1, 5), (7, 1), (1, 7), (3, 2), (2, 3)]
from rv_law import local_rho
amps, Gs, ranks = [], [], []
for s in dids:
    parts = s.split('_'); a = [int(x.replace('m', '-')) for x in parts[1:]]
    E = pari.ellinit(a); Nn = int(pari.ellglobalred(E)[0]); w = int(pari.ellrootno(E)); assert Nn == int(parts[0][1:])
    ranks.append(int(pari.ellanalyticrank(E)[0]))
    pari('global(EE); EE = ellinit(%s)' % str(a))
    ap = np.array([int(x) for x in pari('my(v=vectorsmall(%d),i=0);forprime(p=2,%d,i++;v[i]=ellap(EE,p));v' % (len(primes), P))], float)
    lam = ap / np.sqrt(pr); bad = (Nn % primes == 0)
    L2 = lam * lam - 1
    F = np.where(bad, 1 - 1 / pr ** 2, (1 - 1 / pr) * (1 - L2 / pr))
    c = float(np.exp(np.sum(np.log(F / (1 - 1 / pr)))))
    amp = []
    for (aa, bb) in LINES:
        val = 1.0
        for (n, sg) in ((aa, 1), (bb, -1)):
            x = n
            for p in (2, 3, 5, 7):
                v = 0
                while x % p == 0:
                    x //= p; v += 1
                if v:
                    i = int(np.searchsorted(primes, p))
                    R, _ = local_rho(np.array([lam[i]]), np.array([bad[i]]), float(p), 6)
                    val *= R[sg * v + 2][0]
        amp.append(-w * c * val / np.sqrt(aa * bb))
    amps.append(amp)
    z = np.sort(D[s]); g = z[(z > 20) & (z <= 400)]
    th = np.array([float(t / 2 * mp.log(Nn) - t * mp.log(2 * mp.pi) + mp.im(mp.loggamma(1 + 1j * mp.mpf(t)))) for t in g])
    e = np.exp(2j * th)
    Gs.append([np.sum(e * np.exp(-1j * np.log(aa / bb) * g)) for (aa, bb) in LINES])
amps = np.array(amps) * 380 / (2 * np.pi); Gs = np.array(Gs)
np.savez(SC + 'rv_dev.npz', ids=np.array(dids), amps=amps, G=Gs, rank=np.array(ranks))


def slope_jk(Y, X):
    # rows = curves, columns = lines; through 0; jackknife over curves
    n = Y.shape[0]
    b = np.sum(Y * X) / np.sum(X * X)
    bk = np.array([(np.sum(Y * X) - np.sum(Y[k] * X[k])) / (np.sum(X * X) - np.sum(X[k] * X[k])) for k in range(n)])
    return b, np.sqrt((n - 1) / n * np.sum((bk - bk.mean()) ** 2))


Fz = np.load(SC + 'rv_zeros.npz'); Fl = np.load(SC + 'rv_law.npz')
RAT = [tuple(x) for x in Fz['rat']]; ri = [RAT.index(x) for x in LINES]
GF = Fz['G_b20_400'][:, ri]; AF = Fl['amp'] * 380 / (2 * np.pi)
for nm, Y, X in (('dev r=1', Gs[:, :1].real, amps[:, :1]), ('dev 12 lines', Gs[:, 1:].real, amps[:, 1:]), ('dev 13 lines', Gs.real, amps),
                 ('fresh 13 lines', GF.real, AF), ('all 158 r=1', np.r_[Gs[:, :1].real, GF[:, :1].real], np.r_[amps[:, :1], AF[:, :1]]),
                 ('all 158, 12 lines', np.r_[Gs[:, 1:].real, GF[:, 1:].real], np.r_[amps[:, 1:], AF[:, 1:]]),
                 ('all 158, 13 lines', np.r_[Gs.real, GF.real], np.r_[amps, AF])):
    b, s = slope_jk(Y, X)
    print('%-20s slope %.4f +- %.4f  (%.1f se from 1)' % (nm, b, s, (b - 1) / s))
