"""Completeness of the zero lists of D = 197 and D = -395 on [600, 1000]; own code.
(1) S(t) = N(t) - theta(t)/pi from the stored zeros and the exact theta (mpmath loggamma);
(2) an independent zero search: sign changes of Z(t) = Lambda real form on a fine grid with PARI lfunhardy,
    plus lfunzeros on the same band; compare counts and positions."""
import os, numpy as np, mpmath as mp, sys
import cypari2
pari = cypari2.Pari()
pari.allocatemem(400 * 10**6)
pari.set_real_precision(38)
R = os.path.dirname(os.path.abspath(__file__))+'/../../../../'
zc = np.load(R + 'writeup/data/riemann/lzeros/r134_C.npz')
zb = np.load(R + 'writeup/data/riemann/lzeros/r138_Bp.npz')


def zeros_of(D):
    k = 'D%d' % D
    return zc[k] if k in zc.files else zb[k]


def theta(t, q, a):
    return t / 2 * mp.log(q / mp.pi) + mp.im(mp.loggamma((0.5 + a + 1j * t) / 2))


for D in [int(x) for x in sys.argv[1:]]:
    q = abs(D); a = 0 if D > 0 else 1
    z = zeros_of(D)
    zz = z[(z >= 600) & (z <= 1000)]
    # S(t) at points just after each zero and just before (extremes of the sawtooth), over whole stored range
    zall = np.sort(z[z > 0])
    pts = []
    for i, g in enumerate(zall):
        if g < 590 or g > 1000: continue
        for t in (g - 1e-7, g + 1e-7):
            N = np.searchsorted(zall, t)
            pts.append((t, N - float(theta(t, q, a)) / np.pi))
    pts = np.array(pts)
    Lmean = np.log(q * 800 / (2 * np.pi))
    gaps = np.diff(zz) * np.log(q * zz[:-1] / (2 * np.pi)) / (2 * np.pi)
    ig = np.argmax(gaps)
    print('D', D, 'zeros in band', len(zz), 'S range on band [%.3f, %.3f]' % (pts[:, 1].min(), pts[:, 1].max()),
          'mean S %.3f' % pts[:, 1].mean(), 'largest unfolded gap %.3f at %.4f-%.4f' % (gaps[ig], zz[ig], zz[ig + 1]))
    # independent: PARI lfunzeros on the band
    L = pari('lfuninit(lfuncreate(%d), [1/2, 1010, 0])' % D)
    pz = pari.lfunzeros(L, pari('[600, 1000]'))
    pz = np.array([float(x) for x in pz])
    print('   PARI lfunzeros count', len(pz), 'stored', len(zz), 'max |diff|',
          np.max(np.abs(pz - zz)) if len(pz) == len(zz) else 'n/a')
    # sign changes of the Hardy function on a fine grid inside the largest gap and around the S extreme
    t0, t1 = zz[ig] - 0.5, zz[ig + 1] + 0.5
    ts = np.linspace(t0, t1, 400)
    hz = np.array([float(pari.lfunhardy(L, t)) for t in ts])
    sc = ts[1:][np.sign(hz[1:]) != np.sign(hz[:-1])]
    print('   Hardy sign changes in [%.3f, %.3f]:' % (t0, t1), np.round(sc, 4),
          'stored zeros there', np.round(zz[(zz > t0) & (zz < t1)], 4), 'min |Z| in gap %.3e' % np.min(np.abs(hz[(ts > zz[ig]+0.05) & (ts < zz[ig+1]-0.05)])))
    jm = np.argmax(np.abs(pts[:, 1]))
    print('   max |S| %.3f at t = %.4f' % (pts[jm, 1], pts[jm, 0]))
