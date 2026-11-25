"""Referee check 0 (own code; does not import gl2lib or check_gl2).
(a) root numbers in fresh_set.json vs PARI ellrootno;
(b) functional equation L(s) = w X(s) L(1-s) at complex s OFF the line (PARI lfun, arithmetic normalisation shifted by 1/2),
    X(s) = N^{1/2-s} (2pi)^{2s-1} Gamma(3/2-s)/Gamma(s+1/2) by mpmath;
(c) X(1/2+it)^{-1} = e^{2 i theta(t)} with theta from PARI lngamma (independent of mpmath);
(d) Stirling: |X(s+it)^{-1}| / (N t^2/4pi^2)^{sigma-1/2} -> 1 with error O(1/t), sigma in {-1/4, 5/4};
(e) the first nontrivial zeros of PARI (lfunzeros) agree with the zero file (completeness spot check in [0, 30]).
"""
import os, json, math, cmath, sys
import numpy as np, mpmath as mp, cypari2

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), *(['..'] * 7)))
CV = {c['id']: c for c in json.load(open(REPO + '/writeup/code/study/lowzero/deg2/fresh_set.json'))['curves']}
ZF = np.load(REPO + '/writeup/data/study/r144/fresh_zeros.npz')
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6)
pari.set_real_precision(38)
mp.mp.dps = 30
IDS = ['N307_0_0_1_m8_m9', 'N335_0_0_1_m2_2', 'N373_0_1_1_m2_m2', 'N611_0_0_1_m1_1', 'N870_1_0_0_m5_m3']

def X(N, s):
    s = mp.mpc(s)
    return mp.power(N, 0.5 - s) * mp.power(2 * mp.pi, 2 * s - 1) * mp.gamma(1.5 - s) / mp.gamma(s + 0.5)

worst = {}
for cid in IDS:
    C = CV[cid]; N = C['N']; w = C['w']
    pari('E = ellinit(%s); Lf = lfuncreate(E);' % str(C['ainv']))
    wp = int(pari('ellrootno(E)')); Np = int(pari('ellglobalred(E)[1]'))
    # (b) FE off the line
    fe = 0
    for s in (mp.mpc(0.8, 37.3), mp.mpc(-0.2, 5.1), mp.mpc(0.31, 101.7), mp.mpc(1.2, 13.0)):
        Ls = complex(pari('lfun(Lf, %s + I*%s)' % (mp.nstr(s.real + 0.5, 20), mp.nstr(s.imag, 20))))
        L1s = complex(pari('lfun(Lf, %s + I*%s)' % (mp.nstr(1 - s.real + 0.5, 20), mp.nstr(-s.imag, 20))))
        fe = max(fe, abs(Ls - w * complex(X(N, s)) * L1s) / abs(Ls))
    # (c) theta via PARI lngamma
    th = 0
    for t in (0.7, 9.3, 55.5, 333.3):
        theta = t / 2 * math.log(N) - t * math.log(2 * math.pi) + float(pari('imag(lngamma(1 + I*%r))' % t))
        th = max(th, abs(1 / complex(X(N, mp.mpc(0.5, t))) - cmath.exp(2j * theta)))
    # (d) Stirling ratio
    st = []
    for sig in (-0.25, 1.25):
        for t in (10.0, 100.0, 400.0):
            r = abs(1 / complex(X(N, mp.mpc(sig, t)))) / (N * t * t / (4 * math.pi ** 2)) ** (sig - 0.5)
            st.append((sig, t, r, abs(math.log(r)) * t))
    # (e) zeros of PARI in [0, 30] vs file
    zp = np.array([float(z) for z in pari('lfunzeros(Lf, 30)')])
    zf = np.sort(ZF[cid]); zf = zf[zf <= 30]
    zp = zp[zp > 1e-6]
    zmatch = (len(zp) == len(zf)) and np.max(np.abs(np.sort(zp) - zf)) < 1e-8
    print('%s N=%d (PARI %d) w=%+d (PARI %+d) | FE off line rel %.1e | X^-1 vs e^{2i theta_PARI} %.1e | t*|ln ratio| max %.3f | zeros<=30: PARI %d file %d match %s (max diff %.1e)' % (
        cid, N, Np, w, wp, fe, th, max(x[3] for x in st), len(zp), len(zf), zmatch,
        np.max(np.abs(np.sort(zp)[:len(zf)] - zf[:len(zp)])) if len(zp) and len(zf) else -1))
    assert Np == N and wp == w and fe < 1e-12 and th < 1e-12 and max(x[3] for x in st) < 1.0 and zmatch
print('r0 all PASS')
