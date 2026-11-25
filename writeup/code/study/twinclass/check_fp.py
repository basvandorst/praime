"""Checks of fourpoint.py: (a) the smooth autocorrelation against the exact FFT; (b) every series against direct
local products for sample k; (c) the four-point prediction against the values in r42.json / r43.json."""
import sys; sys.path.insert(0, '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/twinclass')
import math, json, time
import numpy as np
import fourpoint as F
# (a)
for (L0, Q, h, N) in ((10**8, 210, 2, 200000), (10**7, 2310, 4, 38961)):
    C = F.constants(Q, h)
    m = L0 + Q / 2 + Q * np.arange(N); r = C['c'] / (np.log(m) * np.log(m + h))
    M = 1 << int(math.ceil(math.log2(2 * N))); f = np.fft.rfft(r, M); w = np.fft.irfft(f * np.conj(f), M)[:N]
    ws = F.profile_w(L0, Q, h, C['c'], N)(np.arange(N, dtype=float))
    print('(a) w: max rel diff', float(np.max(np.abs(ws[:-10] - w[:-10]) / w[:-10])), 'sum rel diff', float(abs(ws.sum() / w.sum() - 1)))
# (b)
rng = np.random.default_rng(1)
for Q, h in ((210, 2), (2310, 2), (210, 6), (2310, 4), (30030, 2)):
    C = F.constants(Q, h)
    ks = sorted(set([1, 2, 11, 13, 143] + list(rng.integers(1, 40000, 6))))
    worst = 0.0
    for k in ks:
        rho, s0, sp, sm, t1, t2 = F.series(Q, h, k, k + 1, C)
        d = Q * k
        ref = [F.direct(Q, h, d, 'rho'), F.direct(Q, h, d, 'pair'), F.direct(Q, h, d + h, 'pair'), F.direct(Q, h, d - h, 'pair'),
               F.direct(Q, h, d, 'tau'), F.direct(Q, h, d + h, 'tau')]
        got = [rho[0], s0[0], sp[0], sm[0], t1[0], t2[0]]
        worst = max(worst, max(abs(g / r - 1) for g, r in zip(got, ref)))
    print('(b) Q=%d h=%d: worst relative difference over %d k and 6 series: %.2e' % (Q, h, len(ks), worst))
# (c)
r42 = json.load(open('/home/bas/praime/writeup/data/study/r42.json')); r43 = json.load(open('/home/bas/praime/writeup/data/study/r43.json'))
for win, L0, W, Q, ref in (('train', 10**7, 9 * 10**7, 2310, r42['train']['2310']['predicted']),
                           ('train', 10**7, 9 * 10**7, 210, r42['train']['210']['predicted']),
                           ('1e8-1e9', 10**8, 9 * 10**8, 30030, r43['1e8-1e9']['30030']['predicted'])):
    t = time.time(); D = F.decompose(Q, 2, L0, W)
    print('(c) %s Q=%d: own v %.5f, repo %.5f, rel %.2e  (%.1fs)' % (win, Q, D['vX'], ref, D['vX'] / ref - 1, time.time() - t))
