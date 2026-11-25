"""Prime-side reconstruction of X (no zeros).
Identity (exact): sum_{g1,g2} e^{2i th1(g1) - 2i th2(g2)} K(g2 - g1) w(g1) w(g2)
                  = (sigma/sqrt(2 pi)) int e^{-sigma^2 xi^2/2} G1(xi) conj G2(xi) dxi,   G(xi) = sum_g e^{2i th(g) - i xi g} w(g).
Explicit formula + stationary phase (each n): G(xi) ~ R(xi) = -(i^a/sqrt q) e^{xi/2} sum_n chi(n) Lambda(n) e^{-i t_n + i/(24 t_n)}
     * [w(t_n) + i (t_n w''(t_n)/2 + w'(t_n)/2 - w(t_n)/(24 t_n))],   t_n = 2 pi n e^xi / q.
Window: Hann on [100, 594]. xi grid: -7.5 .. 7.5 step 0.005 (trapezoid exact up to aliasing at 2 pi/0.005 = 1257 > 494).
usage: s3_prime.py K0 K1 [rand]  -> functions K0..K1-1 of the selected list; saves GR_<k>.npz (G from zeros, R from primes).
'rand': chi(p) replaced by independent random signs (completely multiplicative), seed = k (noise surrogate)."""
import sys, math, time, os
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import numpy as np
from rlib import *

OUT = __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/gr/'
os.makedirs(OUT, exist_ok=True)
W1, W2 = 100.0, 594.0; TW = W2 - W1
XI = np.round(np.arange(-7.5, 7.5 + 1e-9, 0.005), 6)
P = np.load(FZ + 'pred_OO.npz'); Ds = [int(x) for x in P['Ds']]
cand = [i for i in range(181) if abs(Ds[i]) <= 200]
rng = np.random.default_rng(2026)
SEL = sorted(rng.choice(cand, 25, replace=False).tolist())
k0, k1 = int(sys.argv[1]), int(sys.argv[2]); rand = len(sys.argv) > 3 and sys.argv[3] == 'rand'
t0 = time.time()
qmax = max(abs(Ds[i]) for i in SEL[k0:k1])
NMAX = int(qmax * W2 * math.exp(7.5) / (2 * math.pi)) + 10
pr = primes_upto(NMAX)
# prime powers
pp_n, pp_p = [pr], [pr]
for k in range(2, 30):
    lim = int(round(NMAX ** (1.0 / k))) + 2
    ps = pr[pr <= lim]; v = ps.astype(np.int64) ** k; ok = v <= NMAX
    if not ok.any(): break
    pp_n.append(v[ok]); pp_p.append(ps[ok])
nn = np.concatenate(pp_n); pb = np.concatenate(pp_p); kk = np.concatenate([np.full(len(x), i + 1) for i, x in enumerate(pp_n)])
o = np.argsort(nn); nn, pb, kk = nn[o], pb[o], kk[o]; lam = np.log(pb.astype(float))
print('sieve to %d: %d prime powers (%.0fs)' % (NMAX, len(nn), time.time() - t0), flush=True)
Z = load_group(R120F)


def hann(t):
    x = (t - W1) / TW; c = np.cos(2 * math.pi * x); s = np.sin(2 * math.pi * x)
    w = 0.5 * (1 - c); w1 = math.pi / TW * s; w2 = 2 * math.pi ** 2 / TW ** 2 * c
    inside = (t >= W1) & (t <= W2)
    return np.where(inside, w, 0.0), np.where(inside, w1, 0.0), np.where(inside, w2, 0.0)


for k in range(k0, k1):
    i = SEL[k]; D = Ds[i]; q = abs(D); a = 0 if D > 0 else 1
    tab = real_table(D)
    if rand:
        r = np.random.default_rng(1000 + k)
        sgn = np.where(r.random(len(pr)) < 0.5, -1.0, 1.0)
        sp = dict(zip(pr.tolist()[:0], []))
        idx = np.searchsorted(pr, pb); chip = sgn[idx]
        chip = np.where(q % pb == 0, 0.0, chip)
    else:
        chip = tab[pb % q]
    v = lam * chip ** kk
    m = v != 0; n_ = nn[m].astype(float); v_ = v[m]
    # zeros side
    g = Z['D%d' % D]; g = g[(g > W1) & (g < W2)]
    wg = hann(g)[0]; ph = np.exp(2j * theta(g, q, a)) * wg
    G = np.array([np.sum(ph * np.exp(-1j * x * g)) for x in XI])
    # primes side
    R = np.zeros(len(XI), complex)
    pre = -(1j ** a) / math.sqrt(q)
    for j, x in enumerate(XI):
        ex = math.exp(x)
        lo = np.searchsorted(n_, q * W1 / (2 * math.pi * ex)); hi = np.searchsorted(n_, q * W2 / (2 * math.pi * ex), 'right')
        tn = (2 * math.pi * ex / q) * n_[lo:hi]
        w, w1, w2 = hann(tn)
        amp = w + 1j * (tn * w2 / 2 + w1 / 2 - w / (24 * tn))
        R[j] = pre * math.exp(x / 2) * np.sum(v_[lo:hi] * amp * np.exp(1j * (-tn + 1 / (24 * tn))))
    np.savez(OUT + 'GR_%d%s.npz' % (k, '_rand' if rand else ''), G=G, R=R, D=D, q=q, a=a, i=i, xi=XI)
    big = np.abs(XI) <= 4
    print('k=%d D=%d: rms|G| %.2f, rms|G-R| (|xi|<=4) %.3f, (all) %.3f, max %.3f   (%.0fs)' % (
        k, D, np.sqrt(np.mean(np.abs(G) ** 2)), np.sqrt(np.mean(np.abs(G - R)[big] ** 2)), np.sqrt(np.mean(np.abs(G - R) ** 2)),
        np.abs(G - R).max(), time.time() - t0), flush=True)
