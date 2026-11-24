"""Profinite Cramer-Granville model CG(P): the window starts at a random 120-digit integer x0 (so the sieve pattern is a
random profinite position: no Buchstab / sieve-limit effects), sites with no prime factor <= P are 'rough', a rough site is a
pseudo-prime with probability rho = prod_{p<=P}(1-1/p)^(-1) / LREF (fixed density parameter, LREF = ln x of the real data
it is compared with). Statistic as review44/dirs.py (window's own pooled matrix), q = 3, 5, 7, k = 1, H = 1e4, 1e5, 1e6.
Usage: cgsim4.py P LREF seed nwin q1,q2,.. k1,k2  -> appends to cg4_<P>_<LREF>.json (cgsim3 with q, k lists)"""
import numpy as np, math, json, sys, os, random
sys.path.insert(0, "./")
from common import *
P_ = int(sys.argv[1]); LREF = float(sys.argv[2]); seed = int(sys.argv[3]); NW = int(sys.argv[4])
QL = [int(a) for a in sys.argv[5].split(',')]; KL = [int(a) for a in sys.argv[6].split(',')]
def sp(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i*i::i] = False
    return np.flatnonzero(s)
SP = sp(P_); rho = 1.0 / np.prod(1 - 1.0 / SP) / LREF
assert rho < 1, rho
rng = np.random.default_rng(seed); prng = random.Random(seed)
fn = 'cg4_%d_%.2f.json' % (P_, LREF)
res = json.load(open(fn)) if os.path.exists(fn) else {}
PAD = 4000; L = SPAN + 2 * PAD
for w in range(NW):
    x0 = prng.getrandbits(400) | 1            # random ~120-digit start; offsets t = 0..L-1 stand for x0 - PAD + t
    base = x0 - PAD
    pts = []
    for a in range(0, L, 2 * 10 ** 7):
        b = min(a + 2 * 10 ** 7, L); seg = np.ones(b - a, bool)
        for p in SP.tolist(): seg[(-(base + a)) % p::p] = False
        t = a + np.flatnonzero(seg)
        pts.append(t[rng.random(len(t)) < rho])
    T_ = np.concatenate(pts)                  # offsets of pseudo-primes
    inside0 = (T_[:-1] >= PAD) & (T_[:-1] < PAD + SPAN)
    out = {}
    for q, k in [(q, k) for q in QL for k in KL]:
        r = red(q); C = len(r); lut = -np.ones(q, np.int64); lut[r] = np.arange(C)
        c = lut[(T_ + (base % q)) % q]
        a_, b_ = c[:-k], c[k:]; inside = inside0[:len(a_)] & (a_ >= 0) & (b_ >= 0) if k == 1 else ((T_[:-k] >= PAD) & (T_[:-k] < PAD + SPAN) & (a_ >= 0) & (b_ >= 0))
        ai, bi = a_[inside], b_[inside]; Tm = np.bincount(ai * C + bi, minlength=C * C).reshape(C, C).astype(float)
        P = Tm / Tm.sum(1, keepdims=True); pi = Tm.sum(1) / Tm.sum()
        B, labels, par = basis(q, P, pi)
        blk = {}
        for H in (10**4, 10**5, 10**6):
            NB = SPAN // H; bl = ((T_[:-k][inside] - PAD) // H)
            X = np.bincount((bl * C + ai) * C + bi, minlength=NB * C * C).reshape(NB, C, C).astype(float)
            Na = X.sum(2); E = Na[:, :, None] * P[None]
            z = ((X - E) / np.sqrt(E)).reshape(NB, C * C) @ B
            blk[str(H)] = ((z ** 2).mean(0) * NB / (NB - 1)).tolist()
        out['%d_%d' % (q, k)] = dict(labels=labels, parity=par.tolist(), blk=blk, P=P.tolist(), n=int(inside.sum()))
    res['%d_%d' % (seed, w)] = out
    json.dump(res, open(fn, 'w'))
    print('P=%d seed=%d w=%d  ' % (P_, seed, w) + '  '.join('%s: %.4f %.4f %.4f' % (key, *[np.mean(out[key]['blk'][str(H)]) for H in (10**4, 10**5, 10**6)]) for key in out), flush=True)
