#!/usr/bin/env python3
"""Side numbers: calibration of the a = 1 model, dispersion of Delta S, and the model's own expectation
of the 'perfect row' miss (mean |log2(p_n+1/median)|) by simulation from the exact law (computed)."""
import math, json
import numpy as np
EXPS = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253,
        4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091,
        756839, 859433, 1257787, 1398269, 2976221, 3021377, 6972593, 13466917, 20996011,
        24036583, 25964951, 30402457, 32582657, 37156667, 42643801, 43112609, 57885161,
        74207281, 77232917, 82589933]
d = np.load('s1e9.npz')
eS = dict(zip(d['exp_p'].tolist(), d['exp_SW'].tolist())); eSP = dict(zip(d['exp_p'].tolist(), d['exp_SP'].tolist()))
sq, sW, gx, gW = d['small_q'], d['small_SW'], d['grid_x'], d['grid_SW']
dS = np.array([eS[b] - eS[a] for a, b in zip(EXPS, EXPS[1:])])
dSP = np.array([eSP[b] - eSP[a] for a, b in zip(EXPS, EXPS[1:])])
print(f"sum dS: W {dS.sum():.3f}, a=1 {dSP.sum():.3f} (50 steps); from 127 on: W {dS[11:].sum():.3f}, a=1 {dSP[11:].sum():.3f} (39 steps)")
print(f"var of dS from 127 on: {dS[11:].var(ddof=1):.3f} (Exp(1): 1); mean {dS[11:].mean():.3f}")

def qvec(S0, T):
    """smallest x with S(x) >= S0 + T for an array T (exact primes below 1e6, interpolated above)"""
    tgt = S0 + T
    out = np.empty(T.size)
    small = tgt <= sW[-1]
    out[small] = sq[np.searchsorted(sW, tgt[small], 'left')]
    big = (~small) & (tgt <= gW[-1]); far = tgt > gW[-1]
    j = np.searchsorted(gW, tgt[big], 'left')
    out[big] = gx[j - 1] + 1000 * (tgt[big] - gW[j - 1]) / (gW[j] - gW[j - 1])
    # beyond 10^9 (rare tail draws): extend with the local hazard per doubling e^g (1 + log2(sqrt12)/log2 x)
    h = math.exp(0.5772156649015329) * (1 + math.log2(math.sqrt(12)) / math.log2(1e9))
    out[far] = 1e9 * 2 ** ((tgt[far] - gW[-1]) / h)
    return out

rng = np.random.default_rng(20261115)
R = 4000
sims = np.zeros(R); sims39 = np.zeros(R)
obs = []
for a, b in zip(EXPS, EXPS[1:]):
    med = qvec(eS[a], np.array([math.log(2)]))[0]
    x = qvec(eS[a], rng.exponential(1.0, R))
    m = np.abs(np.log2(x / med))
    sims += m
    if a >= 127: sims39 += m
    obs.append(abs(math.log2(b / med)))
obs = np.array(obs)
sims /= 50; sims39 /= 39
print(f"mean |log2(p_n+1/median)|: observed {obs.mean():.3f} (all 50), model expectation {sims.mean():.3f} "
      f"(sd {sims.std():.3f}); P(model <= observed) = {np.mean(sims <= obs.mean()):.3f}")
print(f"  from 127 on: observed {obs[11:].mean():.3f}, model {sims39.mean():.3f} (sd {sims39.std():.3f}); "
      f"P(model <= observed) = {np.mean(sims39 <= obs[11:].mean()):.3f}")
