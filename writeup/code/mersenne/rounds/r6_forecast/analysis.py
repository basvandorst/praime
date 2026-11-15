#!/usr/bin/env python3
"""Forecast law for the next Mersenne exponent, PIT / coverage / log scores, #53 forecast.

Input: s1e9.npz from sieve_sum.py (exact prime sums of w below 10^9).
Model W (heuristic): each prime q is an exponent independently with prob ~ w(q),
  w(q) = e^g log2(a q)/q, a = 6 (q = 1 mod 4), 2 (q = 3 mod 4); P(next > x) = exp(-(S(x)-S(p))).
Model P (pure Poisson in log log): a = 1, continuous: Lambda(p, x) = e^g log2(x/p).
Writes results.txt and results.json.
"""
import json, math, sys
import numpy as np
from scipy import stats

EG = math.exp(0.5772156649015329)
C12 = math.log2(math.sqrt(12.0))        # log2(sqrt(12)) = 1.79248...
EXPS = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253,
        4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091,
        756839, 859433, 1257787, 1398269, 2976221, 3021377, 6972593, 13466917, 20996011,
        24036583, 25964951, 30402457, 32582657, 37156667, 42643801, 43112609, 57885161,
        74207281, 77232917, 82589933, 136279841]
VER = 83880077       # all exponents below verified (GIMPS milestones, 2026-11-15)
TESTED = 142869329   # all exponents below tested at least once (same page)
LEVELS = [0.05, 0.25, 0.5, 0.75, 0.95]

d = np.load(sys.argv[1] if len(sys.argv) > 1 else 's1e9.npz')
sq, sW, sP = d['small_q'], d['small_SW'], d['small_SP']
gx, gW, gP, gpi = d['grid_x'], d['grid_SW'], d['grid_SP'], d['grid_pi']
NMAX = int(d['N'])
eS = dict(zip(d['exp_p'].tolist(), d['exp_SW'].tolist()))
eSP = dict(zip(d['exp_p'].tolist(), d['exp_SP'].tolist()))
ew = dict(zip(d['exp_p'].tolist(), d['exp_wW'].tolist()))
ewP = dict(zip(d['exp_p'].tolist(), d['exp_wP'].tolist()))
SMALLMAX = int(sq[-1])
out = []
def P(*a):
    s = ' '.join(str(x) for x in a); out.append(s); print(s)

def S(x, W=True):
    """S(x) = sum_{q <= x} w(q); exact for x <= 10^6 or x a multiple of 1000, else linear interp."""
    x = int(x)
    if x <= 10**6:
        k = np.searchsorted(sq, x, 'right'); return float((sW if W else sP)[k - 1])
    g = gW if W else gP
    j = x // 1000
    if x % 1000 == 0: return float(g[j - 1])
    return float(g[j - 1] + (g[j] - g[j - 1]) * (x % 1000) / 1000)

def quant(S0, level, W=True):
    """Smallest x with S(x) >= S0 + L (L = -ln(1-level)); exact prime below 10^6, interpolated above."""
    T = S0 - math.log1p(-level)
    ss = sW if W else sP
    if T <= ss[-1]:
        k = np.searchsorted(ss, T, 'left'); return int(sq[k])
    g = gW if W else gP
    j = int(np.searchsorted(g, T, 'left'))
    if j >= g.size: return float('nan')
    x0 = gx[j - 1]
    return float(x0 + 1000 * (T - g[j - 1]) / (g[j] - g[j - 1]))

def smooth(x):            # smooth part of S_W: e^g [log2 x + log2(sqrt12) ln ln x]
    return EG * (math.log2(x) + C12 * math.log(math.log(x)))

def smooth_quant(p, level):
    L = -math.log1p(-level); t = smooth(p) + L
    y = p * 2 ** (L / EG)
    for _ in range(50):
        f = smooth(y) - t; dfdy = EG / math.log(2) * (1 + C12 * math.log(2) / math.log(y)) / y
        y -= f / dfdy
    return y

R = {}
# ---------------------------------------------------------------- 0. the integral
P("== 0. Expected count in [x, y]: exact prime sum vs integral (computed) ==")
P("Integral (heuristic, prime density 1/ln q, half of primes in each class mod 4):")
P("  E[x,y] = e^g [ log2(y/x) + log2(sqrt 12) * ln(ln y / ln x) ]")
P("  (the a = 1 'log log' form e^g log2(log M_y/log M_x) = e^g log2(y/x) is the first term)")
P(f"{'x':>10} {'y':>10} {'exact':>9} {'integral':>9} {'e^g log2(y/x)':>13} {'exact/int':>9}")
R['integral_check'] = []
for (x, y) in [(10, 100), (100, 1000), (1000, 10**4), (10**4, 10**5), (10**5, 10**6), (10**6, 10**7),
               (1000, 10**7), (10**7, 10**8), (10**8, 10**9), (10**4, 10**9)]:
    ex = S(y) - S(x)
    it = EG * (math.log2(y / x) + C12 * math.log(math.log(y) / math.log(x)))
    lead = EG * math.log2(y / x)
    P(f"{x:>10} {y:>10} {ex:9.4f} {it:9.4f} {lead:13.4f} {ex/it:9.5f}")
    R['integral_check'].append([x, y, ex, it, lead])
P("Offset D(x) = S_exact(x) - e^g[log2 x + log2(sqrt12) ln ln x]:")
Ds = []
for x in [10**2, 10**3, 10**4, 10**5, 10**6, 10**7, 10**8, 10**9]:
    Ds.append([x, S(x) - smooth(x)]); P(f"  D({x:.0e}) = {Ds[-1][1]:.5f}")
R['D'] = Ds

# ---------------------------------------------------------------- 1. transitions
P("\n== 1-3. Transitions p_n -> p_n+1 (computed under the heuristic model) ==")
trans = []
nver = sum(1 for e in EXPS if e <= VER)   # 51 exponents below the verified bound
for n in range(len(EXPS) - 1):
    p, q = EXPS[n], EXPS[n + 1]
    dS = eS[q] - eS[p]; wq = ew[q]
    u = 1 - math.exp(-dS); ulo = 1 - math.exp(-(dS - wq)); mid = (u + ulo) / 2
    Qs = [quant(eS[p], l) for l in LEVELS]
    # model-exact coverage of [Q25,Q75] and [Q05,Q95] (matters only where the law is discrete)
    def cov(lo_l, hi_l):
        a = quant(eS[p], lo_l); b = quant(eS[p], hi_l)
        if isinstance(b, (int, np.integer)) and b <= SMALLMAX:
            k_a = np.searchsorted(sq, a); Fa_minus = 1 - math.exp(-(sW[k_a - 1] - eS[p]))
            Fb = 1 - math.exp(-(S(b) - eS[p]))
            return Fb - Fa_minus
        return hi_l - lo_l
    in50 = (u >= 0.25) and (ulo < 0.75)
    in90 = (u >= 0.05) and (ulo < 0.95)
    med = Qs[2]
    medP = p * 2 ** (math.log(2) / EG)
    # log scores (bits)
    dSP = eSP[q] - eSP[p]; wqP = ewP[q]
    bitsW_d = (dS - wq) / math.log(2) - math.log2(-math.expm1(-wq))
    bitsP_d = (dSP - wqP) / math.log(2) - math.log2(-math.expm1(-wqP))
    dt = math.log2(q / p)
    bitsP_c = -math.log2(EG) + EG * dt / math.log(2)
    LamW = EG * (dt + C12 * math.log(math.log(q) / math.log(p)))
    hW = EG * (1 + C12 * math.log(2) / math.log(q))
    bitsW_c = -math.log2(hW) + LamW / math.log(2)
    trans.append(dict(n=n + 1, p=p, q=q, verified=(q <= VER), dS=dS, u=u, ulo=ulo, mid=mid,
                      Q=Qs, cov50=cov(0.25, 0.75), cov90=cov(0.05, 0.95), in50=in50, in90=in90,
                      med=med, ratio=q / med, ratioP=q / medP, bitsW_d=bitsW_d, bitsP_d=bitsP_d,
                      bitsW_c=bitsW_c, bitsP_c=bitsP_c))
V = [t for t in trans if t['verified']]
P(f"verified transitions: {len(V)} (p_1 = 2 ... p_{len(V)+1} = {V[-1]['q']}); provisional: {trans[-1]['p']} -> {trans[-1]['q']}")
tot = sum(t['dS'] for t in V)
P(f"sum of Delta S over verified transitions = {tot:.3f} (Poisson expectation of the count 50; time-rescaling: Delta S_n ~ Exp(1) iid)")

# rank by u (largest = most surprising long gap) and by two-sided tail 2 min(u, 1-ulo)
order_u = sorted(V, key=lambda t: -t['u'])
for r, t in enumerate(order_u): t['rank_u'] = r + 1
for t in V: t['two'] = min(1.0, 2 * min(t['u'], 1 - t['ulo']))
order2 = sorted(V, key=lambda t: t['two'])
for r, t in enumerate(order2): t['rank2'] = r + 1

def fmt(x):
    if isinstance(x, (int, np.integer)) or x < 1e6: return f"{x:>11.0f}" if not isinstance(x, (int, np.integer)) else f"{x:>11d}"
    return f"{x:>11.4g}"
P(f"\n{'n':>2} {'p_n':>10} {'p_n+1':>10} {'ratio':>6} {'dS':>6} {'u':>6} {'u_lo':>6} {'rk_u':>4} {'2-sid':>6} {'rk2':>3} "
  f"{'Q05':>11} {'Q25':>11} {'median':>11} {'Q75':>11} {'Q95':>11} {'p/med':>6} {'p/medP':>6} {'bitsW':>6} {'bitsP':>6}")
for t in trans:
    tag = '' if t['verified'] else '  (provisional)'
    P(f"{t['n']:>2} {t['p']:>10} {t['q']:>10} {t['q']/t['p']:6.3f} {t['dS']:6.3f} {t['u']:6.3f} {t['ulo']:6.3f} "
      f"{t.get('rank_u',0):>4} {t.get('two',float('nan')):6.3f} {t.get('rank2',0):>3} "
      + ' '.join(fmt(x) for x in t['Q']) + f" {t['ratio']:6.3f} {t['ratioP']:6.3f} {t['bitsW_c']:6.2f} {t['bitsP_c']:6.2f}{tag}")

# coverage
c50 = sum(t['in50'] for t in V); c90 = sum(t['in90'] for t in V)
e50 = sum(t['cov50'] for t in V); e90 = sum(t['cov90'] for t in V)
pb50 = stats.binomtest(c50, len(V), 0.5).pvalue; pb90 = stats.binomtest(c90, len(V), 0.9).pvalue
P(f"\nCoverage over {len(V)} verified transitions:")
P(f"  50% interval [Q25,Q75]: {c50}/{len(V)} = {c50/len(V):.2f}  (nominal 0.50; model-exact incl. discreteness {e50/len(V):.3f}); binomial p = {pb50:.3f}")
P(f"  90% interval [Q05,Q95]: {c90}/{len(V)} = {c90/len(V):.2f}  (nominal 0.90; model-exact {e90/len(V):.3f}); binomial p = {pb90:.3f}")
below = sum(1 for t in V if t['u'] < 0.05); above = sum(1 for t in V if t['ulo'] >= 0.95)
below25 = sum(1 for t in V if t['u'] < 0.25); above75 = sum(1 for t in V if t['ulo'] >= 0.75)
P(f"  misses: below Q05 {below}, above Q95 {above}; below Q25 {below25}, above Q75 {above75}")
R['coverage'] = dict(n=len(V), c50=c50, c90=c90, model50=e50 / len(V), model90=e90 / len(V), p50=pb50, p90=pb90,
                     below05=below, above95=above, below25=below25, above75=above75)

# uniformity
mid = np.array([t['mid'] for t in V]); uu = np.array([t['u'] for t in V])
ks_all = stats.kstest(mid, 'uniform')
big = [t for t in V if t['p'] >= 127]
ks_big = stats.kstest(np.array([t['u'] for t in big]), 'uniform')
dSarr = np.array([t['dS'] for t in V])
P(f"\nUniformity of the PIT values:")
P(f"  KS (mid-PIT, all {len(V)}): D = {ks_all.statistic:.4f}, p = {ks_all.pvalue:.3f}")
P(f"  KS (u, p_n >= 127, {len(big)} steps, discreteness < 0.04): D = {ks_big.statistic:.4f}, p = {ks_big.pvalue:.3f}")
h, _ = np.histogram(mid, bins=10, range=(0, 1))
P(f"  histogram of mid-PIT in 10 bins (expected {len(V)/10:.0f} each): {h.tolist()}")
chi = stats.chisquare(h); P(f"  chi-square on the 10 bins: stat {chi.statistic:.2f}, p = {chi.pvalue:.3f}")
h5, _ = np.histogram(mid, bins=5, range=(0, 1))
P(f"  5 bins: {h5.tolist()}")
R['ks'] = dict(D_all=ks_all.statistic, p_all=ks_all.pvalue, D_big=ks_big.statistic, p_big=ks_big.pvalue, n_big=len(big),
               hist10=h.tolist(), chi_p=chi.pvalue, sum_dS=tot)

P("\nMost surprising steps (two-sided tail 2 min(u, 1-u_lo); rank among 50):")
for t in order2[:8]:
    P(f"  {t['p']} -> {t['q']}: u = {t['u']:.4f}, u_lo = {t['ulo']:.4f}, two-sided {t['two']:.3f}, rank_u {t['rank_u']}, rank2 {t['rank2']}, p/median {t['ratio']:.3f}")
P("Named steps:")
for (p, q) in [(127, 521), (4423, 9689), (216091, 756839), (2203, 2281), (9689, 9941), (42643801, 43112609)]:
    t = next(t for t in V if t['p'] == p)
    P(f"  {p} -> {q}: u = {t['u']:.4f}, rank_u {t['rank_u']}/50, two-sided {t['two']:.3f} (rank {t['rank2']}), p/median {t['ratio']:.3f}, dS = {t['dS']:.3f}")
mx = max(V, key=lambda t: t['dS'])
P(f"  multiplicity: P(largest of 50 iid Exp(1) gaps >= {mx['dS']:.3f}) = "
  f"{1-(1-math.exp(-mx['dS']))**50:.3f}; P(smallest of 50 iid U <= {min(t['u'] for t in V):.4f}) = {1-(1-min(t['u'] for t in V))**50:.3f}")
R['surprise'] = [dict(p=t['p'], q=t['q'], u=t['u'], ulo=t['ulo'], two=t['two'], rank_u=t['rank_u'], rank2=t['rank2'],
                      ratio=t['ratio']) for t in order2]

# 3. median guesses and log scores
lr = np.array([math.log2(t['ratio']) for t in V]); lrP = np.array([math.log2(t['ratioP']) for t in V])
P("\nPerfect-row scoring (guess = forecast median):")
P(f"  p_n+1/median: geometric mean {2**lr.mean():.3f}, median {2**np.median(lr):.3f}, mean |log2 ratio| {np.abs(lr).mean():.3f}"
  f" (pure-Poisson median 1.3098 p_n: geo mean {2**lrP.mean():.3f}, mean |log2| {np.abs(lrP).mean():.3f})")
P(f"  within a factor 1.25 of the median: {sum(1 for x in lr if abs(x) <= math.log2(1.25))}/50; within 2: {sum(1 for x in lr if abs(x) <= 1)}/50")
P(f"  range of ratios: {2**lr.min():.3f} ({min(V, key=lambda t: t['ratio'])['p']}) to {2**lr.max():.3f} ({max(V, key=lambda t: t['ratio'])['p']})")
for lab, sub in [('all 50', V), ('p_n >= 127 (39)', big)]:
    bWd = sum(t['bitsW_d'] for t in sub); bPd = sum(t['bitsP_d'] for t in sub)
    bWc = sum(t['bitsW_c'] for t in sub); bPc = sum(t['bitsP_c'] for t in sub)
    P(f"  log score {lab}: discrete pmf over primes: W {bWd:.2f} bits, P(a=1) {bPd:.2f} bits, W-P {bWd-bPd:+.2f};"
      f"  density of log2 p_n+1: W {bWc:.2f}, pure Poisson {bPc:.2f}, W-P {bWc-bPc:+.2f}  (per step {bWc/len(sub):.3f} vs {bPc/len(sub):.3f})")
    R['logscore_' + ('all' if len(sub) == len(V) else 'big')] = dict(Wd=bWd, Pd=bPd, Wc=bWc, Pc=bPc, n=len(sub))
P(f"  reference: an exactly calibrated forecast has expected density score per step = -log2(h) + log2(e) = "
  f"{-math.log2(EG)+1/math.log(2):.3f} bits for the pure model (h = e^g per doubling)")

# ---------------------------------------------------------------- 4. #53
P("\n== 4. Forecast for the exponent after 136,279,841 (heuristic) ==")
p52 = 136279841; S52 = eS[p52]; ST = S(TESTED - 1)   # 'below T': T itself untested
Qu = [quant(S52, l) for l in LEVELS]; Qc = [quant(ST, l) for l in LEVELS]
pnone = math.exp(-(ST - S52))
P(f"  unconditional: Q05 {Qu[0]:.4g}, Q25 {Qu[1]:.4g}, median {Qu[2]:.4g}, Q75 {Qu[3]:.4g}, Q95 {Qu[4]:.4g}")
P(f"  P(no exponent in (136279841, {TESTED})) = {pnone:.4f} (expected count {ST-S52:.4f})")
P(f"  conditional on none below {TESTED} (tested once, assumed conclusive): Q05 {Qc[0]:.4g}, Q25 {Qc[1]:.4g}, median {Qc[2]:.4g}, Q75 {Qc[3]:.4g}, Q95 {Qc[4]:.4g}")
P(f"  P(next < 2^28 | none below {TESTED}) = {1-math.exp(-(S(2**28)-ST)):.4f}; P(next < 10^9 | ...) = {1-math.exp(-(S(10**9)-ST)):.4f}")
E28 = S(2**28) - ST; E9 = S(10**9) - ST
E28u = S(2**28) - S52; E9u = S(10**9) - S52
Ever = S(TESTED) - S(VER) ;
P(f"  expected exponents in [{TESTED}, 2^28] = {E28:.4f}; in [{TESTED}, 10^9] = {E9:.4f}")
P(f"  (from 136279841 without the tested range: {E28u:.4f} and {E9u:.4f})")
P(f"  expected exponents in [{VER}, {TESTED}) (tested once, 1 found): {Ever:.4f}")
R['next53'] = dict(uncond=Qu, cond=Qc, pnone=pnone, E28=E28, E9=E9, E28u=E28u, E9u=E9u, Ever=Ever, S52=S52, ST=ST)

# ---------------------------------------------------------------- 5. smooth inversion vs exact (for the JS)
P("\n== 5. Closed form vs exact quantiles (for the web function) ==")
worst = {}
pts = sorted(set([int(x) for x in np.unique(np.round(np.logspace(6, math.log10(2.9e8), 300)))] + [e for e in EXPS if e >= 10**6]))
for p in pts:
    S0 = S(p)
    for l in LEVELS:
        ex = quant(S0, l); sm = smooth_quant(p, l)
        rel = sm / ex - 1
        if abs(rel) > abs(worst.get(l, (0, 0))[0]): worst[l] = (rel, p)
for l in LEVELS: P(f"  level {l}: max relative error of the closed form for p in [1e6, 2.9e8]: {worst[l][0]:+.2e} at p = {worst[l][1]}")
R['closed_form_worst'] = {str(l): worst[l] for l in LEVELS}
R['trans'] = [{k: (v if not isinstance(v, (np.integer,)) else int(v)) for k, v in t.items()} for t in trans]
with open('results.json', 'w') as f: json.dump(R, f, indent=1, default=float)
with open('results.txt', 'w') as f: f.write('\n'.join(out) + '\n')
