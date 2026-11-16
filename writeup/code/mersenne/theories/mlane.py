"""Theories 1 and 2: the local factor of the Mersenne heuristic.
P(2^p - 1 has no factor q = 2kp + 1, k <= K) = prod_k (1 - P(q prime | p) * [q = +-1 mod 8] / k),
with P(q prime | p) = 2 * prod_{r odd small} (local factor given p mod r) * prod_{r > R} r(r-2)/(r-1)^2 adjusted for r | k, / ln q.
T1 averages over p mod r for r >= 5 (lanes p mod 12); T2 keeps p mod r for r <= 13 (classes mod 4*3*5*7*11*13).
Tests: (a) actual survival of trial factoring to K by class, for primes p in [1e4, 1e5] (K = 800) and [1e5, 8e5] (K = 100);
(b) the 51 exponents below B against the mod-4 null of Wagstaff (a = 2 / 6)."""
import numpy as np, math, sys, time
t0 = time.time()
EXPS = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213, 19937,
        21701, 23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433, 1257787, 1398269, 2976221, 3021377, 6972593,
        13466917, 20996011, 24036583, 25964951, 30402457, 32582657, 37156667, 42643801, 43112609, 57885161, 74207281,
        77232917, 82589933]
B = 83880077
G = math.exp(0.5772156649015329)
NS = 170_000_000
s = np.ones(NS // 2, bool); s[0] = False
for i in range(3, int(NS ** 0.5) + 1, 2):
    if s[i // 2]: s[i * i // 2::i] = False
def isprime(x):
    x = np.asarray(x, np.int64); r = np.zeros(x.shape, bool); o = (x & 1) == 1
    r[o] = s[x[o] // 2]; r[x == 2] = True; return r
P = np.concatenate([[2], 2 * np.nonzero(s)[0] + 1]).astype(np.int64)
print(f'sieve {time.time() - t0:.0f} s, primes below {NS}: {len(P)}', flush=True)
SMALL = [3, 5, 7, 11, 13]
KMAX = 100_000
k = np.arange(1, KMAX + 1, dtype=np.int64)
spf = np.zeros(KMAX + 1, np.int64)
for i in range(2, KMAX + 1):
    if spf[i] == 0: spf[i::i] = np.where(spf[i::i] == 0, i, spf[i::i])
corr = np.ones(KMAX + 1)
for i in range(2, KMAX + 1):
    n = i; last = 0
    while n > 1:
        r = spf[n]
        if r != last and r > 13: corr[i] *= (r - 1) / (r - 2)
        last = r; n //= r
corr = corr[1:]
C = 0.6601618158468696
for r in SMALL: C /= r * (r - 2) / (r - 1) ** 2
def logF(p, res, K, avg=()):
    """log prod_{k<=K} (1 - P_k/k) for a prime p of size p with p = res[r] mod r for r in res (r=4 and odd r in SMALL);
    primes r in avg are averaged over their r-1 nonzero classes (the product averaged, not the logs)."""
    kk = k[:K]; m4 = res[4]
    allowed = np.isin((2 * ((kk * m4) % 4) + 1) % 8, (1, 7))
    base = 2.0 * C * corr[:K] * allowed / np.log(2.0 * kk * p + 1)
    fixed = np.ones(K)
    for r in SMALL:
        if r in avg: continue
        fixed *= np.where((2 * kk * res[r] + 1) % r == 0, 0.0, r / (r - 1))
    if not avg:
        return float(np.sum(np.log1p(-base * fixed / kk)))
    tot = []
    import itertools
    for combo in itertools.product(*[range(1, r) for r in avg]):
        f = fixed.copy()
        for r, c in zip(avg, combo):
            f *= np.where((2 * kk * c + 1) % r == 0, 0.0, r / (r - 1))
        tot.append(np.sum(np.log1p(-base * f / kk)))
    tot = np.array(tot); mx = tot.max()
    return float(mx + np.log(np.mean(np.exp(tot - mx))))
LANES = [1, 5, 7, 11]
def lane_res(L): return {4: L % 4, 3: L % 3}
print('\n== T1 lane factors (relative to the mean of the four lanes), averaged over p mod 5, 7, 11, 13 ==')
print('K convergence at p = 1e6:')
for K in (1000, 3000, 10000):
    lf = np.array([logF(1e6, lane_res(L), K, avg=(5, 7, 11, 13)) for L in LANES])
    f = np.exp(lf - lf.max()); f /= f.mean()
    print(f'  K {K:>6}: ' + '  '.join(f'lane {L}: {x:.4f}' for L, x in zip(LANES, f)))
sizes = [100, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]
LF = {}
for sz in sizes:
    lf = np.array([logF(sz, lane_res(L), 10000, avg=(5, 7, 11, 13)) for L in LANES])
    LF[sz] = lf
    f = np.exp(lf - lf.max()); f /= f.mean()
    m4 = (f[0] + f[1]) / (f[2] + f[3])
    print(f'p {sz:>8.0e}: ' + '  '.join(f'{L}: {x:.3f}' for L, x in zip(LANES, f)) + f'   mod-4 ratio {m4:.3f} (Wagstaff log(6p)/log(2p) = {math.log(6 * sz) / math.log(2 * sz):.3f})', flush=True)
def lane_factor(p, L):
    x = np.log(np.asarray(p, float)); xs = np.log(sizes)
    allf = np.array([LF[sz] for sz in sizes])
    allf = allf - np.log(np.mean(np.exp(allf), axis=1, keepdims=True))
    return np.exp(np.interp(x, xs, allf[:, LANES.index(L)]))
def null(p):
    p = np.asarray(p, float); a = np.where(p % 4 == 1, 6.0, 2.0)
    return np.minimum(G * np.log2(a * p) / p, 1)
def base(p):
    p = np.asarray(p, float)
    return G * (np.log2(2 * p) + 0.5 * np.log2(3)) / p
print('\n== T1 test (b): the 51 exponents below B by lane (p > 3) ==')
PB = P[(P > 3) & (P < B)]
E = np.array([e for e in EXPS if e > 3])
obs = {L: int(np.sum(E % 12 == L)) for L in LANES}
nl = {}; tl = {}
for L in LANES:
    q = PB[PB % 12 == L]
    nl[L] = float(np.sum(null(q)))
    tl[L] = float(np.sum(base(q) * lane_factor(q, L)))
sc = sum(nl.values()) / sum(tl.values())
for L in LANES: tl[L] *= sc
print('lane  observed  mod-4 null  T1 (normalised to the same total)')
for L in LANES: print(f'{L:>4} {obs[L]:>9} {nl[L]:>11.2f} {tl[L]:>9.2f}')
def multi_ll(o, e):
    tot = sum(e.values()); return sum(o[L] * math.log2(e[L] / tot) for L in LANES)
print(f'split log-likelihood (bits): null {multi_ll(obs, nl):.2f}, T1 {multi_ll(obs, tl):.2f}, gain {multi_ll(obs, tl) - multi_ll(obs, nl):+.2f} bits')
rng = np.random.default_rng(1)
n = sum(obs.values()); pt = np.array([tl[L] for L in LANES]); pt /= pt.sum(); pn = np.array([nl[L] for L in LANES]); pn /= pn.sum()
g = []
for _ in range(20000):
    o = rng.multinomial(n, pt); g.append(np.sum(o * np.log2(pt / pn)))
g = np.array(g); og = multi_ll(obs, tl) - multi_ll(obs, nl)
print(f'if T1 were true, the gain would be {g.mean():+.2f} +- {g.std():.2f} bits; observed {og:+.2f}')
g0 = []
for _ in range(20000):
    o = rng.multinomial(n, pn); g0.append(np.sum(o * np.log2(pt / pn)))
g0 = np.array(g0)
print(f'if the mod-4 null were true: {g0.mean():+.2f} +- {g0.std():.2f}; P(gain <= observed) under T1 {np.mean(g <= og):.3f}, P(gain >= observed) under null {np.mean(g0 >= og):.3f}')
print(f'\n({time.time() - t0:.0f} s)', flush=True)
def survival(lo, hi, K):
    p = P[(P >= lo) & (P < hi)]
    alive = np.ones(len(p), bool)
    for kk in range(1, K + 1):
        q = 2 * kk * p + 1
        e = np.ones(len(p), np.int64); b = np.full(len(p), 2, np.int64); x = p.copy()
        while np.any(x > 0):
            m = (x & 1) == 1
            e = np.where(m, (e * b) % q, e); b = (b * b) % q; x >>= 1
        alive &= ~(e == 1)
    return p, alive
def T2logF(p, K):
    """per-prime log factor with p mod r kept for r in SMALL (actual residues), r > 13 averaged."""
    out = np.empty(len(p))
    cache = {}
    for i, pp in enumerate(p):
        key = (int(pp % 4),) + tuple(int(pp % r) for r in SMALL)
        if key not in cache:
            res = {4: key[0]}; res.update({r: key[j + 1] for j, r in enumerate(SMALL)})
            cache[key] = logF(float(np.exp(np.mean(np.log(p)))), res, K)
        out[i] = cache[key]
    return out
print('\n== T1/T2 test (a): actual survival of trial factoring ==')
for lo, hi, K in ((10_000, 100_000, 800), (100_000, 800_000, 100)):
    p, alive = survival(lo, hi, K)
    print(f'primes in [{lo:,}, {hi:,}), K = {K}: {len(p)} primes, survivors {alive.sum()} ({alive.mean():.4f})')
    pm = float(np.exp(np.mean(np.log(p))))
    lf = np.array([logF(pm, lane_res(L), K, avg=(5, 7, 11, 13)) for L in LANES])
    print('  lane  primes  survival  T1 predicted (absolute)  ratio obs/pred')
    for L, x in zip(LANES, lf):
        m = p % 12 == L
        print(f'  {L:>4} {m.sum():>7} {alive[m].mean():>9.4f} {math.exp(x):>12.4f} {alive[m].mean() / math.exp(x):>16.3f}')
    t2 = np.exp(T2logF(p, K))
    print(f'  T2 per-prime prediction: mean {t2.mean():.4f} vs survival {alive.mean():.4f}; correlation {np.corrcoef(t2, alive)[0, 1]:.3f}')
    qs = np.quantile(t2, np.linspace(0, 1, 11))
    print('  deciles of the T2 prediction: predicted mean / observed survival / n')
    for a, b in zip(qs[:-1], qs[1:]):
        m = (t2 >= a) & (t2 <= b)
        print(f'    {t2[m].mean():.4f} / {alive[m].mean():.4f} / {m.sum()}')
    ll_t2 = np.sum(np.where(alive, np.log2(t2), np.log2(1 - t2)))
    lane_pred = np.exp(np.array([lf[LANES.index(int(x))] for x in p % 12]))
    ll_t1 = np.sum(np.where(alive, np.log2(lane_pred), np.log2(1 - lane_pred)))
    c = alive.mean(); ll_0 = np.sum(np.where(alive, np.log2(c), np.log2(1 - c)))
    print(f'  log-likelihood of who survives (bits): constant {ll_0:.1f}, T1 lanes {ll_t1:.1f}, T2 per prime {ll_t2:.1f}', flush=True)
print(f'\n== T2 test (b): the 51 exponents with the per-prime factor (within each p mod 4 class, normalised) ==')
def t2_weights(q, K=2000):
    keys = np.stack([q % 4] + [q % r for r in SMALL], 1)
    uk, inv = np.unique(keys, axis=0, return_inverse=True)
    out = np.zeros(len(q))
    for sz_lo, sz_hi in ((4, 1e3), (1e3, 1e5), (1e5, 1e7), (1e7, B)):
        m = (q >= sz_lo) & (q < sz_hi)
        if not m.any(): continue
        pm = float(np.exp(np.mean(np.log(q[m]))))
        vals = np.array([logF(pm, {4: int(u[0]), **{r: int(u[j + 1]) for j, r in enumerate(SMALL)}}, K) for u in uk])
        w = np.exp(vals[inv[m]])
        for c4 in (1, 3):
            mm = (q[m] % 4) == c4
            w[mm] /= w[mm].mean()
        out[m] = w
    return out
w = t2_weights(PB)
pn = null(PB); pt2 = pn * w
ise = np.isin(PB, E)
ll0 = np.sum(np.log2(pn[ise])) + np.sum(np.log2(1 - pn[~ise]))
ll2 = np.sum(np.log2(pt2[ise])) + np.sum(np.log2(1 - np.minimum(pt2[~ise], 0.999)))
print(f'expected total: null {pn.sum():.2f}, T2 {pt2.sum():.2f}; log-likelihood gain of T2 over the null: {ll2 - ll0:+.2f} bits')
print('T2 weight at the exponents: ' + ', '.join(f'{int(e)}:{x:.2f}' for e, x in zip(PB[ise], w[ise])))
sims = []
for _ in range(2000):
    lab = rng.random(len(PB)) < pn
    sims.append(np.sum(np.log2(w[lab])))
sims = np.array(sims); obs_s = np.sum(np.log2(w[ise]))
simt = []
for _ in range(2000):
    lab = rng.random(len(PB)) < np.minimum(pt2, 1)
    simt.append(np.sum(np.log2(w[lab])))
simt = np.array(simt)
print(f'sum of log2 weights at the exponents: {obs_s:+.2f}; under the null {sims.mean():+.2f} +- {sims.std():.2f} (P >= obs {np.mean(sims >= obs_s):.3f}); under T2 {simt.mean():+.2f} +- {simt.std():.2f} (P <= obs {np.mean(simt <= obs_s):.3f})')
print(f'({time.time() - t0:.0f} s)')
