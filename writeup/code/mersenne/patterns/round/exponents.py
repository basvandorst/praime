"""Do the 31 Mersenne exponents <= 500,000 concentrate in an arm, angle or ring of the round
layouts more than pseudo-exponents do?  Null: each prime p <= N picked with probability
w(p) = min(1, e^gamma log2(a p)/p), a = 6 (p = 1 mod 4) or 2 (p = 3 mod 4), conditioned on
exactly 31 picks (exact conditional sampling).  Family-wise p by the min-p method."""
import numpy as np, common as C, json, sys
NSETS = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
isp = C.sieve()
P, w = C.null_weights(isp)
null = C.sample_null(P, w, NSETS, 31, seed=20261115)
sets = np.vstack([C.MEXP[None, :], null])          # row 0 = the real exponents
print('null sets', NSETS, ' mean of null sum w:', w.sum().round(2))
FIB = [21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584]
ARMS = {'polar': [6, 19, 25, 44, 333, 710], 'sunflower': FIB, 'wheel': [708]}

# Sacks curve constants: n on curve 4n = y^2 + D, y = round(2 sqrt n)
Pe = np.nonzero(C.sieve(10**5))[0][1:]
def legendre_vec(D, p):
    a = np.mod(D, p).astype(np.int64); e = (p-1)//2; r = np.ones_like(p); b = a.copy()
    while e.any():
        r = np.where(e & 1, (r*b) % p, r); b = (b*b) % p; e >>= 1
    return np.where(a == 0, 0, np.where(r == 1, 1, -1))
def sacks_logC(nums):
    y = np.rint(2*np.sqrt(nums)).astype(np.int64); D = 4*nums - y*y
    out = np.zeros(nums.shape)
    cache = {}
    for i, (d, yy) in enumerate(zip(D.ravel(), y.ravel())):
        if (d, yy % 2) not in cache:
            disc = -d                                # x^2+d: disc -4d = -D ; x^2+x+c: 1-4c = -D
            two = 1.0 if yy % 2 == 0 else 2.0        # odd values on the left curves (c odd)
            leg = legendre_vec(disc, Pe); om = np.where(leg == 0, 1, 1+leg)
            sq = disc >= 0 and int(round(np.sqrt(disc)))**2 == disc
            cache[(d, yy % 2)] = -50.0 if sq else np.log(two) + np.sum(np.log1p(-om/Pe) - np.log1p(-1/Pe))
        out.ravel()[i] = cache[(d, yy % 2)]
    return out, D

def rayleigh(th, k):
    return np.abs(np.exp(1j*k*th).sum(1))**2/th.shape[1]
def max_sector(th, width=np.pi/4):
    s = np.sort(np.mod(th, 2*np.pi), 1); m = s.shape[1]
    ext = np.concatenate([s, s+2*np.pi], 1); best = np.zeros(len(s))
    for i in range(m):
        cnt = (ext[:, i:i+m] < ext[:, i:i+1] + width).sum(1)
        best = np.maximum(best, cnt)
    return best
def max_mult(v):
    out = np.zeros(len(v))
    for j in range(len(v)):
        _, c = np.unique(v[j], return_counts=True); out[j] = c.max()
    return out
def nn_mean(x, y):
    d = np.hypot(x[:, :, None]-x[:, None, :], y[:, :, None]-y[:, None, :])
    d[:, np.arange(x.shape[1]), np.arange(x.shape[1])] = np.inf
    return d.min(2).mean(1)

def pvals(S):
    """S: (1+J) values, larger = more extreme. p for row 0 vs null rows, and for each null row vs
    the other null rows."""
    obs, nul = S[0], S[1:]
    p0 = (1 + np.sum(nul >= obs))/(1 + len(nul))
    srt = np.sort(nul)
    ge = len(nul) - np.searchsorted(srt, nul, 'left')   # count of null >= value (incl. self)
    pn = ge/len(nul)
    return p0, pn

results = {}
allmin = []
for name in C.LAYOUTS:
    L = C.layout(name)
    th = L['theta'][sets]; x = L['x'][sets]; y = L['y'][sets]
    stats = {}
    for k in (1, 2, 3, 4):
        stats[f'rayleigh_k{k}'] = rayleigh(th, k)
    stats['max_sector_45deg'] = max_sector(th)
    stats['nn_dist(-)'] = -nn_mean(x, y)
    if 'ring' in L:
        rv = L['ring'][sets].astype(int)
        vals = np.unique(rv[1:]); E = np.array([(rv[1:] == v).mean(0).sum() if False else (rv[1:] == v).sum()/NSETS for v in vals])
        O = np.stack([(rv == v).sum(1) for v in vals], 1)
        stats['ring_chi2'] = ((O-E)**2/E).sum(1)
        stats['ring_maxcount'] = O.max(1)
    else:
        rr = L['rad'][sets]
        q = np.quantile(rr[1:].ravel(), np.linspace(0, 1, 9)[1:-1])
        b = np.searchsorted(q, rr)
        O = np.stack([(b == i).sum(1) for i in range(8)], 1); E = 31/8
        stats['radial8_chi2'] = ((O-E)**2/E).sum(1)
    if 't' in L:
        stats['count_t<0.25'] = (L['t'][sets] < 0.25).sum(1)
    for m in ARMS.get(name, []):
        stats[f'arm_maxcount_mod{m}'] = max_mult(sets % m)
    if name == 'sacks':
        lc, D = sacks_logC(sets.astype(np.int64))
        stats['curve_mean_logC'] = lc.mean(1)
        stats['curve_maxmult'] = max_mult(D*2 + (np.rint(2*np.sqrt(sets)).astype(int) % 2))
    names = list(stats)
    P0 = []; PN = []
    for s in names:
        p0, pn = pvals(np.asarray(stats[s], float)); P0.append(p0); PN.append(pn)
    P0 = np.array(P0); PN = np.array(PN)
    minobs = P0.min(); minnull = PN.min(0)
    fw = (1 + np.sum(minnull <= minobs))/(1 + NSETS)
    allmin.append((minobs, minnull))
    results[name] = dict(stats={s: dict(obs=float(stats[s][0]), null_mean=float(np.mean(stats[s][1:])), p=float(p)) for s, p in zip(names, P0)},
                         min_p=float(minobs), familywise_p=float(fw))
    best = names[int(np.argmin(P0))]
    print(f"{name:11s} min p = {minobs:.3f} ({best}: obs {stats[best][0]:.2f}, null mean {np.mean(stats[best][1:]):.2f})  family-wise p = {fw:.3f}  [{len(names)} stats]")
# across all layouts: each set's smallest family-wise-free min-p over layouts
mo = np.array([a for a, _ in allmin]); mn = np.stack([b for _, b in allmin])
# calibrate per-layout minima to family-wise scale first, then take the min over layouts
fwn = []
for i in range(len(allmin)):
    srt = np.sort(mn[i]); fwn.append(np.searchsorted(srt, mn[i], 'right')/NSETS)
fwn = np.array(fwn)
fwo = np.array([results[nm]['familywise_p'] for nm in C.LAYOUTS])
overall = (1 + np.sum(fwn.min(0) <= fwo.min()))/(1 + NSETS)
print(f"overall (all 12 layouts, all stats): smallest family-wise p {fwo.min():.3f}; corrected over layouts p = {overall:.3f}")
results['_overall'] = dict(min_layout_fw=float(fwo.min()), overall_p=float(overall), nsets=NSETS)
json.dump(results, open('out/exponents.json', 'w'), indent=1)
