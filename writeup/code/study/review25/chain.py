"""Twin class-chain comparator: per-class one-gap twin law (round 59's z^0 term, per class a mod 30030), first-order
chain on the 1485 twin-start classes, k-step distance distribution; window-averaged like r89. Uses r59.lnS_parts only
for the S-ratios (already cross-checked against my own S to 1e-6)."""
import sys, math, json, time
sys.path.insert(0, '/home/bas/praime/writeup/code/study')
import numpy as np, r59
from scipy import stats
t0 = time.time()
CLS = r59.CLS; PP = r59.PP; idx = np.full(PP, -1); idx[CLS] = np.arange(len(CLS))
lams = [math.log(1e9), math.log(1.5e9), math.log(2e9)]
ND = 200; ds = np.arange(6, 1201, 6)
lnS02 = sum(math.log(1 - 2 / q) - 2 * math.log(1 - 1 / q) for q in r59.MID) + r59.TAIL[2]
F = np.zeros((3, len(CLS), ND))
for j, d in enumerate(ds):
    h, lnB, lns, lnr = r59.lnS_parts(int(d))
    ok = (np.gcd(CLS + d, PP) == 1) & (np.gcd(CLS + d + 2, PP) == 1); cls = CLS[ok]
    if not len(cls): continue
    M = ((np.gcd(cls[:, None] + h[None, :], PP) == 1) & (np.gcd(cls[:, None] + h[None, :] + 2, PP) == 1)).astype(float) if len(h) else np.zeros((len(cls), 0))
    R = np.expm1(lnr)
    for li, lam in enumerate(lams):
        rho0 = 1 / (lam * r59.DEN); p = rho0 ** 2 * np.exp(lns)
        L = M @ np.log1p(-p) if len(h) else np.zeros(len(cls))
        if len(h) > 1:
            w = p / (1 - p); L = L + 0.5 * np.sum((M @ (w[:, None] * w[None, :] * R)) * M, axis=1)
        F[li, idx[cls], j] = rho0 ** 2 * math.exp(lnB - lnS02) * np.exp(L)
print('table', round(time.time() - t0), 's; row sums', F[1].sum(1).min(), F[1].sum(1).max(), flush=True)
# check: class average equals rlaw z^0
r = np.load('rlaw.npz'); print('class-avg vs r89 z0 max rel', np.max(np.abs(F.mean(1) / np.where(r['out'][0] > 0, r['out'][0], 1) - 1)[r['out'][0] > 0]))
# chain: state vector over (class, distance)
dist = {k: np.zeros((3, ND)) for k in (1, 2, 3)}
nxt = (CLS[:, None] + ds[None, :]) % PP; nidx = idx[nxt]     # -1 where not a twin class
for li in range(3):
    V = F[li] / len(CLS)          # V[a', D] = prob mass of being at twin start of class a' after 1 step at distance D: index by start class
    # represent as mass at end class: E[c, D]
    E = np.zeros((len(CLS), ND))
    for j in range(ND):
        m = nidx[:, j] >= 0
        np.add.at(E[:, j], nidx[m, j], V[m, j])
    dist[1][li] = E.sum(0)
    for k in (2, 3):
        E2 = np.zeros_like(E)
        for j in range(ND):          # extend by one gap of distance ds[j2]
            if not E[:, j].any(): continue
            for j2 in range(ND - j - 1):
                jj = j + j2 + 1      # distance ds[j] + ds[j2] = 6(j+1) + 6(j2+1) -> index j + j2 + 1
                w = E[:, j] * F[li, :, j2]
                m = (w > 0) & (nidx[:, j2] >= 0)
                np.add.at(E2[:, jj], nidx[m, j2], w[m])
        E = E2; dist[k][li] = E.sum(0)
    print('lam', li, round(time.time() - t0), 's', flush=True)
u = np.linspace(lams[0], lams[-1], 401); wt = np.exp(u) / u ** 2
cnt = json.load(open('counts.json'))
BANDS = [(lo, lo + 24) for lo in range(6, 600, 30)] + [(lo, lo + 54) for lo in range(606, 1200, 60)]
for k in (1, 2, 3):
    A = np.array([np.sum(np.interp(u, lams, dist[k][:, j]) * wt) / wt.sum() for j in range(ND)]); A /= A.sum()
    c = cnt[f'{k}_in']; N = sum(c.values()); x = 0; dof = 0; zm = 0
    for lo, hi in BANDS:
        o = sum(c.get(str(d), 0) for d in range(lo, hi + 1, 6)); e = N * A[(ds >= lo) & (ds <= hi)].sum()
        if e >= 50: x += (o - e) ** 2 / e; dof += 1; zm = max(zm, abs(o - e) / math.sqrt(e))
    print('k', k, 'class-chain chi2 %.1f dof %d maxz %.1f q999 %.1f' % (x, dof - 1, zm, stats.chi2.ppf(.999, dof - 1)), flush=True)
