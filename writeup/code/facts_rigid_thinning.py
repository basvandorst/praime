import numpy as np, math
from facts_maxgap_is import PR, logprod
X, L = 10**9, 10**8; r = math.isqrt(X + L); rng = np.random.default_rng(11)
real = np.ones(L, bool); S = np.ones(L, bool)
for q in PR[PR <= r]:
    q = int(q); real[(-X) % q::q] = False; S[rng.integers(q)::q] = False
rho = 1.0 / (math.log(X) * math.exp(logprod(1, r)))
pos = np.nonzero(S)[0]
def thin(B):
    if B == 0: return pos[rng.random(len(pos)) < rho]                       # independent coin
    blk = pos // B; pri = rng.random(len(pos)); o = np.lexsort((pri, blk)); b = blk[o]
    start = np.searchsorted(b, b, 'left'); rank = np.arange(len(b)) - start
    n = np.bincount(blk, minlength=L // B + 1); k = np.floor(rho * n + rng.random(len(n))).astype(int)
    keep = rank < k[b]; return np.sort(pos[o][keep])
def stats(P):
    ind = np.zeros(L, bool); ind[P] = True; out = []
    for H in (10**4, 10**5, 10**6):
        s = ind[:L//H*H].reshape(-1, H).sum(1).astype(float); t = np.arange(len(s)); c = np.polyfit(t, s, 1)
        out.append(np.var(s - np.polyval(c, t)) / s.mean())
    g = np.diff(P); h = np.bincount(g, minlength=400)[:400] / len(g)
    return out, h, (g >= 100).mean(), (g >= 150).mean()
Rv, Rh, R100, R150 = stats(np.nonzero(real)[0])
print(f"real            var/mean H=1e4,1e5,1e6: {Rv[0]:.3f} {Rv[1]:.3f} {Rv[2]:.3f}   P(g>=100)={R100:.5f} P(g>=150)={R150:.6f}")
for B in (0, 1000, 10**4, 10**5):
    v, h, p100, p150 = stats(thin(B))
    chi = sum((h[g] - Rh[g])**2 / Rh[g] for g in range(2, 120, 2)) * len(np.nonzero(real)[0])
    print(f"thin B={B:>6}: var/mean {v[0]:.3f} {v[1]:.3f} {v[2]:.3f}   P(g>=100)={p100:.5f} P(g>=150)={p150:.6f}  chi2(gaps 2..118, 59 bins)={chi:.0f}", flush=True)
