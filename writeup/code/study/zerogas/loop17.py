import sys, numpy as np
from sine import p1
rng = np.random.default_rng(17); N = int(sys.argv[1]); M = int(sys.argv[2]); K = N//2
ks = np.arange(1, K + 1); G = 20; blk = M//G
sV = np.zeros(G); sP = np.zeros((G, K)); sPP = np.zeros((G, K)); sVP = np.zeros((G, K)); sC = np.zeros((G, K))
ph = -2j*np.pi*ks[:, None]/N
for it in range(M):
    Z = (rng.standard_normal((N, N)) + 1j*rng.standard_normal((N, N)))/np.sqrt(2); Q, R = np.linalg.qr(Z)
    x = np.sort(np.angle(np.linalg.eigvals(Q*(np.diag(R)/np.abs(np.diag(R)))[None, :])) % (2*np.pi))*N/(2*np.pi)
    s = np.diff(np.concatenate([x, [x[0] + N]])); V = float(np.mean(s*s))
    P = np.abs(np.exp(ph*x[None, :]).sum(axis=1))**2
    g = it//blk
    if g >= G: break
    sV[g] += V; sP[g] += P; sPP[g] += P*P; sVP[g] += V*P
    sC[g] += np.cos(2*np.pi*np.outer(ks, s)/N).mean(axis=1)
def slope(idx):
    n = blk*len(idx); mV = sV[idx].sum()/n; mP = sP[idx].sum(0)/n
    return (sVP[idx].sum(0)/n - mV*mP)/(sPP[idx].sum(0)/n - mP*mP), 1 - sC[idx].sum(0)/n
sl, oc = slope(np.arange(G)); cl = oc/(np.pi**2*ks**2)
# jackknife over the 20 groups, for bands of modes
bands = [(0.03, 0.08), (0.08, 0.13), (0.13, 0.18), (0.18, 0.25), (0.25, 0.33), (0.33, 0.42), (0.42, 0.5)]
u = ks/N
print('N %d, %d matrices: mean spacing variance %.5f (sine process 0.17999)' % (N, blk*G, sV.sum()/(blk*G) - 1))
for lo, hi in bands:
    sel = (u >= lo) & (u < hi)
    if not sel.any(): continue
    w = cl[sel]**2*ks[sel]**2                      # weight: closure slope^2 * Var(P) (P exponential with mean k): inverse-variance for the ratio
    est = lambda sl_: float((sl_[sel]/cl[sel]*w).sum()/w.sum())
    a = est(sl); jk = np.array([est(slope(np.delete(np.arange(G), j))[0]) for j in range(G)])
    e = np.sqrt((G - 1)/G*((jk - jk.mean())**2).sum())
    print('  u in [%.2f, %.2f) (k = %d..%d): a = slope/closure = %.3f +- %.3f' % (lo, hi, ks[sel][0], ks[sel][-1], a, e), flush=True)
