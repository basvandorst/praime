"""Review 26: round 90 checks. Own law (law_win.npz, 17 lam on [ln1e7, ln1e8]) against r90.json and own counts."""
import json, math
import numpy as np
from scipy.interpolate import CubicSpline
J = json.load(open('/home/bas/praime/writeup/data/study/r90.json'))
C = json.load(open('/home/bas/scratch/study-review26/count.json'))
Z = np.load('/home/bas/scratch/study-review26/law_win.npz'); dist = Z['dist']; lams = Z['lams']; Ds = Z['D']
QS = (3, 5, 7, 11, 13); KS = (1, 2, 3, 4)


def rep(dk):
    return {q: dk[Ds % q == 0].sum() / dk.sum() for q in QS}


def zs(th):
    out = []
    for k in KS:
        n = C['data'][str(k)]['n']
        for q in QS:
            m = C['data'][str(k)][str(q)] / n; se = math.sqrt(m * (1 - m) / n)
            out.append((th[k][q] - m) / se)
    return np.array(out)


# 0. data identical?
for k in KS:
    n = C['data'][str(k)]['n']
    for q in QS:
        assert abs(C['data'][str(k)][str(q)] / n - J['data'][str(k)][str(q)]['p']) < 1e-15, (k, q)
print('counts: identical to r90.json for all 20 values; primes', C['primes'], J['primes'])
# 1. own law at the three lam vs r90 theory_mid_lam
i3 = [0, 8, 16]
mx = 0
for k in KS:
    own = rep(dist[8, k - 1])
    for q in QS:
        mx = max(mx, abs(own[q] - J['theory_mid_lam'][str(k)][str(q)]) / J['data'][str(k)][str(q)]['se'])
print('own law at mid lam vs r90 theory_mid_lam: max |diff|/se = %.3f' % mx)
# 2. window averages
u = np.linspace(lams[0], lams[-1], 4001); wt = np.exp(u) / u
def avg_lin3():
    return {k: rep(np.array([np.sum(np.interp(u, lams[i3], dist[i3, k - 1, j]) * wt) / wt.sum() for j in range(len(Ds))])) for k in KS}
def avg_spline():
    out = {}
    for k in KS:
        cs = CubicSpline(lams, dist[:, k - 1, :], axis=0)
        out[k] = rep((cs(u) * wt[:, None]).sum(0) / wt.sum())
    return out
A3 = avg_lin3(); AS = avg_spline()
mx = 0; mx2 = 0
for k in KS:
    for q in QS:
        se = J['data'][str(k)][str(q)]['se']
        mx = max(mx, abs(A3[k][q] - J['theory'][str(k)][str(q)]) / se)
        mx2 = max(mx2, abs(A3[k][q] - AS[k][q]) / se)
print('own 3-lam linear average vs r90 theory: max |diff|/se = %.3f' % mx)
print('3-lam linear vs 17-lam spline average: max |diff|/se = %.3f' % mx2)
for name, th in (('r90 theory', {k: {q: J['theory'][str(k)][str(q)] for q in QS} for k in KS}), ('own lin3', A3), ('own spline', AS)):
    z = zs(th); print('%-12s max|z| %.2f  sum z2 %.1f' % (name, np.abs(z).max(), (z * z).sum()))
# 3. single lam
for lam in (lams[0], lams[8], 17.656, lams[-1]):
    th = {}
    for k in KS:
        cs = CubicSpline(lams, dist[:, k - 1, :], axis=0); th[k] = rep(cs(lam))
    z = zs(th); print('single lam %.3f: max|z| %.2f sum z2 %.1f  (k,q of max: %s)' % (lam, np.abs(z).max(), (z * z).sum(), divmod(int(np.abs(z).argmax()), 5)))
# weighted mean of lam and of 1/lam
print('prime-weighted mean lam %.3f, 1/mean(1/lam) %.3f' % ((u * wt).sum() / wt.sum(), wt.sum() / (wt / u).sum()))
# per-value z table, spline average
z = zs(AS); print('spline-avg z:', np.round(z, 2).reshape(4, 5))
