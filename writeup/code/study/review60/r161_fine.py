"""fine-bin (width 0.05) shape test on (0.4,3]: amplitude fitted, full jackknife covariance over fields; also the 0.6364 constant."""
import json, numpy as np
H = "/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/review_159_161"
pred = json.load(open("/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/cubic/predictions_cubic.json"))
E = np.array(pred['edges'])
for sgn in (-1, 1):
    key = '%+d' % sgn; z = np.load(H + "/r161_all_%+d.npz" % sgn); P = np.array(pred['pred'][key])
    sel = E[:-1] >= 0.4 - 1e-9
    for nm, (nu_, de_, PP) in {"fine 52": (z['NUM'][:, sel], z['DEN'][:, sel], P[sel]),
                               "fine 52 w/o principal-form primes": (z['NUMx'][:, sel], z['DENx'][:, sel], None)}.items():
        n = len(nu_); N = nu_.sum(0); Dn = de_.sum(0); m = N / Dn; mj = (N - nu_) / (Dn - de_)
        dm = mj - mj.mean(0); cov = (n - 1) / n * dm.T @ dm; sig = np.sqrt(np.diag(cov))
        if PP is None:
            print("  %s: mean %.3f (diag-weighted), chi2 vs own mean %.1f" % (nm, np.sum(m / sig**2) / np.sum(1 / sig**2), np.sum(((m - np.sum(m / sig**2) / np.sum(1 / sig**2)) / sig)**2)))
            continue
        Ci = np.linalg.inv(cov); s = (PP @ Ci @ m) / (PP @ Ci @ PP); r = m - s * PP
        rd = (m - s * PP) / sig
        print("D %s 0 %s: s(GLS) %.3f; chi2 full cov %.1f, diag %.1f (52 bins, 51 dof); bins with |resid| > 2.5 sigma: %s" % (
            '<' if sgn < 0 else '>', nm, s, r @ Ci @ r, np.sum(rd**2), [round(float(E[:-1][sel][i]), 2) for i in np.nonzero(np.abs(rd) > 2.5)[0]]))
def primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0].astype(float)
l = primes(10 ** 7)
print("prod (1 - l^-2 + l^-5) = %.6f" % np.exp(np.sum(np.log1p(-l ** -2 + l ** -5))))
