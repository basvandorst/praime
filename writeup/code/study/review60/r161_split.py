"""Amplitude of the non-principal part: data minus the measured principal-form-prime part, against the prediction minus
its n^2/4 spikes (coefficient 1, 1/sqrt|y - n^2/4|, right side for D < 0, left side for D > 0), on (0.4,3]; jackknife over fields."""
import json, numpy as np
H = "/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/review_159_161"
pred = json.load(open("/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/cubic/predictions_cubic.json"))
E = np.array(pred['edges']); h = np.diff(E)
for sgn in (-1, 1):
    key = '%+d' % sgn; z = np.load(H + "/r161_all_%+d.npz" % sgn); P = np.array(pred['pred'][key])
    Kp = np.zeros(60)
    for n in range(1, 4):
        y0 = n * n / 4
        if sgn < 0: Kp += np.where(E[1:] > y0, (2 * np.sqrt(np.clip(E[1:] - y0, 0, None)) - 2 * np.sqrt(np.clip(E[:-1] - y0, 0, None))) / h, 0)
        else: Kp += np.where(E[:-1] < y0, (2 * np.sqrt(np.clip(y0 - E[:-1], 0, None)) - 2 * np.sqrt(np.clip(y0 - E[1:], 0, None))) / h, 0)
    sel = E[:-1] >= 0.4 - 1e-9
    num = z['NUM'] - z['NUMp'].sum(0); den = z['DEN']; Q = (P - Kp)
    n_ = len(num); N = num.sum(0); Dn = den.sum(0); m = N / Dn; mj = (N - num) / (Dn - den)
    sig = np.sqrt((n_ - 1) / n_ * ((mj - mj.mean(0)) ** 2).sum(0)); w = 1 / sig ** 2
    s = np.sum((Q * m * w)[sel]) / np.sum((Q * Q * w)[sel]); sj = (mj[:, sel] * (Q * w)[sel]).sum(1) / np.sum((Q * Q * w)[sel])
    print("D %s 0: non-principal part on (0.4,3]: amplitude %.3f +- %.3f; share of |pred| in principal spikes (rms over bins) %.2f" % (
        '<' if sgn < 0 else '>', s, np.sqrt((n_ - 1) / n_ * np.sum((sj - sj.mean()) ** 2)), np.sqrt(np.mean(Kp[sel] ** 2)) / np.sqrt(np.mean(P[sel] ** 2))))
print("two-component fit of the full data on (0.4,3]: m = s1 Kp + s2 (P - Kp), diag jackknife weights; and the projections")
for sgn in (-1, 1):
    key = '%+d' % sgn; z = np.load(H + "/r161_all_%+d.npz" % sgn); P = np.array(pred['pred'][key])
    Kp = np.zeros(60)
    for n in range(1, 4):
        y0 = n * n / 4
        if sgn < 0: Kp += np.where(E[1:] > y0, (2 * np.sqrt(np.clip(E[1:] - y0, 0, None)) - 2 * np.sqrt(np.clip(E[:-1] - y0, 0, None))) / h, 0)
        else: Kp += np.where(E[:-1] < y0, (2 * np.sqrt(np.clip(y0 - E[:-1], 0, None)) - 2 * np.sqrt(np.clip(y0 - E[1:], 0, None))) / h, 0)
    sel = E[:-1] >= 0.4 - 1e-9; num, den = z['NUM'], z['DEN']
    n_ = len(num); N = num.sum(0); Dn = den.sum(0); m = N / Dn; mj = (N - num) / (Dn - den)
    sig = np.sqrt((n_ - 1) / n_ * ((mj - mj.mean(0)) ** 2).sum(0)); w = (1 / sig ** 2)[sel]
    X = np.vstack([Kp[sel], (P - Kp)[sel]]).T
    A = np.linalg.inv(X.T @ (X * w[:, None])) @ (X * w[:, None]).T
    s = A @ m[sel]; sj = mj[:, sel] @ A.T; cs = (n_ - 1) / n_ * (sj - sj.mean(0)).T @ (sj - sj.mean(0))
    PP = np.sum(P[sel] ** 2 * w)
    print("  D %s 0: s1 = %.3f +- %.3f, s2 = %.3f +- %.3f, corr %.2f;  <P,Kp>/<P,P> = %.2f, <P,P-Kp>/<P,P> = %.2f" % (
        '<' if sgn < 0 else '>', s[0], np.sqrt(cs[0, 0]), s[1], np.sqrt(cs[1, 1]), cs[0, 1] / np.sqrt(cs[0, 0] * cs[1, 1]),
        np.sum(P[sel] * Kp[sel] * w) / PP, np.sum(P[sel] * (P - Kp)[sel] * w) / PP))
