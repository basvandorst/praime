"""disc4.py -- discovery: (a) the low band [2, 300] (long scales, where primes are most rigid) for all patterns;
(b) scrambled control: Poisson pseudo-pairs with the same smooth mean and the measured variance -> all coefficients 0;
(c) the transfer function primes -> twins in bands of gamma (cross-spectrum), the 'spectrum' view."""
import numpy as np
import tzlib, theory, models

D = models.Data(["disc.npz"])
j0, j1 = D.J(10 ** 6), D.J(2 * 10 ** 9, right=True) - 1

print("(a) low band [2, 300]: zeta coefficient and joint chi loudness")
for name in ("h2", "h4", "h6", "h8", "h10", "h12", "h30", "SG", "T026", "T046"):
    Y, S = D.Ypat(name)
    S_, regs = models.build(D, name, null_set=[])
    b, cov, d = tzlib.fit(Y, [r[1] for r in regs], D.u, j0, j1, 2, 300, models.kforms(name))
    e = np.sqrt(np.diag(cov))
    print("  %-5s " % name + "   ".join("%s %.3f+-%.3f (pred %.3f)" % (r[0].split(" (")[0], x, s, r[2]) for r, x, s in zip(regs, b, e)))

print("\n(b) scrambled control: Poisson pseudo-twins, mean S/ln^2 per integer, weight ln^2, 5 seeds")
rng = np.random.default_rng(12345)
S = theory.singular_const("h2")
Xf = D.Xf
width = np.diff(np.concatenate([[Xf[0] - 1], Xf]))
L = D.u
S_, regs = models.build(D, "h2")
for seed in range(5):
    lam = S * width / L ** 2
    cnt = rng.poisson(lam)
    cellw = cnt * L ** 2
    cumw = np.cumsum(cellw)
    Yp = (cumw - S * Xf) / D.sq
    b, cov, d = tzlib.fit(Yp, [r[1] for r in regs], D.u, j0, j1, 2, 6000, 2)
    e = np.sqrt(np.diag(cov))
    print("  seed %d: zeta %.3f+-%.3f  chi-loudness %.3f+-%.3f  null-set max|z| %.2f" % (seed, b[0], e[0], b[1], e[1], np.max(np.abs(b[2:] / e[2:]))))

print("\n(c) transfer function H(gamma) = <F_pair conj F_theta> / <|F_theta|^2> in bands, twins h=2 after removing the predicted chi part")
Y, S = D.Ypat("h2")
S_, regs = models.build(D, "h2", null_set=[])
ychi = regs[1][1]
Yres = Y - ychi                               # subtract the predicted character part (loudness 1)
uu = D.u[j0:j1 + 1]; w = uu[:-1] ** -1.0
Dp = np.fft.rfft(np.diff(Yres[j0:j1 + 1]) * w)
Dt = np.fft.rfft(np.diff(D.Yzeta()[j0:j1 + 1]) * w)
gam = 2 * np.pi * np.fft.rfftfreq(len(uu) - 1, d=tzlib.DU)
edges = [2, 30, 60, 100, 200, 400, 800, 1500, 3000, 6000]
for a, bnd in zip(edges[:-1], edges[1:]):
    m = (gam >= a) & (gam < bnd)
    num = np.sum(Dp[m] * np.conj(Dt[m])); den = np.sum(np.abs(Dt[m]) ** 2)
    H = num / den
    # error: residual power / regressor power
    res = Dp[m] - H * Dt[m]
    err = np.sqrt(np.sum(np.abs(res) ** 2) / (m.sum() - 1) / den / 2)
    coh = np.abs(num) ** 2 / den / np.sum(np.abs(Dp[m]) ** 2)
    print("  gamma [%4d,%4d): bins %5d  H = %6.3f %+6.3fi  +- %.3f   (pred 2S = %.3f)  coherence^2 %.3f" % (a, bnd, m.sum(), H.real, H.imag, err, 2 * S, coh))
