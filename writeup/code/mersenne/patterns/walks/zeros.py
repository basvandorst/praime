"""Explicit-formula reconstruction of the three walks from zeros (heuristic check of the explanation).
  M(x)/sqrt x  ~ 2 Re sum_{g>0} x^{ig} / (rho zeta'(rho))                         (zeta zeros)
  L(x)/sqrt x  ~ 1/zeta(1/2) + 2 Re sum_{g>0} zeta(2 rho) x^{ig} / (rho zeta'(rho))
  R(x) log x / sqrt x ~ 1 + 2 Re sum_{g>0} x^{ig} / rho                          (zeros of L(s,chi_-4))
Truncated at K zeros; compared with the exact series on 3000 log-spaced x in [1000, 500000]."""
import json, math
import numpy as np
import mpmath as mp
from common import series, N

mp.mp.dps = 20
K = 150
zz = []
for k in range(1, K + 1):
    r = mp.zetazero(k)
    zp = mp.zeta(r, derivative=1)
    z2 = mp.zeta(2 * r)
    zz.append((float(r.imag), complex(1 / (r * zp)), complex(z2 / (r * zp))))
g = np.array([a for a, _, _ in zz]); cM = np.array([b for _, b, _ in zz]); cL = np.array([c for _, _, c in zz])

# zeros of L(s, chi_-4): sign changes of the real completed function on the critical line
chi = [0, 1, 0, -1]
def Lam(t):
    s = mp.mpf(0.5) + 1j * mp.mpf(t)
    return (mp.mpf(4) / mp.pi) ** ((s + 1) / 2) * mp.gamma((s + 1) / 2) * mp.dirichlet(s, chi)
T = 120.0
ts = np.arange(0.5, T, 0.05)
vals = [float(mp.re(Lam(t))) for t in ts]
imag_max = max(abs(float(mp.im(Lam(t))) / (abs(float(mp.re(Lam(t)))) + 1e-300)) for t in ts[::200])
gch = []
for i in range(len(ts) - 1):
    if vals[i] == 0 or vals[i] * vals[i + 1] < 0:
        a, b = ts[i], ts[i + 1]
        fa = vals[i]
        for _ in range(40):
            m = (a + b) / 2
            fm = float(mp.re(Lam(m)))
            if fa * fm <= 0:
                b = m
            else:
                a, fa = m, fm
        gch.append((a + b) / 2)
gch = np.array(gch)
cR = 1 / (0.5 + 1j * gch)

S = series(N)
xs = np.unique(np.geomspace(1000, N, 3000).astype(np.int64))
lx = np.log(xs.astype(float))
def recon(coef, gam, base):
    return base + 2 * np.real(np.exp(1j * np.outer(lx, gam)) @ coef)

out = {"K_zeta": K, "gamma_K": float(g[-1]), "chi4_zeros_found": len(gch), "chi4_first": [round(x, 6) for x in gch[:8]],
       "chi4_imag_rel_check": imag_max}
actual = {"M": S["M"][xs] / np.sqrt(xs), "L": S["L"][xs] / np.sqrt(xs), "R": S["R"][xs] * lx / np.sqrt(xs)}
for key, coef, gam, base in [("M", cM, g, 0.0), ("L", cL, g, -0.6847652800500), ("R", cR, gch, 1.0)]:
    a = actual[key]
    for k in [10, 50, len(gam)]:
        f = recon(coef[:k], gam[:k], base)
        if key == "M":
            f = f - 2 / np.sqrt(xs)
        if key == "R":
            # finite-x mean: (1/2) li(sqrt x) * log x / sqrt x instead of 1
            f = f - 1 + np.array([float(mp.li(math.sqrt(x))) / 2 for x in xs]) * lx / np.sqrt(xs)
        corr = float(np.corrcoef(a, f)[0, 1])
        out[f"{key}_K{k}"] = dict(corr=corr, rms_resid=float(np.sqrt(np.mean((a - f) ** 2))),
                                  rms_actual_about_mean=float(np.std(a)))
    # predicted log-scale sd (all K zeros): sqrt(2 sum |c|^2)
    out[f"{key}_pred_sd_trunc"] = float(np.sqrt(2 * np.sum(np.abs(coef) ** 2)))
json.dump(out, open("zeros_out.json", "w"), indent=1)
np.savez("zeros_recon.npz", xs=xs, **actual)
for k, v in out.items():
    print(k, v)
