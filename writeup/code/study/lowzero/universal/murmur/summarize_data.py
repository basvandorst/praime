"""Summarize data_prime.npz and data_squarefree.npz (from trace_density_data.py)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from densities import H_int, GAMMA

D = os.path.dirname(os.path.abspath(__file__))
import glob
fam = {}
for k in ("prime", "squarefree"):
    acc = None
    for fn in sorted(glob.glob(os.path.join(D, f"data_{k}_*.npz"))):
        z = dict(np.load(fn))
        print("chunk", os.path.basename(fn), "levels", int(z["nlev"]), "last N", int(z["last"]))
        if acc is None:
            acc = z
        else:
            for key in ("T", "dim", "zub", "prime", "T2", "cnt", "trans", "wsum", "nlev"):
                acc[key] = acc[key] + z[key]
    fam[k] = acc
ref = {(10.0, 2.0): {"prime": 0.582042 - 0.007767j, "squarefree": 0.841700 - 0.009061j},
       (14.0, 2.0): {"prime": 0.639320 + 0.007696j, "squarefree": 0.894776 + 0.009868j}}
main = {"prime": 6 / np.pi ** 2, "squarefree": GAMMA / (2 * np.pi ** 2)}

e = fam["prime"]["edges"]
mid = 0.5 * (e[1:] + e[:-1])
print("Binned natural averages E[eps_f a_f(p)] = sum T / sum dim  (bin width 0.25 in y = p/N); err = standard error of the bin mean")
print("  y_mid | PRIME levels: data   err    M_prime  M_zub | SQUAREFREE levels: data   err    M_zub  M_prime")
stats = {}
for k, d in fam.items():
    mean = d["T"] / d["dim"]
    var = d["T2"] / d["dim"] - mean ** 2           # dim-weighted variance of T/dim per sample
    neff = d["cnt"]
    err = np.sqrt(np.maximum(var, 0) / np.maximum(neff, 1))
    stats[k] = (mean, err, d["prime"] / d["dim"], d["zub"] / d["dim"])
for i in range(len(mid)):
    if i % 2 == 1 and mid[i] > 3:
        continue
    a, b = stats["prime"], stats["squarefree"]
    print(f"  {mid[i]:5.3f} | {a[0][i]:13.3f} {a[1][i]:6.3f} {a[2][i]:8.3f} {a[3][i]:6.3f} | {b[0][i]:17.3f} {b[1][i]:6.3f} {b[3][i]:6.3f} {b[2][i]:8.3f}")
print()
for k in ("prime", "squarefree"):
    mean, err, mp_, mz = stats[k]
    sel = err > 0
    chi_p = np.sum(((mean - mp_) / err)[sel] ** 2)
    chi_z = np.sum(((mean - mz) / err)[sel] ** 2)
    print(f"{k:10s}: {int(fam[k]['nlev'])} levels, total dim {fam[k]['wsum']:.0f}; bins {sel.sum()}:"
          f"  chi2 vs M_prime = {chi_p:.1f},  chi2 vs M_zub = {chi_z:.1f};"
          f"  RMS(data-M_prime) = {np.sqrt(np.mean((mean - mp_) ** 2)):.3f}, RMS(data-M_zub) = {np.sqrt(np.mean((mean - mz) ** 2)):.3f}")
print("\nTransform of the data = natural-weight family average of the prime part of sum_gamma H(gamma) e^{2 i theta_f(gamma)} (times eps_f),")
print("divided by -int H/2pi; compare with the transform of the density and with the main term (r = 1 line):")
for j, w in enumerate(fam["prime"]["windows"]):
    w = tuple(float(v) for v in w)
    norm = -H_int(*w) / (2 * np.pi)
    for k in ("prime", "squarefree"):
        val = fam[k]["trans"][j] / fam[k]["wsum"] / norm
        print(f"   t0={w[0]:4.1f} sigma={w[1]}: {k:10s} data {val.real:+.4f}{val.imag:+.4f}i   density transform {ref[w][k].real:+.4f}{ref[w][k].imag:+.4f}i"
              f"   main term {main[k]:.4f}")

# Noise estimate for the transform: binned version  -(1/2pi) sum_bins mean_bin * int_bin G(y)/y dy, with the bin standard errors.
from densities import G_kernel
print("\nBinned transform (same as above but from bin means) with 1-sigma noise from the bin standard errors:")
for j, w in enumerate(fam["prime"]["windows"]):
    w = tuple(float(v) for v in w)
    norm = -H_int(*w) / (2 * np.pi)
    for k in ("prime", "squarefree"):
        mean, err, mp_, mz = stats[k]
        Ib = []
        for i in range(len(mid)):
            yy = np.linspace(max(e[i], 1e-4), e[i + 1], 400)
            g = G_kernel(yy, *w, nt=1201)
            Ib.append(np.trapezoid(g / yy, yy))
        Ib = np.array(Ib)
        val = -(1 / (2 * np.pi)) * np.sum(mean * Ib) / norm
        mod = -(1 / (2 * np.pi)) * np.sum((mp_ if k == "prime" else mz) * Ib) / norm
        sig = (1 / (2 * np.pi)) * np.sqrt(np.sum((err * np.abs(Ib)) ** 2)) / abs(norm)
        print(f"   t0={w[0]:4.1f}: {k:10s} binned data {val.real:+.4f}{val.imag:+.4f}i  +- {sig:.4f}   binned model {mod.real:+.4f}{mod.imag:+.4f}i")
