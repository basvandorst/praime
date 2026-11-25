"""Family periodogram: M(u) = mean over curves of sum_{0<g<=150} e^{2i theta(g)} e^{-iug} - leakage of the per-curve main term
(-eps c/2pi) int_0^150 e^{-iut} dt; |M(u)| against the noise from curve-level sign flips. Coherent lines r = e^u show as peaks.
Genus 2 (55 curves) and, for comparison, Sym^3 (on (0, 50], 75 curves)."""
import math, sys
import numpy as np
sys.path.insert(0, '.')
from r1_offsets import theta, load
rows = load()
us = np.arange(-4, 4.0001, 0.01)
rng = np.random.default_rng(3)
for cls, T in (('g2', 150.0), ('sym3', 50.0)):
    S = [r for r in rows if r['cls'] == cls]
    Z = []
    for r in S:
        g = r['g']; g = g[(g > 0) & (g <= T)]; ph = np.exp(2j * theta(g, r['Q'], r['mu']))
        F = (ph[None, :] * np.exp(-1j * us[:, None] * g[None, :])).sum(1)
        lk = np.where(np.abs(us) < 1e-12, T, (np.exp(-1j * us * T) - 1) / (-1j * np.where(np.abs(us) < 1e-12, 1, us)))
        Z.append(F + r['eps'] * r['c'] / (2 * math.pi) * lk)
    Z = np.array(Z); M = Z.mean(0)
    nulls = np.array([np.abs((Z * rng.choice([-1, 1], size=(len(S), 1))).mean(0)) for _ in range(300)])
    q99 = np.quantile(nulls, 0.99, axis=0); sd = np.sqrt((np.abs(Z - M) ** 2).sum(0) / (len(S) - 1) / len(S))
    z = np.abs(M) / sd
    print('%s: |M(u)|/se, top peaks (u, e^u, |M|, z, M):' % cls)
    idx = np.argsort(-z)
    shown = []
    for i in idx:
        if any(abs(us[i] - s) < 0.08 for s in shown): continue
        shown.append(us[i])
        print('   u %+.2f  r %.3f  |M| %.2f  z %.1f  M %+.2f%+.2fi  (null 99%% %.2f)' % (us[i], math.exp(us[i]), abs(M[i]), z[i], M[i].real, M[i].imag, q99[i]))
        if len(shown) >= 12: break
    print('   fraction of u-grid with z > 3: %.3f (Gaussian 0.011 for |complex| z>3: %.4f)' % (np.mean(z > 3), math.exp(-4.5)))
    np.savez('r13_%s.npz' % cls, us=us, M=M, sd=sd)
