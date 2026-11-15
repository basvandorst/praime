"""Tests of the frozen harness: Poisson-binomial vs a plain DP, and planted features."""
import numpy as np
import harness as H

rng = np.random.default_rng(1)
# 1. Poisson-binomial: harness discovery PMF vs direct DP on a random class
v = rng.random(len(H.P))
cls = (0, 0.37)
vD, wD = v[H.DMASK], H.W[H.DMASK]
pmf, E = H.disc_pmf(vD, wD, [cls])
m = H.in_class(vD, cls)
d = np.zeros(H.K + 1); d[0] = 1
for wi in wD[m]:
    d[1:] = d[1:] * (1 - wi) + d[:-1] * wi
    d[0] *= (1 - wi)
print("PB discovery max abs diff vs DP:", np.abs(pmf[0] - d).max(), " E", E[0], wD[m].sum())
vV, wV = v[H.VMASK], H.W[H.VMASK]
st = H.class_stats(vV, wV, [cls])
pv = H.pmf_from_sums(st[:, 1], st[:, 2:])[0]
mv = H.in_class(vV, cls)
d = np.zeros(H.K + 1); d[0] = 1
for wi in wV[mv]:
    d[1:] = d[1:] * (1 - wi) + d[:-1] * wi
    d[0] *= (1 - wi)
print("PB validation max abs diff vs DP:", np.abs(pv - d).max(), " tail P(X>=20):", d[20:].sum(), H.survival(pv[None])[0, 20])

# 2. planted feature: h = frac(k * 0.75487766624669), class h <= 0.3 has rate x m
h = np.mod(H.P * 0.75487766624669276, 1.0)
inC = h <= 0.3
for mult in (2.0, 3.0, 5.0):
    sets = [H.REAL]
    for i in range(200):
        wp = np.minimum(H.W * np.where(inC, mult, 1.0), 1.0)
        sets.append(np.flatnonzero(rng.random(len(H.P)) < wp))
    r = H.evaluate(h, sets)
    pv = r["pv"][1:]
    print(f"planted x{mult}: chosen class == planted-side in {np.mean([H.describe(r['classes'][int(c)])[:2]=='<=' for c in r['cls'][1:]]):.2f}; "
          f"median p_v {np.median(pv):.2e}; share p_v<0.05 {np.mean(pv<0.05):.2f}, <1e-3 {np.mean(pv<1e-3):.2f}, <1e-4 {np.mean(pv<1e-4):.2f}; median bits {np.median(r['bits'][1:]):.2f}")

# 3. planted on the real labels: a feature that holds on every exponent and on 30% of the rest
f = (np.isin(H.P, H.EXPS) | inC).astype(float)
r = H.evaluate(f, [H.REAL])
print("planted on real labels: class", H.describe(r["classes"][int(r["cls"][0])]), "Od/Ed %d/%.2f Ov/Ev %d/%.2f p_v %.2e bits %.1f" %
      (r["Od"][0], r["Ed"][0], r["Ov"][0], r["Ev"][0], r["pv"][0], r["bits"][0]))
# 4. a null feature on the real labels
r = H.evaluate(h, [H.REAL])
print("null feature on real labels: Od/Ed %d/%.2f Ov/Ev %d/%.2f p_v %.3f" % (r["Od"][0], r["Ed"][0], r["Ov"][0], r["Ev"][0], r["pv"][0]))
