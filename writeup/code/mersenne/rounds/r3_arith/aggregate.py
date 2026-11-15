"""Family-wise calibration over all scored features and pre-registration of the top 5 (validation only)."""
import glob, os, pickle, json, hashlib, time
import numpy as np
import harness as H

res = {}
for f in sorted(glob.glob(os.path.join(H.HERE, "results", "*.pkl"))):
    g = os.path.basename(f)[:-4]
    for nm, r in pickle.load(open(f, "rb")).items():
        res[nm] = (g, r)
names = list(res)
PV = np.array([res[n][1]["pv"] for n in names])          # (nf, NS+1)
BT = np.array([res[n][1]["bits"] for n in names])
nf, ns1 = PV.shape
print(f"features scored: {nf}; label sets: 1 real + {ns1-1} synthetic")
by_g = {}
for n in names:
    by_g[res[n][0]] = by_g.get(res[n][0], 0) + 1
print("by group:", by_g)

minp = PV.min(axis=0)
maxb = BT.max(axis=0)
real_min, syn_min = minp[0], minp[1:]
fw_p = (1 + np.sum(syn_min <= real_min)) / ns1
bar5 = np.quantile(syn_min, 0.05)
print("\nFamily-wise (best validation p over all features):")
print("  synthetic min-p quantiles 1%%/5%%/50%%: %.2e / %.2e / %.2e" % tuple(np.quantile(syn_min, [0.01, 0.05, 0.5])))
print("  real min p_v = %.3g  -> family-wise p = %.3f (share of null sets with a smaller best p)" % (real_min, fw_p))
print("  5%% family-wise bar for one feature: p_v < %.2e" % bar5)
real_maxb, syn_maxb = maxb[0], maxb[1:]
print("  best bits: real %.2f; null 50%%/95%% %.2f / %.2f; share of null sets >= real %.3f" %
      (real_maxb, np.quantile(syn_maxb, .5), np.quantile(syn_maxb, .95), np.mean(syn_maxb >= real_maxb)))
n05 = (PV < 0.05).sum(axis=0)
print("  features with p_v < 0.05: real %d; null median %d, 95%% %d; share null >= real %.3f" %
      (n05[0], np.median(n05[1:]), np.quantile(n05[1:], .95), np.mean(n05[1:] >= n05[0])))
sumb = BT.sum(axis=0)
print("  sum of bits over all features: real %.1f; null median %.1f, 95%% %.1f" % (sumb[0], np.median(sumb[1:]), np.quantile(sumb[1:], .95)))

order = np.argsort(PV[:, 0], kind="stable")
print("\nTop 15 features by validation p (real labels):")
print("%-30s %-8s %-16s %-10s %-11s %-9s %-7s %s" % ("feature", "group", "class", "D O/E", "V O/E", "p_v", "bits", "fw rank"))
for i in order[:15]:
    g, r = res[names[i]]
    c = r["classes"][int(r["cls"][0])]
    rank = np.mean(syn_min <= PV[i, 0])
    print("%-30s %-8s %-16s %2d/%-7.2f %2d/%-8.2f %-9.3g %-7.2f %.3f" % (names[i], g, H.describe(c), r["Od"][0], r["Ed"][0],
          r["Ov"][0], r["Ev"][0], PV[i, 0], BT[i, 0], rank))

prereg_path = os.path.join(H.HERE, "prereg.json")
if not os.path.exists(prereg_path):
    top = []
    for i in order[:5]:
        g, r = res[names[i]]
        c = r["classes"][int(r["cls"][0])]
        top.append({"feature": names[i], "group": g, "class": [int(c[0]), float(c[1])], "class_text": H.describe(c),
                    "D": [int(r["Od"][0]), float(r["Ed"][0]), float(r["pd"][0])],
                    "V": [int(r["Ov"][0]), float(r["Ev"][0]), float(r["pv"][0])], "bits": float(BT[i, 0])})
    doc = {"time": time.strftime("%Y-%m-%d %H:%M:%S"), "harness_sha256": hashlib.sha256(open(os.path.join(H.HERE, "harness.py"), "rb").read()).hexdigest(),
           "n_features": nf, "familywise_bar_5pct": float(bar5), "real_min_pv": float(real_min), "familywise_p": float(fw_p), "top5": top}
    json.dump(doc, open(prereg_path, "w"), indent=1)
    print("\nPRE-REGISTERED", prereg_path, hashlib.sha256(open(prereg_path, "rb").read()).hexdigest())
else:
    print("\nprereg.json exists already (not overwritten)")
