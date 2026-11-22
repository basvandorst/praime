"""Zero-sum prediction: max over u = ln y in [ln 100, ln X] of |sum_{|gamma|<T} e^{i gamma u}/(1/4+gamma^2)|."""
import json, math, sys
import numpy as np
Z = json.load(open(sys.argv[1])); S = json.load(open(sys.argv[2]))
X = float(sys.argv[3]) if len(sys.argv) > 3 else 1e8
u = np.arange(math.log(100), math.log(X), 0.0005)
print('%-14s %6s %8s %8s %8s %8s %8s %8s %8s' % ('chi', 'gam1', 'w1/C', 'wT/C', 'predT', 'pred/C', 'obs', 'obs/C', 'ln y*'))
rows = {}
for name, z in Z.items():
    g = np.array(z['zeros']); w = 1 / (0.25 + g ** 2)
    C = S[name]['C']
    s = np.exp(1j * np.outer(u, g)) @ w
    i = int(np.argmax(np.abs(s)))
    g1 = g[np.argmin(np.abs(g))]
    w1 = w.max() * (2 if z['real'] else 1)
    rows[name] = dict(gamma1=float(g1), w1_over_C=float(w1 / C), wT_over_C=float(w.sum() / C), pred=float(abs(s[i])),
                      pred_frac=float(abs(s[i]) / C), pred_lny=float(u[i]), obs=S[name]['max'], obs_frac=S[name]['frac'],
                      obs_lny=math.log(S[name]['at']))
    r = rows[name]
    print('%-14s %6.3f %8.3f %8.3f %8.4f %8.3f %8.4f %8.3f %5.2f/%5.2f' % (name, abs(g1), r['w1_over_C'], r['wT_over_C'], r['pred'], r['pred_frac'], r['obs'], r['obs_frac'], r['pred_lny'], r['obs_lny']))
json.dump(rows, open('predict.json', 'w'), indent=1)
