"""Estimator bias of I.migap_lag1 (Miller-Madow) and T4.H_cond (plug-in) at the sample sizes of the blocks:
draw N iid pairs from the theory's joint law J (multinomial), apply the same estimators, 400 repetitions.
Prints, per window: block N, true value, mean estimate, bias, and the bias-corrected block z."""
import json, os, math
import numpy as np
import probes as PR
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
th = json.load(open(os.path.join(D, 'theory_r3.json')))
rows = [json.loads(l) for l in open(os.path.join(D, 'results_r3.jsonl'))]
rng = np.random.default_rng(20261122)
def Hb(c):
    p = c / c.sum(); p = p[p > 0]; return float(-np.sum(p * np.log2(p)))
out = {}
for win in th:
    blocks = sorted([r for r in rows if r['window'] == win and r['block'] >= 0], key=lambda r: r['block'])
    res = {}
    for s in ('I.migap_lag1', 'T4.H_cond'):
        rb = []
        for b in blocks:
            J = np.clip(np.array(th[win][str(b['block'])]['_joint']), 0, None); J /= J.sum()
            N = b['stats']['_npairs']
            est = []
            for _ in range(40):
                cnt = rng.multinomial(N, J.ravel()).reshape(20, 20)
                if s == 'I.migap_lag1':
                    x = np.repeat(np.arange(20), cnt.sum(1)); 
                    # rebuild pair lists for the probes' estimator
                    xs, ys = np.nonzero(cnt); w = cnt[xs, ys]
                    est.append(PR._mi(np.repeat(xs, w), np.repeat(ys, w), 20, 20))
                else:
                    est.append(Hb(cnt.ravel().astype(float)) - Hb(cnt.sum(1).astype(float)))
            bias = float(np.mean(est)) - th[win][str(b['block'])][s]
            rb.append(b['stats'][s] - th[win][str(b['block'])][s] - bias)
        rb = np.array(rb); z = rb.mean() / (rb.std(ddof=1) / math.sqrt(len(rb)))
        res[s] = {'z_corrected': float(z), 'mean_resid_corrected': float(rb.mean())}
    out[win] = res
    print(win, {k: round(v['z_corrected'], 2) for k, v in res.items()}, flush=True)
json.dump(out, open(os.path.join(D, 'bias_r3.json'), 'w'), indent=1)
