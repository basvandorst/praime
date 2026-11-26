# Reads the scans of round 165 and prints the curve t*(r), the feasible set and its components, the local maxima.
# Usage: report.py <jsonl> [tau]    (the answer key is not read here)
import sys, json, numpy as np
rows = sorted((json.loads(l) for l in open(sys.argv[1])), key=lambda d: d['r'])
tau = float(sys.argv[2]) if len(sys.argv) > 2 else 0.002
r = np.array([d['r'] for d in rows]); t = np.array([d['t_lp'] for d in rows])
print('%d points, r %.2f..%.2f, max |t_lp - t_low| %.1e, total %.0f s' % (len(r), r[0], r[-1], max(abs(d['t_lp'] - d['t_low']) for d in rows), sum(d['secs'] for d in rows)))
feas = t >= -tau; comp = []; i = 0
while i < len(r):
    if feas[i]:
        j = i
        while j + 1 < len(r) and feas[j + 1]: j += 1
        k = i + int(np.argmax(t[i:j + 1])); comp.append((r[i], r[j], r[k], t[k])); i = j + 1
    else: i += 1
print('feasible components (t* >= -%g): ' % tau + '; '.join('[%.2f, %.2f] max %.5f at %.2f' % (a, b, tk, rk) for a, b, rk, tk in comp))
mx = [(r[i], t[i]) for i in range(1, len(r) - 1) if t[i] > t[i - 1] and t[i] >= t[i + 1]]
print('local maxima: ' + '; '.join('%.2f (%.4f)' % m for m in mx))
step = max(1, len(r)//76)
print('curve (every %d-th point): ' % step + ' '.join('%.1f:%.3f' % (a, b) for a, b in zip(r[::step], t[::step])))
