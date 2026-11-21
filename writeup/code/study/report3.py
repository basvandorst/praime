"""Round 3 verdicts: block z for every statistic with a theory value; T2 = Pearson chi^2 of the 361-cell joint law."""
import json, os, math
import numpy as np
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')
th = json.load(open(os.path.join(D, 'theory_r3.json')))
rows = [json.loads(l) for l in open(os.path.join(D, 'results_r3.jsonl'))]
out = {}
for win in th:
    P = [r for r in rows if r['window'] == win]
    whole = [r for r in P if r['block'] == -1][0]
    blocks = sorted([r for r in P if r['block'] >= 0], key=lambda r: r['block'])
    out[win] = {}
    print('\n### window %s (%d primes)\n' % (win, whole['members']))
    print('| statistic | primes | theory | rel. diff | z (blocks) | verdict |\n|---|---|---|---|---|---|')
    for s, T in th[win]['-1'].items():
        if s.startswith('_') or s not in whole['stats']:
            continue
        x = whole['stats'][s]
        rb = np.array([b['stats'][s] - th[win][str(b['block'])][s] for b in blocks])
        se = rb.std(ddof=1) / math.sqrt(len(rb))
        z = float(rb.mean() / se) if se > 0 else float('nan')
        v = 'fits' if abs(z) <= 3 else 'misses'
        out[win][s] = {'x': x, 'T': T, 'rel': (x - T) / abs(T) if T else None, 'zT': z, 'verdict': v}
        print('| %s | %.6g | %.6g | %+.2e | %.1f | %s |' % (s, x, T, (x - T) / abs(T) if T else 0, z, v))
    # T2: chi^2 over the 361 cells of the joint law, whole window
    obs = np.array(whole['stats']['_cells'], float)
    J = np.array(th[win]['-1']['_joint'])[:19, :19].ravel()
    exp = J * whole['stats']['_npairs']
    big = exp >= 5
    o2 = np.append(obs[big], obs[~big].sum()); e2 = np.append(exp[big], exp[~big].sum())
    if e2[-1] < 5:
        o2, e2 = o2[:-1], e2[:-1]
    chi2 = float(np.sum((o2 - e2) ** 2 / e2)); dof = len(o2)
    zc = (chi2 / dof - 1) / math.sqrt(2 / dof)
    out[win]['T2.chi2_joint'] = {'chi2': chi2, 'dof': dof, 'chi2_per_dof': chi2 / dof, 'z': zc,
                                 'sum_obs': float(o2.sum()), 'sum_exp': float(e2.sum()),
                                 'verdict': 'fits' if abs(zc) <= 3 else 'misses'}
    print('| T2.chi2_joint | chi2/dof = %.3f (dof %d) | 1 | | %.1f | %s |' % (chi2 / dof, dof, zc,
                                                                            out[win]['T2.chi2_joint']['verdict']))
json.dump(out, open(os.path.join(D, 'verdicts_r3.json'), 'w'), indent=1)
