"""power.py: expected values and standard errors of the registered readings on the reserved design, from the frozen
predictions and a noise model calibrated on the development set ONLY (dev_out.json).

Noise model (multiplicative, log scale): for a statistic s_ijk (raw: ln(M/P_full); CV: ln(M/P_full) - z, z = sum of the
three pair log-ratios), s_ijk = a_i + a_j + a_k + e_ijk with character effects a ~ N(0, sig_a^2) and e ~ N(0, sig_e^2);
(sig_a, sig_e) fitted on the 56 dev triples by least squares with character dummies (method of moments for sig_a).
Pairs: the same with two characters.  Hypothesis H_c: ln M = ln P_c + noise.  1000 simulated data sets per hypothesis;
reported: mean reading, sd of the reading over simulations (= expected se), and the sd for sig_a doubled.
usage: power.py PRED.json DEV_OUT.json
"""
import sys, json, math, itertools
import numpy as np


def comps(y, tuples, chars):
    X = np.array([[1.0] + [1.0 if c in t else 0.0 for c in chars] for t in tuples])
    b = np.linalg.lstsq(X, y, rcond=None)[0]; res = y - X @ b
    dof = len(y) - np.linalg.matrix_rank(X); se2 = float((res ** 2).sum() / dof)
    a = b[1:] - b[1:].mean(); noise = se2 * float(np.mean(np.diag(np.linalg.pinv(X.T @ X))[1:]))
    return math.sqrt(se2), math.sqrt(max(float(a.var(ddof=1)) - noise, 0.0))


def calibrate(dev):
    DEV = [1001, -1003, 1005, -1048, 1068, -1195, 1253, -1267]
    tri = [tuple(r['Ds']) for r in dev['triple']]
    rp = {tuple(r['Ds']): r['meas'] / r['pred'] for r in dev['pair']}
    y = np.log([r['meas'] / r['pred']['full'] for r in dev['triple']])
    z = np.array([math.log(rp[(a, b)]) + math.log(rp[(a, c)]) + math.log(rp[(b, c)]) for a, b, c in tri])
    out = {'raw': comps(y, tri, DEV), 'cv': comps(y - z, tri, DEV)}
    yp = np.log(list(rp.values())); out['pair'] = comps(yp, list(rp.keys()), DEV)
    return out


def ols_slope(x, y):
    x = x - x.mean(); y = y - y.mean()
    return float(np.sum(x * y) / np.sum(x * x))


def main(argv):
    P = json.load(open(argv[1])); dev = json.load(open(argv[2]))
    cal = calibrate(dev)
    print('noise calibrated on dev (sig_e, sig_a): raw %.3f %.3f | CV %.3f %.3f | pairs %.3f %.3f' % (cal['raw'] + cal['cv'] + cal['pair']))
    chars = P['sets']['reserved']; ci = {c: i for i, c in enumerate(chars)}
    tri = [tuple(r['key']) for r in P['triples']]; TT = {tuple(r['key']): r for r in P['triples']}
    I = np.array([[ci[a], ci[b], ci[c]] for a, b, c in tri])
    comps_ = ['full', 'diag', 'kirk', 'no3', 'ind', 'lead']
    LP = {c: np.log([TT[t][c] for t in tri]) for c in comps_}
    print('reserved: %d characters, %d triples, %d pairs' % (len(chars), len(tri), len(P['pairs'])))
    print('ln(P_c/P_full) over the reserved triples: ' + '  '.join('%s mean %+.3f sd %.3f' % (c, np.mean(LP[c] - LP['full']), np.std(LP[c] - LP['full'])) for c in comps_[1:]))
    a3 = np.array([TT[t]['A3_over_Aij'] for t in tri]); print('three-body factor A3/prod Aij: median %.3f, 5-95%% %.3f .. %.3f, min %.3f max %.3f' % (np.median(a3), *np.percentile(a3, [5, 95]), a3.min(), a3.max()))
    rng = np.random.default_rng(1125)
    NS = 1000
    for label, (se_, sa_) in [('CV', cal['cv']), ('raw', cal['raw'])]:
        for mult in (1.0, 2.0):
            sa = sa_ * mult
            print('\n%s statistic, sig_e %.3f, sig_a %.3f%s' % (label, se_, sa, '  (sig_a doubled: sensitivity)' if mult > 1 else ''))
            print('   reading                         ' + ''.join('%12s' % ('H_' + c) for c in comps_) + '   sd(H_full)')
            res = {c: {'level': [], 'slope': [], **{'g_' + d: [] for d in comps_[1:]}} for c in comps_}
            for s in range(NS):
                a = rng.normal(0, sa, len(chars)); e = rng.normal(0, se_, len(tri))
                noise = a[I].sum(axis=1) + e
                noise = noise - (sa ** 2 * 3 + se_ ** 2) / 2          # E[exp(noise)] = 1: unbiased ratios
                for h in comps_:
                    y = LP[h] + noise
                    res[h]['level'].append(np.mean(np.exp(y - LP['full'])))
                    res[h]['slope'].append(ols_slope(LP['full'], y))
                    for d in comps_[1:]:
                        res[h]['g_' + d].append(ols_slope(LP['full'] - LP[d], y - LP[d]))
            for key in ['level', 'slope'] + ['g_' + d for d in comps_[1:]]:
                print('   %-30s' % key + ''.join('%12.4f' % np.mean(res[h][key]) for h in comps_) + '   %.4f' % np.std(res['full'][key]))
    # does the delete-one-character jackknife (used by moment_test.py) estimate the sampling sd? (CV noise, H_full)
    se_, sa_ = cal['cv']; jl, jg, vl, vg = [], [], [], []
    for s in range(40):
        a = rng.normal(0, sa_, len(chars)); e = rng.normal(0, se_, len(tri))
        y = LP['full'] + a[I].sum(axis=1) + e - (3 * sa_ ** 2 + se_ ** 2) / 2
        r = np.exp(y - LP['full']); u = LP['full'] - LP['no3']; w = y - LP['no3']
        vl.append(r.mean()); vg.append(ols_slope(u, w))
        rl, rg = [], []
        for ch in range(len(chars)):
            k = ~np.any(I == ch, axis=1)
            rl.append(r[k].mean()); rg.append(ols_slope(u[k], w[k]))
        n = len(chars)
        jl.append(math.sqrt((n - 1) / n * np.sum((np.array(rl) - np.mean(rl)) ** 2)))
        jg.append(math.sqrt((n - 1) / n * np.sum((np.array(rg) - np.mean(rg)) ** 2)))
    print('\njackknife check (40 simulated CV data sets under H_full): level: sd over sims %.4f, mean jackknife se %.4f;'
          ' gamma_no3: sd %.4f, mean jackknife se %.4f' % (np.std(vl), np.mean(jl), np.std(vg), np.mean(jg)))
    # pairs: level reading
    se_, sa_ = cal['pair']; pairs = [tuple(r['key']) for r in P['pairs']]
    Ip = np.array([[ci[a], ci[b]] for a, b in pairs]); lv = []
    for s in range(NS):
        a = rng.normal(0, sa_, len(chars)); e = rng.normal(0, se_, len(pairs))
        n = a[Ip].sum(axis=1) + e - (2 * sa_ ** 2 + se_ ** 2) / 2
        lv.append(np.mean(np.exp(n)))
    print('\npairs: pooled ratio expected 1, se %.4f (sig_e %.3f, sig_a %.3f)' % (np.std(lv), se_, sa_))
    PPd = {tuple(r['key']): r for r in P['pairs']}
    print('pairs: mean P_diag/P_full %.4f (sd %.4f); mean P_lead/P_full %.4f' % (np.mean([PPd[p]['diag'] / PPd[p]['full'] for p in pairs]),
          np.std([PPd[p]['diag'] / PPd[p]['full'] for p in pairs]), np.mean([PPd[p]['lead'] / PPd[p]['full'] for p in pairs])))
    zl = [r['zetaL'] for r in P['singles']]
    print('zeta*L: mean P_diag/P_full %.4f (range %.4f .. %.4f)' % (np.mean([z['diag'] / z['full'] for z in zl]), min(z['diag'] / z['full'] for z in zl), max(z['diag'] / z['full'] for z in zl)))


if __name__ == '__main__':
    main(sys.argv)
