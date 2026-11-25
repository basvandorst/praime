"""Expected values and errors of the registered readings (noise model from the used set R120 only; frozen predictions).
Regression: y_i = stat_i - ignored_i on [1, I(2|D), ..., I(13|D)]; predictions treated the same way (regression-matched).
Noise per L-function: sd = sd_R120 * sqrt(744/n_i) (n_i = expected spacings), Var and <r~> correlated (corr from R120)."""
import json, math, numpy as np
P = json.load(open('predictions_spacing.json')); NZ = P['noise']; PR = P['primes']
ALT = ('law', 'law2', 'slow', 'bblm', 'ignored'); out = {}
def design(D): return np.array([[1.0] + [float(abs(d) % p == 0) for p in PR] for d in D])
print('per set: predicted per-prime coefficients (p = 2, 3, 5; y = stat - ignored), expected se, expected gamma under each variant (gamma = amplitude on the law pattern, weights 1/se^2)')
for name in ('B', 'Dn', 'Fn', 'C', 'Bp', 'R120'):
    S = P[name]; X = design(S['D']); n = np.array(S['n_spacings_expected']); A = np.linalg.inv(X.T @ X)
    res = {}
    for st, sd0 in (('var', NZ['sd_var']), ('r', NZ['sd_r'])):
        sd = sd0 * np.sqrt(NZ['n_ref'] / n); Hm = A @ X.T; cov = (Hm * sd ** 2) @ Hm.T      # heteroscedastic OLS covariance
        ign = np.array(S['pred']['ignored'][st])
        coef = {v: (Hm @ (np.array(S['pred'][v][st]) - ign)) for v in ALT}
        se = np.sqrt(np.diag(cov))[1:4]; c = coef['law'][1:4]; W = np.diag(1 / se ** 2); C3 = cov[1:4, 1:4]
        den = c @ W @ c; seg = math.sqrt(c @ W @ C3 @ W @ c) / den
        gam = {v: float(c @ W @ coef[v][1:4] / den) for v in ALT}
        chi = {v: float((coef[v][1:4] - c) @ np.linalg.inv(C3) @ (coef[v][1:4] - c)) for v in ALT}
        lev_se = float(np.sqrt(np.mean(sd ** 2) / len(sd)))
        res[st] = dict(coef={v: coef[v][1:7].tolist() for v in ALT}, se=np.sqrt(np.diag(cov))[1:7].tolist(), gamma=gam, gamma_se=seg, chi2_vs_law=chi, level_se=lev_se,
                       level={v: float(np.mean(np.array(S['pred'][v][st]) - np.array(S['pred']['law'][st]))) for v in ALT})
        print('%-4s %-3s law %s  se %s | gamma se %.3f | gamma: %s | chi2(3) of alt vs law: %s | level se %.5f'
              % (name, st, ' '.join('%+.5f' % x for x in c), ' '.join('%.5f' % x for x in se), seg,
                 ' '.join('%s %.2f' % (v, gam[v]) for v in ALT), ' '.join('%s %.0f' % (v, chi[v]) for v in ALT if v != 'law'), lev_se))
    out[name] = res
# pooled over the five reserved sets and over the two height bands
for tag, sets in (('pooled 5 reserved', ('B', 'Dn', 'Fn', 'C', 'Bp')), ('band 100-600', ('B', 'Dn', 'Fn')), ('band 600-1000', ('C', 'Bp'))):
    for st in ('var', 'r'):
        w = np.array([1 / out[s][st]['gamma_se'] ** 2 for s in sets]); se = 1 / math.sqrt(w.sum())
        print('%-18s %-3s gamma se %.3f; alternatives at (gamma_alt - 1)/se: %s' % (tag, st, se, ' '.join('%s %+.1f' % (v, (np.mean([out[s][st]['gamma'][v] for s in sets]) - 1) / se) for v in ALT if v != 'law')))
        out[tag + '_' + st] = dict(gamma_se=se)
# height dependence: law ratio of the p = 2, 3, 5 amplitude between bands, same discriminants (B -> Bp, R120 -> C)
for a, b in (('B', 'Bp'), ('R120', 'C')):
    for st in ('var', 'r'):
        ca = np.array(out[a][st]['coef']['law'][:3]); cb = np.array(out[b][st]['coef']['law'][:3])
        print('height dependence %s -> %s (%s): law coefficient ratio p=2,3,5: %s' % (a, b, st, ' '.join('%.3f' % x for x in cb / ca)))
json.dump(out, open('power.json', 'w'), indent=1, default=float)
