# Build predictions_cubic.json (frozen predictions, comparators, expected holdout errors, decisions) from
# pred_cache.json (theory) and the discovery counts (errors and the one fitted comparator amplitude only).
import os, sys, json, math, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from finiteX import epsp, jack, coarse, cTT, KC
pc = json.load(open(os.path.join(HERE, 'pred_cache.json'))); edges = np.array(pc['edges']); yc = (edges[:-1] + edges[1:]) / 2
HI = (8, 60); LO = (0, 4); G = 4
kap = KC[1] / (12 * 1.2020569031595942)
def count(sgn, X):  # Roberts/TT: number of fields with X <= |D| < 2X
    C = 1 / (4 * 1.2020569031595942) if sgn < 0 else 1 / (12 * 1.2020569031595942); K = math.sqrt(3) * kap if sgn < 0 else kap
    return C * X + K * ((2 * X)**(5 / 6) - X**(5 / 6))
out = {'what': 'S3 cubic fields: E[sqrt|D| a_K(p)] in bins of y = p/|D|, p <= 3|D| unramified, fields X <= |D| < 2X by sign',
       'edges': edges.tolist(), 'coarse_group': G, 'range_hi_bins': list(HI), 'range_lo_bins': list(LO),
       'pred': {}, 'numerical_error': {}, 'finiteX_factor': {}, 'comparators': {}, 'sigma_hold': {}, 'power': {}, 'decisions': {}}
for sgn in (-1, 1):
    k = '%+d' % sgn; z = np.load(os.path.join(HERE, 'cubic_disc_%s.npz' % k)); num, den, D = z['num'], z['den'], np.abs(z['D'])
    P = np.array(pc[k]['pred']); out['pred'][k] = P.tolist(); out['numerical_error'][k] = pc[k]['numerr']
    out['finiteX_factor'][k] = {'disc': 1 / (1 + epsp(sgn, 5000.0)), 'hold': 1 / (1 + epsp(sgn, 20000.0))}
    nd = len(D); nh = count(sgn, 20000.0)
    m, cov, sig, s, se = jack(num, den, P)
    scale = np.sqrt((nd / nh) * np.log(yc * 30000.0) / np.log(yc * D.mean()))
    out['sigma_hold'][k] = (sig * scale).tolist()
    sl = slice(*HI); U = np.array(pc[k]['single_line_unit']); w = 1 / sig[sl]**2
    a1 = float(np.sum(U[sl] * m[sl] * w) / np.sum(U[sl]**2 * w))
    e_h = epsp(sgn, 20000.0)
    comps = {'zero': 0 * P, 'single_line_zeta3': np.array(pc[k]['single_line_zeta3']), 'single_line_fitted_on_discovery': a1 * U,
             'gl2_prime_level_transplant': np.array(pc[k]['gl2_prime_level']), 'tt_secondary_only': np.array(pc[k]['tt_secondary']) / (1 + e_h),
             'lines_without_prime_square_term': np.array(pc[k]['lines_no_square_term'])}
    out['comparators'][k] = {n: v.tolist() for n, v in comps.items()}
    out.setdefault('comparator_notes', {'single_line_fitted_on_discovery': 'a * 2 pi K(4 pi sqrt y), a fitted on discovery y in (0.4, 3]',
        'tt_secondary_only': 'c y^(-1/3)/(1+eps_hold), c = -1.02495 (D<0), -1.77526 (D>0)', 'single_line_zeta3': 'r = 1 line alone, amplitude zeta(3)',
        'gl2_prime_level_transplant': 'GL(2) prime-level natural-weight lines (s/d)^2 with the kernel of the sign', 'lines_without_prime_square_term': 'prediction + 1/(2 sqrt y)'})
    out['comparator_notes']['fitted_amplitude_%s' % k] = a1
    # expected holdout errors of the decision statistics, from the discovery jackknife scaled by the field counts
    mc, covc, _, _, _ = jack(coarse(num[:, sl]), coarse(den[:, sl]), coarse(P[sl]) / G)
    sc2 = (nd / nh) * np.log(0.5 * 30000) / np.log(0.5 * D.mean())
    covh = covc * sc2
    _, _, _, s_hi, se_hi = jack(num[:, sl], den[:, sl], P[sl]); _, _, _, s_lo, se_lo = jack(num[:, slice(*LO)], den[:, slice(*LO)], P[slice(*LO)])
    pw = {'n_fields_hold_expected': nh, 'n_fields_disc': nd, 's_hi_err_hold': se_hi * math.sqrt(sc2), 's_lo_err_hold': se_lo * math.sqrt(sc2),
          'disc_s_hi': s_hi, 'disc_s_hi_err': se_hi, 'disc_s_lo': s_lo, 'disc_s_lo_err': se_lo}
    for n, C in comps.items():
        d = coarse(P[sl] - C[sl]) / G; pw['expected_delta_chi2_' + n] = float(d @ np.linalg.solve(covh, d))
    pw['separation_finiteX_vs_limit_sigma'] = (out['finiteX_factor'][k]['hold'] - 1) / pw['s_lo_err_hold']
    out['power'][k] = pw
    out['decisions'][k] = {'primary': sgn < 0, 'amp_z': 3.0, 'amp_signal': 10.0, 'shape_chi2_max': 32.91, 'delta_min': 25.0,
                           'lo_z': 3.0, 'lo_excess_z': 3.0, 'lo_primary': True}
out['decision_text'] = ('Per sign. HI = y in (0.4, 3] (52 bins; 13 coarse bins of width 0.2), LO = y in (0, 0.2]. '
  'A (amplitude): s_HI = <P,m>/<P,P> (weights 1/sigma^2, jackknife over fields): |s_HI - 1| <= 3 sigma and s_HI/sigma >= 10. '
  'B (shape): coarse chi2 of m - s_HI P with the jackknife covariance <= 32.91 (99.9% of chi2, 12 dof). '
  'C (comparators): coarse chi2(comparator) - chi2(P) >= 25 for every comparator. '
  'D (finite X at small y): s_LO / f_hold within 3 sigma of 1 (f = 1/(1+eps\'), eps\' from the Roberts secondary term) and (s_LO - 1) >= 3 sigma. '
  'PRIMARY CONFIRMED iff A, B, C for D < 0 and D for both signs pass; A-C for D > 0 are secondary (discovery: s_HI = 1.19 +- 0.07).')
out['falsifier'] = ('If for D < 0 on y in (0.4, 3] at 20000 <= |D| < 40000 the amplitude of the limit formula misses 1 by more than 3 sigma, '
  'or the shape chi2 exceeds 32.91, the predicted murmuration density (averaged one-point law with lines, Davenport-Heilbronn local masses) is falsified '
  'at this conductor; if s_LO misses 1/(1+eps\') by more than 3 sigma for either sign, the finite-X reading of the small-y excess is falsified.')
json.dump(out, open(os.path.join(HERE, 'predictions_cubic.json'), 'w'), indent=1)
for k in ('-1', '+1'): print(k, json.dumps({a: (round(b, 3) if isinstance(b, float) else b) for a, b in out['power'][k].items()}), ' f =', out['finiteX_factor'][k])
