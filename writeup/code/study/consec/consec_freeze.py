"""consec_freeze.py: writes predictions_consec.json (the frozen holdout predictions) from the discovery readings and the
null laws (consec_null.py: LOS; consec_bks.py: block dispersion).  Nothing here reads any holdout data.
usage: consec_freeze.py readings_disc.json predictions_consec.json"""
import sys, os, json, math, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from consec_null import bk
from consec_bks import COMP, CHIQ, predict

D = json.load(open(sys.argv[1])); OUT = sys.argv[2]
N_DISC = D['nprimes']
RUNS = {
    'main': dict(range=[100000000, 1000000000], schemes=['C23', 'C31', 'C47', 'C39', 'G'], nprimes_expected=50847534 - 5761455,
                 command='consec_tables.py hold_main 1e8 1e9 --cases C23,C31,C47,C39,G   (repeat until "all segments done")',
                 analyse='consec_analyze.py hold_main readings_hold_main.json', bins=[[100000000, 300000000], [300000000, 1000000000]]),
    'ec': dict(range=[100000000, 400000000], schemes=['E37s', 'E37m2', 'E37m3'], nprimes_expected=21336326 - 5761455,
               command='consec_tables.py hold_ec 1e8 4e8 --cases E37   (repeat until "all segments done")',
               analyse='consec_analyze.py hold_ec readings_hold_ec.json', bins=[[100000000, 300000000], [300000000, 400000000]]),
}
P = dict(name='consecutive primes and non-abelian labels', created_from='discovery range [0, 1e8) only',
         H0=('Conditional on the abelian (congruence) part of the labels, the non-abelian parts of the labels of p_n and '
             'p_(n+k) are independent, each with its Chebotarev / Sato-Tate law, for every k, every gap and every order in '
             '1/log x (relative corrections O(x^(-1/2+eps)), heuristic); equivalently N2 = LOS on the abelian part x '
             'independent Chebotarev draws inside the fibres.'),
         H1='none derived (no non-abelian main term exists at fixed shifts; see consec.md 2.3)',
         runs=RUNS, decisions={}, readings={}, secondary={}, watch={}, power={})
# ------------------------------------------------------------------ H0 decisions (main falsifier)
P['decisions'] = {
    'D1a': 'pooled contingency chi^2 of the fibre labels inside abelian fibres (consec_analyze "indep"), all schemes, k = 1, 2, 3: '
           'pass if p >= 1e-3; H0 falsified if p < 1e-5',
    'D1b': 'every fibre-interaction cell (consec_analyze "corr", k = 1, 2, 3, all function pairs, centred r): |z| < 4.5 passes; '
           'H0 falsified if any |z| >= 5',
    'D1c': 'within part of the per-block (1e4) fibre transition tables (consec_analyze "w4"): inverse-variance pooled mean over the '
           'schemes within 4 errors of 1 (H0); each scheme within 4.5 errors; H0 falsified if the pooled mean is > 5 errors from 1',
    'D2': 'gap-conditioned lag-1 interaction: each scheme and function chi^2 p >= 1e-3; each gap = 0 mod q contrast |z| < 4',
    'D3': 'non-abelian label against the abelian datum of the neighbour (consec_analyze "nbr"): every |z| < 4.5',
    'D4': 'run lengths in the fibre subsequence against the geometric law: each scheme total p >= 1e-3',
    'D5': 'first-third interaction given the middle label (consec_analyze "triple"): each scheme p >= 1e-3',
    'D6': 'short-shift dispersion (column part at H = 1e2, 1e3): every non-abelian function R within 4 errors of 1; inverse-variance '
          'pooled R per H within 4 errors of 1 (a BKS-C2-type constant would give about 1 - 0.42/log x = 0.979)',
    'H0_verdict': 'H0 stands if D1a, D1b, D1c, D2, D3, D4, D5, D6 pass; H0 is falsified only by the falsifier clauses of D1a, D1b, D1c; '
                  'a failure of D2-D6 without a falsifier is reported as an open anomaly',
}
# expected pair counts and errors (scaled from the discovery counts)
for run, R in RUNS.items():
    sc = R['nprimes_expected'] / N_DISC
    for nm in R['schemes']:
        E = D[nm]; ent = {}
        for key, v in E['corr'].items():
            n = v['n'] * sc; ent[key] = dict(H0_r=0.0, n_expected=n, sd_expected=1 / math.sqrt(n), r_5sigma=5 / math.sqrt(n))
        P['readings'][nm + '|corr'] = ent
        P['power'][nm] = dict(scale=sc, lag1_fibre_pairs=E['corr']['k1_f0_f0']['n'] * sc,
                              detectable_r_at_5sd=5 / math.sqrt(E['corr']['k1_f0_f0']['n'] * sc),
                              abelian_b1_for_comparison=E['ab_b']['1']['b'] if 'ab_b' in E else None,
                              w4_within_err_expected=E['w4']['within_err'] / math.sqrt(sc))
# ------------------------------------------------------------------ secondary: known laws (LOS on the abelian part; BKS column)
for run, R in RUNS.items():
    x0, x1 = R['range']
    for nm in R['schemes']:
        if nm not in CHIQ or nm == 'C39': continue
        q = CHIQ[nm]; b1, b2, b3 = bk(x0, x1, 1, q), bk(x0, x1, 2, q), bk(x0, x1, 3, q)
        rho1 = (b2 - b1 * b1) / (1 - b1 * b1)
        P['secondary'][nm + '|abelian'] = dict(
            b_LOS={'1': b1, '2': b2, '3': b3},
            b_disc_over_LOS={k: D[nm]['ab_b'][k]['b'] / bk(3, 1e8, int(k), q) for k in ('1', '2', '3')},
            rule='measured b_k / LOS b_k in [0.6, 1.5] (LOS error term O(log^-7/4) is 15-25% of b_1 here; discovery ratios given)',
            w4a_within_LOS=1 + 2 * rho1, w4a_rule='abelian sign table within part within 4 errors + 0.015 of 1 + 2 rho_1, '
                                                'rho_1 = (b_2 - b_1^2)/(1 - b_1^2) (LOS triple term, consec.md 2.5)')
    bks = {}
    for nm in R['schemes']:
        if nm not in COMP: continue
        comps = COMP[nm] + ([[(1, CHIQ[nm], 1.0)]] if nm in CHIQ else [])
        for H in (100, 1000, 10000, 100000, 1000000):
            for b0, b1_ in R['bins'] + [R['range']]:
                for i, comp in enumerate(comps):
                    f = 'f%d' % i if i < len(COMP[nm]) else 'chi'
                    bks['%s|H%d_%d_%d|%s' % (nm, H, b0, b1_, f)] = dict(law=predict(comp, H, b0, b1_), comp=comp)
    P['secondary']['bks_' + run] = bks
P['secondary']['bks_rule'] = ('column / dispersion R = sum T^2 / sum V: for H = 1e5 and 1e6 in each height bin with >= 50 blocks, '
                              '|R - law| < 4 sqrt(err^2 + 0.03^2) for every function; reported, not part of the H0 verdict '
                              '(law = Montgomery-Soundararajan for d = 1, Bui-Keating-Smith Theorem C1 made local, min(1, .) kink; '
                              'H = 1e4 is near the crossover and only described)')
# ------------------------------------------------------------------ watch items from the discovery range (all below 4 sd)
P['watch'] = {
    'W1': dict(what='C31 run lengths of the class "0 roots" (fibre code 1): chi^2/dof', disc=[D['C31']['runs']['1']['chi2'], D['C31']['runs']['1']['dof']],
               real_effect='chi^2 - dof grows like the sample: (40.0 - 23) x 7.83 + 23 = 156 on 23 dof (p ~ 1e-20)',
               fluctuation='chi^2 ~ chi^2(dof), p uniform', rule='real if p < 1e-4, fluctuation if p > 0.01'),
    'W2': dict(what='C23 first-third interaction given the middle label: chi^2 on 3 dof', disc=[D['C23']['triple']['chi2'], D['C23']['triple']['dof']],
               real_effect='chi^2 ~ 3 + (8.9 - 3) x 7.83 = 49', fluctuation='chi^2(3)', rule='real if p < 1e-4, fluctuation if p > 0.01'),
    'W3': dict(what='E37m3 lag-2 interaction of f0 (a_p = 0 mod 3 indicator) with itself: z', disc=D['E37m3']['corr']['k2_f0_f0']['z'],
               real_effect='z = 2.37 x sqrt(2.70) = 3.9 (same r)', fluctuation='z ~ N(0, 1)', rule='real if z > 3, fluctuation if |z| < 2'),
    'W4': dict(what='G sector runs, sector 2 (fibre code 1): chi^2/dof', disc=[D['G']['runs']['1']['chi2'], D['G']['runs']['1']['dof']],
               real_effect='chi^2 ~ 8 + (20.3 - 8) x 7.83 = 104', fluctuation='chi^2(8)', rule='real if p < 1e-4, fluctuation if p > 0.01'),
}
json.dump(P, open(OUT, 'w'), indent=1, sort_keys=True)
print('written', OUT, 'sha256', hashlib.sha256(open(OUT, 'rb').read()).hexdigest())
