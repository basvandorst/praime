"""Freeze the predictions of the twin-class law for the fresh decade D = [2e10, 2e11) and its nine sub-windows
S_s = [2e10 (1 + s), 2e10 (2 + s)), and (dry run) for the seen decade [1e8, 1e9) with sub-windows of 1e8.
Exact Hardy-Littlewood sums for Q = 2310, 30030; closed form for Q = 210 (reported readings only).
Resumable: results accumulate in predictions_partial.json; when complete, predictions_twinclass.json is written with
the frozen error model and the decision rules, and its sha256 printed.   usage: freeze.py [fresh|dry]"""
import sys, os, json, time, hashlib, math
sys.path.insert(0, '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/twinclass')
import predict
which = sys.argv[1] if len(sys.argv) > 1 else 'fresh'
if which == 'fresh':
    base, step = 2 * 10 ** 10, 2 * 10 ** 10
    out_name = 'predictions_twinclass.json'
else:
    base, step = 10 ** 8, 10 ** 8
    out_name = 'predictions_dry.json'
WINDOWS = {'D': (base, 9 * step)}
WINDOWS.update({'S%d' % s: (base + s * step, step) for s in range(9)})
part_name = 'predictions_partial_%s.json' % which
part = json.load(open(part_name)) if os.path.exists(part_name) else {}
t0 = time.time()
for w, (L0, W) in WINDOWS.items():
    for Q in (30030, 2310, 210):
        key = '%s|%d' % (w, Q)
        if key in part:
            continue
        if time.time() - t0 > 420:
            print('time slice used; run again', flush=True); json.dump(part, open(part_name, 'w'), indent=1); sys.exit(0)
        t = time.time()
        part[key] = predict.predict_closed(L0, W, Q) if Q == 210 else predict.predict(L0, W, Q)
        json.dump(part, open(part_name, 'w'), indent=1)
        print(key, '%.1fs' % (time.time() - t), {k: round(v, 4) for k, v in part[key].items() if isinstance(v, float) and k in ('pX', 'pZ', 'pR', 'beta_HL')}, flush=True)
FROZEN = {
    'description': 'Twin counts across classes mod 210/2310/30030: linear response to the class prime counts (round 155) '
                   '+ Hardy-Littlewood remainder. Windows: D = [%d, %d) and nine equal sub-windows.' % (base, base + 9 * step),
    'windows': {w: list(v) for w, v in WINDOWS.items()},
    'predictions': part,
    'error_model': {
        'beta_30030_inflation': 1.17, 'u2Z_30030_inflation': 1.08, 'R_level3_inflation': 1.0, 'R_level2_inflation': 1.04,
        'source': '27 seen sub-windows (1e7..1e10): sd/nominal of beta 1.17, of u2Z ratio 1.08, of R level ratios 0.95/1.04/0.99',
        'cross_sd': 'nominal 2 * within_pred * sqrt((1 - sum m_s^2/m^2) / (2 (K-1))) (classes independent within and across sub-windows)'},
    'alternatives': {
        'beta_independent': 0.0, 'beta_one_member': 0.5,
        'R_independent_ratio_30030': 'u2R = u2X + u2Z under independence: (pX + pZ)/pR',
        'u2Z_poisson_primes': 'prime class counts Poisson: pZ_poisson = 2 kbar (1 - r1) -> ratio about 1.4-2.5',
        'A_twin_cross_z': {'2310': -1.26, '30030': -1.35, 'joint': -1.85,
                           'source': 'mean z of the twin remainder cross part on the decades 1e7, 1e8, 1e9 (seen)'}},
    'decisions': {
        'P1': 'beta: mean over S0..S8 of beta_hat(30030) - mean beta_HL within 3 se, se = 1.17 * rms(se_beta_nominal)/3; '
              'rivals beta = 0 and beta = 0.5 each rejected iff more than 5 se away',
        'P2': 'fine level: mean over S0..S8 of the R level-3 ratio (30030 within 2310) within 3 * sqrt(2/1350)/3 of 1',
        'P3': 'linear part: mean over S0..S8 of u2Z(30030)/pZ within 3 * 1.08 * sqrt(2/1484)/3 of 1',
        'P4': 'open point: Z = (z_2310 + z_30030)/sqrt(2), z = (cross_meas - cross_HL)/cross_sd of the remainder R on D; '
              'Z <= -3: excess anti-correlation of twin class deviations between sub-windows CONFIRMED (four-point law fails '
              'for twins in decade windows); Z >= -0.92 (midpoint to A_twin, -1.85/2): chance favoured; otherwise OPEN',
        'verdict': 'LAW CONFIRMED iff P1, P2, P3 pass; P4 is reported with its own three-way reading',
        'reported_only': 'u2X and u2R ratios at levels 1-3 for D and each S_s; mod 210 uses the closed form (+-0.01)'}}
s = json.dumps(FROZEN, indent=1, sort_keys=True)
open(out_name, 'w').write(s)
print('wrote', out_name, 'sha256', hashlib.sha256(s.encode()).hexdigest())
