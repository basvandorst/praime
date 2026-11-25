"""Power and expected readings on the fresh set: deg2_test.py on synthetic data = hypothesis_i + noise drawn from the
development noise model (noise_model.json, scaled by expected zeros), for each hypothesis and SEEDS noise draws.
Reports per reading: mean over draws, sd over draws (the true error), mean bootstrap se (what the test will quote).
usage: power.py PRED.json NOISE.json OUTDIR SEEDS"""
import sys, os, json, subprocess
import numpy as np
P, NM, OUT, K = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]); os.makedirs(OUT, exist_ok=True)
H = os.path.dirname(os.path.abspath(__file__)); py = sys.executable
HYP = ['full', 'k1', 'deg1', 'k1p2', 'w1p2', 'offzeta', 'gue']
res = {}
for h in HYP:
    for s in range(K):
        f = os.path.join(OUT, '%s_%d.json' % (h, s))
        if not os.path.exists(f):
            subprocess.run([py, '-I', os.path.join(H, 'deg2_test.py'), P, 'none', f, '--synthetic', h, '--noise', NM,
                            '--seed', str(1000 + s), '--nboot', '300'], check=True, capture_output=True)
        res.setdefault(h, []).append(json.load(open(f)))
keys = ['gamma_A', 'gamma_deg1', 'gamma_k1', 'gamma_k1p2', 'gamma_k1p3', 'gamma_w1p2', 'gamma_w1p3', 'gamma_offzeta']
lines = []
def say(s): print(s); lines.append(s)
r0 = res['full'][0]; L = r0['l']
say('fresh set: %d curves; per hypothesis %d noise draws; noise model %s' % (r0['n'], K, NM))
say('%-14s ' % 'reading' + ' '.join('%-24s' % h for h in HYP))
for k in keys:
    row = []
    for h in HYP:
        v = np.array([x['gammas'][k][0] for x in res[h]]); b = np.mean([x['gammas'][k][1] for x in res[h]])
        row.append('%7.3f (sd %.3f, se %.3f)' % (v.mean(), v.std(ddof=1) if K > 1 else 0, b))
    say('%-14s ' % k + ' '.join('%-24s' % c for c in row))
for k in ('chi2_p2_full', 'chi2_p2_gue', 'chi2_p2_deg1', 'chi2_p2_k1', 'chi2_p3_full', 'chi2_p3_gue', 'chi2_p3_deg1', 'chi2_p3_k1'):
    say('%-14s ' % k + ' '.join('%-24s' % ('%9.1f' % np.mean([x['chi2'][k] for x in res[h]])) for h in HYP))
say('T1 level: expected se per length (full): ' + ' '.join('%.4f' % v for v in np.mean([x['level_se'] for x in res['full']], 0)))
say('T1 Hotelling p (full, per draw): ' + ' '.join('%.3f' % x['p_T2'] for x in res['full']))
for h in HYP[1:]:
    lv = np.mean([x['level'] for x in res[h]], 0)
    say('   level under %-8s (data - law): ' % h + ' '.join('%+.4f' % v for v in lv))
lr = res['full'][0]['lessrigid']
say('less rigid than GUE: %d curves; law mean (Sigma2 - GUE) at l = 2..10 step 2: ' % lr['n'] + ' '.join('%+.4f' % lr['law'][L.index(l)] for l in (2.0, 4.0, 6.0, 8.0, 10.0))
    + '; se ' + ' '.join('%.4f' % np.mean([x['lessrigid']['se'][L.index(l)] for x in res['full']]) for l in (2.0, 4.0, 6.0, 8.0, 10.0)))
say('decisions passed (fraction of draws):')
for h in HYP:
    d = {}
    for x in res[h]:
        for kk, v in x['decisions'].items():
            if isinstance(v, bool): d[kk] = d.get(kk, 0) + v / K
    say('  %-8s ' % h + ' '.join('%s=%.2f' % (kk.split('_', 1)[0] + '_' + kk.split('_', 1)[1][:14], v) for kk, v in d.items()))
open(os.path.join(OUT, 'power_out.txt'), 'w').write('\n'.join(lines) + '\n')
