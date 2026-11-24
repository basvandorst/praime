"""Exploratory summaries for the registration (seen data only): reads the report JSONs written by map.py.
- combined map over the two heights (inverse-variance) per modulus, lag and direction type;
- error calibration from the permutation and block-shuffle controls (z of every slope should be N(0,1));
- the chain null per type (bias of the slope from O(1/N) terms);
- the mod-5 mixed (rc) excess at the two heights, its scaling fits and the expected value and power at the fresh height."""
import json, math, sys
import numpy as np

def load(f): return json.load(open(f))

NEW, OLD = load('rep_new.json'), load('rep_old.json')
lam_new, lam_old = NEW['lam'], OLD['lam']
lam_f = float(np.mean([math.log(3.0e10 + i * 1e8 + 5e7) for i in range(24)]))
lam_resp_f = math.log(2.9e10 + 5e7)
print('heights: lam new %.4f (24 windows), old %.4f (9 windows), fresh %.4f (24 windows), fresh response stretch %.4f'
      % (lam_new, lam_old, lam_f, lam_resp_f))

def Phi(x): return 0.5 * (1 + math.erf(x / math.sqrt(2)))

# ---------------------------------------------------------------- 1. combined map
print('\n1. EXPLORATORY MAP (seen data): excess = measured slope - LR slope per unit ln H, per direction, by type')
print('   cell     type  ndir | new (lam 22.69)            | old (lam 27.02)            | combined                 | LR new   y13/y17')
sel = []
rows = {}
for cell in NEW['cells']:
    q, k = map(int, cell.split('_'))
    for t, a in NEW['cells'][cell]['types'].items():
        b = OLD['cells'][cell]['types'][t]
        w1, w2 = 1 / a['err_tot'] ** 2, 1 / b['err_tot'] ** 2
        m = (w1 * a['excess'] + w2 * b['excess']) / (w1 + w2); e = 1 / math.sqrt(w1 + w2)
        rows[(q, k, t)] = dict(new=(a['excess'], a['err_tot']), old=(b['excess'], b['err_tot']), comb=(m, e), ndir=a['ndir'],
                               lr=a['lr'], lr13=a['lr13'], lr_cross=a['lr_cross'], slope=a['slope'])
        flag = ' <==' if abs(m / e) > 2.5 else ''
        if abs(m / e) > 2.5 and k in (1, 2): sel.append((q, k, t))
        print('   %2d k=%-2d %-5s %3d | %+.4f +- %.4f (%+5.1f) | %+.4f +- %.4f (%+5.1f) | %+.4f +- %.4f (%+5.1f) | %+.4f  %.2f%s'
              % (q, k, t, a['ndir'], a['excess'], a['err_tot'], a['excess'] / a['err_tot'], b['excess'], b['err_tot'], b['excess'] / b['err_tot'],
                 m, e, m / e, a['lr'], a['lr13'] / a['lr'], flag))
print('   selected (|combined z| > 2.5, k = 1, 2):', sel)

# ---------------------------------------------------------------- 2. error calibration
print('\n2. ERROR CALIBRATION: z = slope/err for every (q, k, type) and every group in the controls (expected N(0,1))')
for name in ('perm', 'bshuf'):
    R = load('rep_%s.json' % name); zt = []; zg = []
    for cell, c in R['cells'].items():
        for t, a in c['types'].items(): zt.append(a['slope'][0] / a['slope'][1])
        for g in c['groups']: zg.append(g['slope'][0] / g['slope'][1])
    zt, zg = np.array(zt), np.array(zg)
    print('   %-6s types: n %3d  mean z %+.2f  rms z %.2f  |z|>2: %d (exp %.1f)   groups: n %3d  mean %+.2f  rms %.2f  |z|>2: %d (exp %.1f)'
          % (name, len(zt), zt.mean(), math.sqrt((zt ** 2).mean()), (abs(zt) > 2).sum(), 0.0455 * len(zt),
             len(zg), zg.mean(), math.sqrt((zg ** 2).mean()), (abs(zg) > 2).sum(), 0.0455 * len(zg)))
    if name == 'perm':
        zw = []
        for cell, c in R['cells'].items():
            for t, a in c['types'].items(): zw.append((a['W'][0][0] - 1) / a['W'][0][1])
        zw = np.array(zw)
        print('          level at 1e4 (W - 1)/err over types: mean %+.2f rms %.2f' % (zw.mean(), math.sqrt((zw ** 2).mean())))

# ---------------------------------------------------------------- 3. chain null
try:
    CH = load('rep_chain.json')
    print('\n3. CHAIN NULL (%d simulated windows, first-order chain with the primes\' pooled matrix; no long range)' % len(CH['windows']))
    print('   cell     type  ndir   W(1e4)-1           slope (= the excess of a chain)   z')
    zs = []
    for cell, c in CH['cells'].items():
        q, k = map(int, cell.split('_'))
        for t, a in c['types'].items():
            z = a['slope'][0] / a['slope'][1]; zs.append(z)
            print('   %2d k=%-2d %-5s %3d  %+.5f +- %.5f   %+.5f +- %.5f   %+5.1f' % (q, k, t, a['ndir'], a['W'][0][0] - 1, a['W'][0][1], a['slope'][0], a['slope'][1], z))
    zs = np.array(zs); print('   all type cells: rms z %.2f, |z| > 2: %d of %d' % (math.sqrt((zs ** 2).mean()), (abs(zs) > 2).sum(), len(zs)))
except FileNotFoundError:
    print('\n3. chain null: rep_chain.json missing')

# ---------------------------------------------------------------- 4. the primary reading
print('\n4. PRIMARY: mod 5, k = 1, type rc (the four mixed directions), excess per unit ln H per direction')
r = rows[(5, 1, 'rc')]
(e1, s1), (e2, s2) = r['new'], r['old']
print('   new %+.4f +- %.4f   old %+.4f +- %.4f   combined %+.4f +- %.4f' % (e1, s1, e2, s2, *r['comb']))
for name, f in (('constant', lambda l: 1.0), ('1/lam', lambda l: 1 / l), ('1/lam^2', lambda l: 1 / l ** 2)):
    # least squares amplitude A with excess = A f(lam)
    A = (e1 * f(lam_new) / s1 ** 2 + e2 * f(lam_old) / s2 ** 2) / (f(lam_new) ** 2 / s1 ** 2 + f(lam_old) ** 2 / s2 ** 2)
    sA = 1 / math.sqrt(f(lam_new) ** 2 / s1 ** 2 + f(lam_old) ** 2 / s2 ** 2)
    chi2 = ((e1 - A * f(lam_new)) / s1) ** 2 + ((e2 - A * f(lam_old)) / s2) ** 2
    print('   scaling %-9s: expected at lam %.2f: %+.4f +- %.4f  (chi2 %.2f on 1 dof)' % (name, lam_f, A * f(lam_f), sA * f(lam_f), chi2))
# also the "new height alone" scaled
for name, f in (('1/lam', lambda l: 1 / l), ('1/lam^2', lambda l: 1 / l ** 2)):
    print('   from the new height alone, %s: %+.4f' % (name, e1 * f(lam_f) / f(lam_new)))
sig_f = s1                                                   # 24 windows, N-independent per-window noise
print('   expected error at the fresh height: %.4f (the 24-window error of the new height; the per-block noise does not depend on N)' % sig_f)
for E in (-0.0060, -0.0080, -0.0085, -0.0090):
    pw = Phi((-2 * sig_f - E) / sig_f); pf = Phi((sig_f - E) / sig_f) - Phi((-sig_f - E) / sig_f)
    print('   if the true excess is %+.4f: P(pass: excess <= -2 err) = %.3f, P(falsifier: |excess| < 1 err) = %.4f' % (E, pw, pf))
E0 = 0.0
print('   if the true excess is 0: P(pass) = %.3f, P(falsifier) = %.3f' % (Phi(-2), Phi(1) - Phi(-1)))

# secondary readings: expected z at the fresh height
print('\n5. SECONDARY candidates (|combined z| > 2.5), with the new-height error as the expected fresh error:')
for key in sel:
    r = rows[key]; m, e = r['comb']; sf = r['new'][1]
    pk = Phi((-2 * sf - m) / sf) if m < 0 else 1 - Phi((2 * sf - m) / sf)
    print('   q=%d k=%d %-5s combined %+.4f +- %.4f (%+.1f); new %+.4f +- %.4f; old %+.4f +- %.4f; P(same sign at >= 2 err | combined value) = %.2f'
          % (key[0], key[1], key[2], m, e, m / e, *r['new'], *r['old'], pk))
json.dump({'%d_%d_%s' % k: v for k, v in rows.items()}, open('map_combined.json', 'w'), indent=0)
