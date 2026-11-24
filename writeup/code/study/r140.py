"""Round 140 (preregister/round140.txt): the registered readings, applied to the map.py report of the 24 windows at
2.40e11 (../../data/study/r140/rep_fresh.json) and the frozen model values (allorder/reg_numbers.json).
usage: r140.py"""
import os, json, math
H = os.path.dirname(os.path.abspath(__file__)); R = json.load(open(os.path.join(H, '..', '..', 'data', 'study', 'r140', 'rep_fresh.json')))
F = json.load(open(os.path.join(H, 'allorder', 'reg_numbers.json')))['cells']; C = R['cells']
assert sorted(R['windows']) == [240000000000 + i * 100000000 for i in range(24)]
out = {}; chi = 0.0; P14 = []
print('cell   measured slope      model CG(1e5)       z      comparator y=17   first order   | W(1e4) measured   model')
for cell in ('7_1', '9_1', '11_1', '13_1', '11_2', '13_2'):
    a = C[cell]['all']; s, e = a['slope']; m, me = F[cell]['cg_slope']; z = (s - m) / e; lr = a['lr']; fo = lr + 10.625 * (lr - a['lr13'])
    w, we = a['W'][0]; wm, wme = F[cell]['cg_W1e4']
    out[cell] = {'slope': s, 'err': e, 'model': m, 'model_mc': me, 'z': z, 'lr17': lr, 'first_order': fo, 'W1e4': w, 'W_err': we, 'W_model': wm, 'pass': bool(abs(z) <= 2.5), 'W_z': (w - wm) / we}
    chi += z * z
    if cell.endswith('_1'): P14.append((s, e, lr, fo, m))
    print('%-5s  %+.5f +- %.5f  %+.5f +- %.5f  %+.2f   %+.5f          %+.5f      | %.4f +- %.4f   %.4f (z %+.1f)' % (cell, s, e, m, me, z, lr, fo, w, we, wm, (w - wm) / we))
def pooled(f):
    wt = [1 / e ** 2 for _, e, _, _, _ in P14]; v = sum(w * f(t) for w, t in zip(wt, P14)) / sum(wt); return v, 1 / math.sqrt(sum(wt))
x1, e1 = pooled(lambda t: t[0] - t[2]); x2, e2 = pooled(lambda t: t[0] - t[3]); x3, e3 = pooled(lambda t: t[0] - t[4])
out['P_pass'] = [out[c]['pass'] for c in ('7_1', '9_1', '11_1', '13_1', '11_2', '13_2')]; out['J_chi2'] = chi; out['J_pass'] = bool(chi <= 16.8)
out['F1_pooled_excess_over_comparator'] = [x1, e1]; out['F1_fires'] = bool(abs(x1) <= e1)
out['F2_pooled_minus_first_order'] = [x2, e2]; out['pooled_minus_model'] = [x3, e3]; out['F2_fires'] = bool(abs(x2) <= e2 and x3 < -3 * e3)
out['levels_pass'] = [bool(abs(out[c]['W_z']) <= 2.5) for c in ('7_1', '11_1', '13_1')]
out['law_falsified'] = bool((not out['J_pass'] and sum(not p for p in out['P_pass'][:4]) >= 2) or out['F2_fires'])
print('P1-P6 within 2.5 errors: %s;  J: chi2 %.1f on 6 (limit 16.8): %s' % (out['P_pass'], chi, 'pass' if out['J_pass'] else 'FAIL'))
print('F1 pooled slope - comparator(y = 17), k = 1: %+.5f +- %.5f (%.1f errors): the comparator suffices: %s' % (x1, e1, x1 / e1, out['F1_fires']))
print('F2 pooled slope - first order: %+.5f +- %.5f (%.1f errors); pooled slope - model: %+.5f +- %.5f (%.1f errors): first order holds: %s' % (x2, e2, x2 / e2, x3, e3, x3 / e3, out['F2_fires']))
print('levels W(1e4) for 7_1, 11_1, 13_1 within 2.5 errors: %s;  law falsified: %s' % (out['levels_pass'], out['law_falsified']))
for cell in ('5_1', '7_2', '9_2', '5_2'):
    a = C[cell]['all']; print('no reading  %-4s slope %+.5f +- %.5f   comparator %+.5f' % (cell, a['slope'][0], a['slope'][1], a['lr']))
json.dump(out, open(os.path.join(H, '..', '..', 'data', 'study', 'r140', 'readings.json'), 'w'), indent=0)
