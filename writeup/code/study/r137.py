"""Round 137 (preregister/round137.txt): the within part of the transition tables on 24 windows [6.0e9 + i 1e8, + 1e8),
from the outputs of `r135b.py new1..new4` (six windows each), against the frozen predictions in within/data/.
Readings R1-R5 and the falsifier as registered. usage: r137.py"""
import os, json, math
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(H, '..', '..', 'data', 'study')
J = [json.load(open(os.path.join(D, 'r135b_new%d.json' % i))) for i in (1, 2, 3, 4)]
P = json.load(open(os.path.join(H, 'within', 'data', 'predictions.json'))); RT = json.load(open(os.path.join(H, 'within', 'data', 'reg_T.json')))
assert sum(len(j['windows']) for j in J) == 24 and [w for j in J for w in j['windows']] == [6 * 10 ** 9 + i * 10 ** 8 for i in range(24)]
CELLS = [(q, k, BL) for q, k in ((3, 1), (5, 1), (7, 1), (3, 2), (10, 2)) for BL in (10 ** 4, 10 ** 5, 10 ** 6)]
out = {}
for q, k, BL in CELLS:
    key = 'win_%d_%d_%d' % (q, k, BL); p = P['%d_%d_%d' % (5 if q == 10 else q, k, BL)]
    pw = np.array([v for j in J for v in j[key]['per_window']]); pe = np.array([v for j in J for v in j[key]['per_window_err']])
    m = float(pw.mean()); sb = float(math.sqrt((pe ** 2).sum()) / 24); ss = float(pw.std(ddof=1) / math.sqrt(24)); sm = max(sb, ss)
    out[(q, k, BL)] = {'meas': m, 'sig_block': sb, 'sig_scatter': ss, 'sig_meas': sm, 'sig_tot': math.hypot(sm, p['law_err']), 'law': p['law'], 'law_err': p['law_err'],
                       'MQ': p['MQ'], 'old': p['old_meas'], 'lead': p['lead']}
    c = out[(q, k, BL)]; c['z'] = (m - c['law']) / c['sig_tot']
print('cell           measured   sig_meas (block, scatter)   law +- err      z     | 1       sieve chain  leading  old height')
for (q, k, BL), c in out.items():
    print('(%2d,%d,%.0e)  %.4f     %.4f  (%.4f, %.4f)     %.4f +- %.4f  %+5.2f | %+6.1f  %.4f (%+.1f)  %.4f   %.4f' % (q, k, BL, c['meas'], c['sig_meas'], c['sig_block'], c['sig_scatter'], c['law'], c['law_err'], c['z'],
          (c['meas'] - 1) / c['sig_meas'], c['MQ'], (c['meas'] - c['MQ']) / c['sig_tot'], c['lead'], c['old']))
R = {}; c1 = [out[(q, 1, 10 ** 4)] for q in (3, 5, 7)]
R['R1_pass'] = [bool(abs(c['z']) <= 2.5) for c in c1]; R['R1_rejects_1'] = [bool(c['meas'] <= 1 - 3 * c['sig_meas']) for c in c1]
R['R1_rejects_sieve'] = [bool(c['meas'] <= c['MQ'] - 2.5 * c['sig_tot']) for c in c1]
T = sum(w * (c['meas'] - c['old']) for w, c in zip(RT['weights'], c1)) / RT['den']; sT = math.sqrt(sum((w * c['sig_meas']) ** 2 for w, c in zip(RT['weights'], c1))) / RT['den']
R['R2_T'] = T; R['R2_T_meas_err'] = sT; R['R2_pass'] = bool(0.40 <= T <= 1.30); R['falsified'] = bool(T <= 0)
Tlead = sum(w * (c['meas'] - c['lead']) for w, c in zip(RT['weights'], c1)) / RT['den']
S = float(np.mean([out[(q, 1, 10 ** 5)]['meas'] - out[(q, 1, 10 ** 4)]['meas'] for q in (5, 7)]))
sS = math.sqrt(sum(out[(q, 1, B)]['sig_meas'] ** 2 for q in (5, 7) for B in (10 ** 4, 10 ** 5))) / 2
R['R3_S'] = S; R['R3_S_meas_err'] = sS; R['R3_pass'] = bool(-0.0227 <= S <= -0.0009); R['R3_rejects_no_tail'] = bool(S <= -0.0002 - 3 * sS)
a, b = out[(10, 2, 10 ** 4)], out[(3, 2, 10 ** 4)]
R['R4_pass'] = [bool(abs(a['z']) <= 2.5), bool(abs(b['z']) <= 2.5)]; R['R4_10_2_rejects_1'] = bool(a['meas'] >= 1 + 3 * a['sig_meas']); R['R4_3_2_rejects_sieve'] = bool(b['meas'] <= b['MQ'] - 2.5 * b['sig_tot'])
R['R5_pass'] = [bool(abs(out[c]['z']) <= 2.5) for c in ((3, 1, 10 ** 5), (5, 1, 10 ** 5), (7, 1, 10 ** 5), (10, 2, 10 ** 5))]
print('R1 level (3,5,7 at k = 1, H = 1e4): within 2.5 sig_tot %s; rejects 1 %s; rejects the sieve chain %s' % (R['R1_pass'], R['R1_rejects_1'], R['R1_rejects_sieve']))
print('R2 height contrast T = %.3f +- %.3f (measurement); law 0.851 +- 0.181, no height dependence 0, leading order -0.42, sieve chain -0.63; pass %s; falsified %s' % (T, sT, R['R2_pass'], R['falsified']))
print('R3 slope S = %.4f +- %.4f (measurement); law -0.0109; pass %s; rejects no tail %s' % (S, sS, R['R3_pass'], R['R3_rejects_no_tail']))
print('R4 k = 2: (10,2) and (3,2) within 2.5 sig_tot %s; (10,2) rejects 1 %s; (3,2) rejects the sieve chain %s' % (R['R4_pass'], R['R4_10_2_rejects_1'], R['R4_3_2_rejects_sieve']))
print('R5 the four cells at H = 1e5 within 2.5 sig_tot: %s' % R['R5_pass'])
print('verdict: law confirmed (R1-R4 pass) %s' % bool(all(R['R1_pass']) and R['R2_pass'] and R['R3_pass'] and all(R['R4_pass'])))
json.dump({'cells': {'%d_%d_%d' % k: v for k, v in out.items()}, 'readings': R}, open(os.path.join(D, 'r137.json'), 'w'), indent=0)
