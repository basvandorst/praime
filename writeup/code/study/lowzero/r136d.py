"""Round 136, the registered decisions D1-D3 from r136_C.json and r136_Dn.json. usage: r136d.py"""
import os, json, math
H = os.path.dirname(os.path.abspath(__file__)); C = json.load(open(os.path.join(H, 'r136_C.json'))); Dn = json.load(open(os.path.join(H, 'r136_Dn.json')))
def ivw(key):
    w = [1 / t[key + '_se'] ** 2 for t in (C, Dn)]; return sum(wi * t[key] for wi, t in zip(w, (C, Dn))) / sum(w), 1 / math.sqrt(sum(w))
r2, e2 = ivw('R2'); r1, e1 = ivw('R1')
print('D1 pooled slope, inverse-variance mean of C and Dn: %.4f +- %.4f; line 0.966: %s; distance from 0.932: %.1f errors, from 1: %.1f errors' % (r2, e2, 'amplitude 0.932 and L + 0.3 rejected' if r2 >= 0.966 else 'the BK amplitude rejected', (r2 - 0.932) / e2, (r2 - 1) / e2))
print('   the same on gamma (R1): %.4f +- %.4f' % (r1, e1))
d2 = {t['test']: (t['between'] - t['within'], t['null_sd_between_minus_within']) for t in (C, Dn)}
print('D2 between - within: ' + '; '.join('%s %+.4f (null sd %.4f, %.1f sd)' % (k, v[0], v[1], v[0] / v[1]) for k, v in d2.items()))
print('D3 conductor trend: R1(Dn) = %.3f, R1(C) = %.3f: %s' % (Dn['R1'], C['R1'], 'trend supported' if Dn['R1'] <= 0.95 and C['R1'] >= 0.97 else ('trend loses its support' if Dn['R1'] >= 0.97 else 'undecided')))
json.dump({'D1_R2': r2, 'D1_R2_se': e2, 'D1_R1': r1, 'D1_R1_se': e1, 'D2': d2, 'D3': [Dn['R1'], C['R1']]}, open(os.path.join(H, 'r136d.json'), 'w'))
