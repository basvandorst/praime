"""expected errors of every registered decision from the frozen noise model (no zero used). usage: b5_power.py predictions_blind.json"""
import sys, json, math
import numpy as np
P = json.load(open(sys.argv[1])); L = P['L']
def pooled(keys):
    num = sum(abs(complex(*L[k]['pred'])) ** 2 for k in keys)
    return math.sqrt(sum(abs(complex(*L[k]['pred'])) ** 2 * L[k]['sd']['along'] ** 2 for k in keys)) / num, \
           math.sqrt(sum(abs(complex(*L[k]['pred'])) ** 2 * L[k]['sd']['across'] ** 2 for k in keys)) / num
neg = [k for k, e in L.items() if e['eps'][0] < -0.5 and abs(e['eps'][1]) < 1e-9]
cpx = [k for k, e in L.items() if abs(e['eps'][1]) > 1e-9]
print('D6 eps = -1 subset: n = %d, se %.4f (H-noroot predicts -1: separation %.0f se)' % (len(neg), pooled(neg)[0], 2 / pooled(neg)[0]))
se, sei = pooled(cpx); num = sum(abs(complex(*L[k]['pred'])) ** 2 for k in cpx)
hc = sum(complex(*L[k]['alt']['H_conj']) * complex(*L[k]['pred']).conjugate() for k in cpx) / num
hn = sum(complex(*L[k]['alt']['H_noroot']) * complex(*L[k]['pred']).conjugate() for k in cpx) / num
print('D7 complex eps: n = %d, se(Re s) %.4f, se(Im s) %.4f; H-conj gives s = %.3f%+.3fi, H-noroot s = %.3f%+.3fi' % (len(cpx), se, sei, hc.real, hc.imag, hn.real, hn.imag))
X = []; W = []
for k, e in L.items():
    unit = e['T2'] / (2 * math.pi); X.append([1, e['c'] - 1]); W.append((unit / e['sd']['along']) ** 2)
X = np.array(X); W = np.array(W); cov = np.linalg.inv(X.T @ (W[:, None] * X))
print('D5 amplitude follows c: se(alpha) %.4f, se(beta) %.4f (H-c1: beta = 0 at %.0f se)' % (math.sqrt(cov[0, 0]), math.sqrt(cov[1, 1]), 1 / math.sqrt(cov[1, 1])))
for f in ('T1', 'T2', 'T3', 'T4a', 'T4b'):
    ks = [k for k in L if L[k]['cls'] == f]; s = pooled(ks)[0]
    print('D1 %s: n = %d, se %.4f; PASS band 1 +- %.3f; falsifier band 1 +- %.3f; H0 at %.0f se' % (f, len(ks), s, 3 * s, 4 * s, 1 / s))
s = pooled(list(L))[0]
print('D2 all primitive: n = %d, se %.4f; PASS band 1 +- %.3f; falsifier 1 +- %.3f; D3 H0 at %.0f se' % (len(L), s, 3 * s, 4 * s, 1 / s))
print('D4 lines: se %.4f' % P['expected_se']['lines'])
for rd in ('M1', 'M2'):
    print('D8/D9 %s: se %.3f (zero hypothesis at %.1f se if the mechanism holds)' % (rd, P['expected_se']['T5_' + rd], 1 / P['expected_se']['T5_' + rd]))
T5 = [e for e in P['products'].values() if e['cls'] == 'T5']
print('D9 M3: se of the mean %.2f' % (np.mean([e['readings']['M3']['sd'] for e in T5]) / math.sqrt(len(T5))))
print('D10 Artin law on the 29 fields with c != 0: se %.4f' % P['expected_se']['T5_R_artin'])
z = [e for e in T5 if abs(e['rho']['c']) < 1e-9]
for e in z: print('D10b c_rho = 0 field d = %d: sd %.2f; a c = 1 law would give %.1f (%.1f sd)' % (e['d'], e['readings']['R_artin']['sd'], -e['T2'] / (2 * math.pi), e['T2'] / (2 * math.pi) / e['readings']['R_artin']['sd']))
g = [e for e in P['products'].values() if e['cls'] == 'T5g']
num = sum(e['readings']['half_phase']['pred'] ** 2 for e in g); var = sum(e['readings']['half_phase']['pred'] ** 2 * e['readings']['half_phase']['sd'] ** 2 for e in g)
print('D11 T5g half phase: pooled se %.4f; product phase sd %s' % (math.sqrt(var) / num, [round(e['readings']['product_phase']['sd'], 1) for e in g]))
print('D12 pairs: se of the mean about sqrt(2 x 800)/sqrt(10) = %.1f' % (math.sqrt(1600) / math.sqrt(10)))
