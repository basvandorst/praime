import json, math, numpy as np
from scipy.stats import chi2, t as st, kstest
# --- 4. twin spreads per class mod 210
pred = {'train': 0.46, 'A': 0.71, 'B': 0.81, 'C': 0.91}; meas = {'train': 0.40, 'A': 0.18, 'B': 1.00, 'C': 0.53}
rel = math.sqrt(2/14); X2 = 0; ps = []
for w in pred:
    z = (meas[w]-pred[w])/(rel*pred[w]); X2 += z*z
    Y = 14*meas[w]/pred[w]; lo = chi2.cdf(Y, 14); p2 = 2*min(lo, 1-lo); ps.append(p2)
    print(w, 'gauss z %.2f' % z, 'chi2_14 stat %.2f  lower-tail %.4f  two-sided %.4f' % (Y, lo, p2))
print('gauss chi2 %.2f p %.3f' % (X2, chi2.sf(X2, 4)))
# likelihood-ratio / deviance for scaled chi2_14 model: D = sum 14*(r - 1 - ln r), r = meas/pred
D = sum(14*(meas[w]/pred[w] - 1 - math.log(meas[w]/pred[w])) for w in pred)
print('deviance (exact chi2_14 model) %.2f on 4 dof, p %.4f' % (D, chi2.sf(D, 4)))
F = -2*sum(math.log(p) for p in ps); print('Fisher combination of two-sided p: %.2f on 8 dof, p %.4f' % (F, chi2.sf(F, 8)))
# --- 3. round 34
for name, f in (('w1', 'r10_var.json'), ('w2', 'r10_var_3000000000.json'), ('w3', 'r10_var_5000000000.json')):
    d = json.load(open('/home/bas/praime/writeup/data/study/' + f))
    for H in ('10000',):
        m = np.array(d[H]['meas']); th = np.array(d[H]['theory']); r = m - th
        rel_ex = d[H]['mean_meas']/d[H]['mean_theory'] - 1
        se_emp = r.std(ddof=1)/math.sqrt(10)
        # expected block sd of a sample variance from n = 1e8/H intervals: rel sqrt(2/(n-1))
        n = 10**8//int(H); se_th = math.sqrt(2/(n-1))*th.mean()/math.sqrt(10)
        print(name, H, 'rel %.2f%%  t9 %.2f  se_emp %.5f  se(sqrt(2/n)) %.5f  z_th %.2f  block sd/expected %.2f' % (100*rel_ex, d[H]['t9'], se_emp, se_th, r.mean()/se_th, r.std(ddof=1)/(se_th*math.sqrt(10))))
# --- 2. extrapolation of U
d = json.load(open('/home/bas/praime/writeup/data/study/r32b.json'))
U3 = np.array(d['U3']); U4 = np.array(d['U4']); r = U4/U3
for f, lab in ((1.0, 'higher orders = order 4 again'), (0.72, 'r=0.42'), (1.27, 'r=0.56'), (1.5, 'r=0.6')):
    U = U4*r**f; print(lab, 'mean %.3f se %.3f  z %.2f  KS p %.3f  median %.2f' % (U.mean(), U.std(ddof=1)/math.sqrt(61), (U.mean()-1)/(U.std(ddof=1)/math.sqrt(61)), kstest(U, 'expon').pvalue, np.median(U)))
print('order-3 on 61 gaps: z of mean %.2f ; KS p %.3f' % ((U3.mean()-1)/(U3.std(ddof=1)/math.sqrt(61)), kstest(U3, 'expon').pvalue))
print('share of excess removed by order 4: %.2f' % ((U3.mean()-U4.mean())/(U3.mean()-1)))
# --- 5. scoreboard sums
rows = [(96,4,1.4),(68,0,1.0),(45,3,0.7),(20,0,0.3),(4,0,0.0),(8,0,0.1),(9,1,0.1),(36,0,0.5),(8,1,0.0),(56,1,0.8),(4,0,0.0),(36,1,0.5),(18,1,0.3),(37,0,0.6),(56,17,0.8),(55,1,0.8),(145,2,2.2),(2,1,0.0),(2,1,0.0),(9,0,0.1),(2,0,0.0),(11,0,0.0)]
print('sum tests %d misses %d expected %.1f' % tuple(map(sum, zip(*rows))))
