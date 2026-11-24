import json, math, numpy as np
from scipy.stats import norm
A = json.load(open('agg.json'))
P = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/within/data/predictions.json'))
RT = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/within/data/reg_T.json'))
OLDJ = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/data/study/r135c.json'))['all']
OLD = [2e11, 3e11, 4e11, 5.5e11, 6e11, 7e11, 8e11, 9e11, 1e12]; NEW = [6e9 + i*1e8 for i in range(24)]
lo = np.mean([math.log(x+5e7) for x in OLD]); ln = np.mean([math.log(x+5e7) for x in NEW])
ilo = np.mean([1/math.log(x+5e7) for x in OLD]); iln = np.mean([1/math.log(x+5e7) for x in NEW])
print('lam old %.4f new %.4f; ratio of mean 1/lam %.4f, of mean 1/lam^2 %.4f' % (lo, ln, iln/ilo, np.mean([math.log(x+5e7)**-2 for x in NEW])/np.mean([math.log(x+5e7)**-2 for x in OLD])))
MEASOLD = {3: (0.9549, 0.0062), 5: (0.9844, 0.0019), 7: (0.9859, 0.0016)}   # law.py MEAS: mean, max(block, scatter)
def cell(q, k, H):
    a = A['win_%d_%d_%d' % (q, k, H)]; return a['mean'], max(a['sb'], a['ss']), np.array(a['pw'])
# comps for 'notail' at the new height (centre)
cmp = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/within/data/comp_22.6973.json'))
print('\n(1) level cells: z against law with sig_meas only and with sig_tot; comparators with sig_meas')
for q, k in ((3,1),(5,1),(7,1),(3,2),(10,2)):
    for H in (10**4, 10**5, 10**6):
        m, s, _ = cell(q, k, H); p = P['%d_%d_%d' % (5 if q == 10 else q, k, H)]
        st = math.hypot(s, p['law_err'])
        tail = cmp[str(5 if q == 10 else q)]['tail_k%d' % k][str(H)]
        print('(%2d,%d,%.0e) meas %.4f+-%.4f law %.4f+-%.4f z_meas %+5.1f z_tot %+5.1f | lead %.4f z %+5.1f | sieve %.4f z %+5.1f | law-no-tail %.4f z %+5.1f | 1: z %+5.1f'
              % (q, k, H, m, s, p['law'], p['law_err'], (m-p['law'])/s, (m-p['law'])/st, p['lead'], (m-p['lead'])/s, p['MQ'], (m-p['MQ'])/s, p['law']-tail, (m-p['law']+tail)/s, (m-1)/s))
print('\n(2) R1 pass probabilities: P(|X - law| <= 2.5 sig_tot) when X ~ N(alt, sig_meas)')
for q in (3,5,7):
    m, s, _ = cell(q, 1, 10**4); p = P['%d_1_10000' % q]; st = math.hypot(s, p['law_err']); L = p['law']
    tail = cmp[str(q)]['tail_k1']['10000']
    alts = dict(law=L, lead=p['lead'], sieve=p['MQ'], old_measured=MEASOLD[q][0], no_tail=L-tail, one=1.0, tail_x2=L+tail)
    print('  q=%d band [%.4f, %.4f]: ' % (q, L-2.5*st, L+2.5*st) + '  '.join('%s %.4f P=%.2f' % (n, v, norm.cdf((L+2.5*st-v)/s) - norm.cdf((L-2.5*st-v)/s)) for n, v in alts.items()))
print('\n(3) height contrast T')
w = RT['weights']; den = RT['den']
newm = {q: cell(q, 1, 10**4) for q in (3,5,7)}
T = sum(wq*(newm[q][0]-MEASOLD[q][0]) for wq, q in zip(w, (3,5,7)))/den
sTnew = math.sqrt(sum((wq*newm[q][1])**2 for wq, q in zip(w, (3,5,7))))/den
sTold = math.sqrt(sum((wq*MEASOLD[q][1])**2 for wq, q in zip(w, (3,5,7))))/den
print('  T = %.3f; sigma from new measurement %.3f, from the old measured values %.3f, both %.3f' % (T, sTnew, sTold, math.hypot(sTnew, sTold)))
for name, f in (('deficit ~ 1/lam', iln/ilo), ('deficit ~ 1/lam^2', np.mean([math.log(x+5e7)**-2 for x in NEW])/np.mean([math.log(x+5e7)**-2 for x in OLD]))):
    Ts = sum(wq*((1-(1-MEASOLD[q][0])*f)-MEASOLD[q][0]) for wq, q in zip(w, (3,5,7)))/den
    print('  pure scaling of the OLD measured deficits, %s (factor %.4f): T = %.3f; measured T - this = %.3f = %.1f sig_new, %.1f sig_both' % (name, f, Ts, T-Ts, (T-Ts)/sTnew, (T-Ts)/math.hypot(sTnew, sTold)))
# leading order's own change and law's own change
for name, key_new, key_old in (('leading order own change', 'lead', None), ('law own change', 'law', 'old_law')):
    if key_old is None:
        ch = {q: (1-1/((len([a for a in range(q) if math.gcd(a,q)==1])-1)*ln)) - (1-1/((len([a for a in range(q) if math.gcd(a,q)==1])-1)*lo)) for q in (3,5,7)}
    else:
        ch = {q: P['%d_1_10000' % q]['law'] - P['%d_1_10000' % q]['old_law'] for q in (3,5,7)}
    Tc = sum(wq*ch[q] for wq, q in zip(w, (3,5,7)))/den
    print('  %s: changes %s -> T = %.3f' % (name, {q: round(v, 4) for q, v in ch.items()}, Tc))
print('  per-q change new - old measured: ' + '  '.join('q=%d %+.4f +- %.4f (1/lam scaling %+.4f, law own %+.4f)' % (q, newm[q][0]-MEASOLD[q][0], math.hypot(newm[q][1], MEASOLD[q][1]), -(1-MEASOLD[q][0])*(iln/ilo-1), P['%d_1_10000' % q]['law']-P['%d_1_10000' % q]['old_law']) for q in (3,5,7)))
