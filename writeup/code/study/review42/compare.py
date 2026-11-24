import sys, json, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import mypred
T = np.load(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/tables.npz')
P = mypred.Pred((T['lnh'], T['lnPhiD'], T['lnPhiO']))
d = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/bk/predictions_per_lfunction.json'))
lgrid = np.array(d['l'])
want = [0.5, 1.0, 3.5, 7.0, 7.5]
idx = [int(np.argmin(abs(lgrid - w))) for w in want]
Ds = [401, 409, -419, -431, -404, 412, -420, 476, -408, 408, 424, -424, 440, -440, -435, 465, 645, -555, -616, 616, -455, 533, 680, -680]
B = d['B']
maxd = 0
for D in Ds:
    i = B['D'].index(D)
    bp = mypred.badprimes(D)
    for mode in ('full', 'diag'):
        mine8 = mypred.band_average(P, lgrid[idx], D, 100, 600, nodes=8, mode=mode)
        mine16 = mypred.band_average(P, lgrid[idx], D, 100, 600, nodes=16, mode=mode)
        theirs = np.array(B[mode][i])[idx]
        diff = mine16 - theirs
        maxd = max(maxd, np.abs(diff).max())
        print(f"{D:5d} {str(bp):14s} {mode:4s} theirs " + ' '.join(f'{x:.5f}' for x in theirs) +
              " | mine16-theirs " + ' '.join(f'{x:+.1e}' for x in diff) + " | 16-8 " + ' '.join(f'{x:+.0e}' for x in mine16 - mine8))
    # off part comparison
print('max abs diff', maxd)
# also full - diag (off part) comparisons
print('off-part (full-diag) comparison:')
for D in Ds:
    i = B['D'].index(D)
    th = np.array(B['full'])[i][idx] - np.array(B['diag'])[i][idx]
    mi = mypred.band_average(P, lgrid[idx], D, 100, 600, nodes=16, mode='full') - mypred.band_average(P, lgrid[idx], D, 100, 600, nodes=16, mode='diag')
    print(f"{D:5d} theirs " + ' '.join(f'{x:+.5f}' for x in th) + ' | mine ' + ' '.join(f'{x:+.5f}' for x in mi))
