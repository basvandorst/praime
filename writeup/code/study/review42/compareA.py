import sys, json, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import mypred
T = np.load(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/tables.npz')
P = mypred.Pred((T['lnh'], T['lnPhiD'], T['lnPhiO']))
d = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/bk/predictions_per_lfunction.json'))
lg = np.array(d['l']); idx = [int(np.argmin(abs(lg - w))) for w in (0.5, 1, 3.5, 7, 7.5)]
A = d['A']; mx = 0
for N in (101, 115, 124, 129, 203, 236, 165 if 165 in A['N'] else 305, 393, 172, 335):
    i = A['N'].index(N)
    for mode in ('full', 'diag'):
        mine = mypred.band_average(P, lg[idx], N, 100, 500, nodes=16, mode=mode)
        diff = mine - np.array(A[mode][i])[idx]; mx = max(mx, abs(diff).max())
print('test A: max |mine - frozen| over 10 conductors x 5 lengths x 2 laws = %.1e' % mx)
