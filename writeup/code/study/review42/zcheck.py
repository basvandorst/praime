import numpy as np, json, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
from mydata import theta_count, unfold, LZ
d = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/bk/predictions_per_lfunction.json'))
z = np.load(LZ + 'r134_B.npz')
worst = []
for D in d['B']['D']:
    g = np.sort(z['D%d' % D]); a = 0 if D > 0 else 1
    mid = 0.5 * (g[1:] + g[:-1])
    dev = np.arange(1, len(g)) - theta_count(mid, abs(D), a)    # = S(t) + const ~ mean 0 between zeros
    sel = (mid > 100)
    blocks = [dev[(mid > lo) & (mid <= lo + 50)].mean() for lo in range(100, 600, 50)]
    worst.append((min(blocks), max(blocks), D, float(np.diff(unfold(g[g > 100], abs(D))).min()), len(g)))
worst.sort()
print('most negative 50-unit block mean of S(t):', [(round(w[0], 2), w[2]) for w in worst[:5]])
print('most positive:', sorted([(round(w[1], 2), w[2]) for w in worst])[-5:])
print('smallest unfolded gaps:', sorted([(round(w[3], 4), w[2]) for w in worst])[:5])
