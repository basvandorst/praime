import numpy as np
z=np.load(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/data/study/r144/fresh_zeros.npz')
tot=0
for k in z.keys():
    g=np.sort(z[k]); d=np.diff(g)
    for thr in (1e-3,1e-2,0.03):
        pass
    i=np.nonzero(d<1e-3)[0]
    if len(i): print(k, len(g), [(round(g[j],9), '%.2e'%d[j]) for j in i]); tot+=1
    m=d.min()
print('curves with near-duplicates', tot)
print('smallest gaps overall', sorted((np.diff(np.sort(z[k])).min(),k) for k in z.keys())[:6])
