import numpy as np, json, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
from analysis import *
Dall = [D for D in Ds if D != 1153]
med = np.median([abs(D) for D in Dall])
subsets = {'D>0': [D for D in Dall if D > 0], 'D<0': [D for D in Dall if D < 0],
           '|D|<=med': [D for D in Dall if abs(D) <= med], '|D|>med': [D for D in Dall if abs(D) > med]}
print('median |D|', med)
for nm, S in subsets.items():
    a = arrays(S); r = readings(a); se = boot(a, 600, seed=7)
    z = r['m']/se[5:]
    print(f"{nm:9s} n={len(S):3d} T2 p={r['pT']:.3f} max|m/se|={np.max(np.abs(z)):.2f} at l={ls[np.argmax(np.abs(z))]}  "
          f"gA={r['gA']:.4f}+-{se[0]:.4f} gQ={r['gQ']:.4f}+-{se[1]:.4f} goff={r['goff']:.4f}+-{se[2]:.4f} g2={r['g2']:.5f}+-{se[3]:.5f}")
    print('          chi2 p=2,3,5 vs law', [round(r['chi2'][p]['full'], 1) for p in (2, 3, 5)])
