#!/usr/bin/env python3
"""t1_power.py -- expected pooled se of the T1 cells and the pilot-calibrated lower-order alternative H_lo."""
import json, math
import numpy as np
cells = json.load(open('reg/predictions_frozen.json'))['cells']
Lp = math.log(10007)
for m, sl in ((2, (1.0,)), (3, (0.7, 0.8, 0.9, 1.0)), (4, (0.5, 0.6, 0.7))):
    for s in sl:
        cl = [c for c in cells.values() if abs(c['s'] - s) < 1e-9]
        w = np.array([1 / c[f'se{m}_exp'] ** 2 for c in cl])
        se = 1 / math.sqrt(w.sum())
        Spil = cl[0][f'Dpilot{m}'] + cl[0][f'U{m}']
        Hlo = np.sum(w * np.array([Spil * Lp / c['logq'] for c in cl])) / w.sum()
        g6 = [c for c in cl if c['q'] > 500000]; g3 = [c for c in cl if c['q'] < 50000]
        se6 = 1 / math.sqrt(sum(1 / c[f'se{m}_exp'] ** 2 for c in g6)); se3 = 1 / math.sqrt(sum(1 / c[f'se{m}_exp'] ** 2 for c in g3))
        print(f'T1 m={m} s={s}: U={cl[0][f"U{m}"]:+.6f} pooled se={se:.6f} | pilot S={Spil:+.6f} -> H_lo pooled={Hlo:+.6f} (z={Hlo/se:.1f}) | se(F3)={se3:.6f} se(F6)={se6:.6f}')
