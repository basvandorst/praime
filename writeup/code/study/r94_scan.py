"""Round 94's selection scan (committed after the twenty-ninth review): R_zero(X, 0.1) from r93.py at 120 values of X,
log-spaced from 1.2e9 to 4e10; the three registered X were the maximum (3.900e9), the minimum (3.268e9) and 2e10.
usage: r94_scan.py"""
import os, json, math
import numpy as np
import r93

xs = np.exp(np.linspace(math.log(1.2e9), math.log(4e10), 120))
vals = [r93.r_zero(float(X), 0.1)['R_zero'] for X in xs]
out = {'X': xs.tolist(), 'R_zero_delta_0.1': vals, 'argmax_X': float(xs[int(np.argmax(vals))]), 'argmin_X': float(xs[int(np.argmin(vals))])}
print(out['argmax_X'], max(vals), out['argmin_X'], min(vals))
json.dump(out, open(os.path.join(r93.OUT, 'r94_scan.json'), 'w'), indent=1)
