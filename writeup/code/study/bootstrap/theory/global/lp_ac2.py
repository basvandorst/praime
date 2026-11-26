import numpy as np, sys
sys.path.insert(0, '/tmp/claude-1000/prover_boot')
from lp_ac import run, bfun, cellint
for x0 in [1.1, 1.3, 1.5, 1.7]:
    for (U, n) in [(np.log(10), 60), (np.log(10), 120), (np.log(30), 120)]:
        r, edges, ts, bs, base, M = run(x0, U, n, 50.0, 60.0, 0.02)
        val = -r.fun if r.status == 0 else float('nan')
        print("x0=%.1f U=%.2f n=%d  max min(b-F) = %.6f" % (x0, np.exp(U), n, val)); sys.stdout.flush()
