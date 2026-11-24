"""Numbers for the draft registration: model predictions CG(1e5) at the proposed windows' lambda, expected measurement
errors (taken equal to round 139's fresh errors: the per-block noise does not depend on N), approximate comparator and
first-order values at that height (scaled from lam 24.16 by (24.16/lam)^1.9, round 139's observed V scaling), and the
power of each reading.  Output reg_numbers_output.txt / reg_numbers.json."""
import json, math, numpy as np, glob
from statistics import NormalDist
N = NormalDist()
WIN = [2.40e11 + i * 1e8 for i in range(24)]
lam = float(np.mean([math.log(x + 5e7) for x in WIN]))
F = json.load(open('../../../../writeup/data/study/r139/rep_fresh.json'))['cells']
FO = json.load(open('firstorder.json'))['cells']
J = json.load(open('cg4_100000_26.21.json')); n = len(J)
HS = ('10000', '100000', '1000000')
scale = (24.16 / lam) ** 1.9
CELLS = ['7_1', '9_1', '11_1', '13_1', '11_2', '13_2', '5_1', '7_2', '9_2', '5_2']
print('proposed windows [2.40e11 + i 1e8, +1e8), i = 0..23; lam = %.3f; model CG(1e5) with %d windows at LREF 26.21' % (lam, n))
out = {'lam': lam, 'n_model': n, 'cells': {}}
print('cell   CG slope +- MC        CG W(1e4) +- MC     meas err  LR17~   first-order~   (CG-LR)/err  (CG-FO)/err  P(pass|law) P(F1 fires|law)')
for cell in CELLS:
    W = {H: np.array([np.mean(v[cell]['blk'][H]) for v in J.values()]) for H in HS}
    sl = (W['1000000'] - W['10000']) / math.log(100)
    s, ms = sl.mean(), sl.std(ddof=1) / math.sqrt(n)
    w4, mw = W['10000'].mean(), W['10000'].std(ddof=1) / math.sqrt(n)
    err = F[cell]['all']['slope'][1]; werr = F[cell]['all']['W'][0][1]
    lr = F[cell]['all']['lr'] * scale; fo = FO[cell]['s_inf'] * scale
    # reading: |s_meas - s_CG| <= 2.5 err ; under the law s_meas ~ N(s_CG_true, err), s_CG_true ~ N(s, ms)
    tot = math.hypot(err, ms); p_pass = 2 * N.cdf(2.5 * err / tot) - 1
    out['cells'][cell] = dict(cg_slope=[s, ms], cg_W1e4=[w4, mw], err=err, W_err=werr, lr=lr, fo=fo, p_pass=p_pass)
    print('%-5s  %+.5f +- %.5f   %.4f +- %.4f    %.5f  %+.5f  %+.5f       %+.1f        %+.1f        %.3f' % (cell, s, ms, w4, mw, err, lr, fo, (s - lr) / err, (s - fo) / err, p_pass))
json.dump(out, open('reg_numbers.json', 'w'), indent=1)
