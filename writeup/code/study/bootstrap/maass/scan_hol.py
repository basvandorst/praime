# Controls of rounds 165 and 168: holomorphic weight kappa as a real parameter (level 1).  Usage: scan_hol.py <X> <k_lo> <k_hi> <step> [center]
import sys, json, os, time, numpy as np
from maass import Problem
X = float(sys.argv[1]); lo, hi, step = map(float, sys.argv[2:5]); center = 'center' in sys.argv; N = None
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'maass', 'hol_X%g%s.jsonl' % (X, '_center' if center else ''))
done = set()
if os.path.exists(out):
    for line in open(out): done.add(round(json.loads(line)['r'], 4))
for kap in [round(lo + i*step, 4) for i in range(int(round((hi - lo)/step)) + 1)]:
    if kap in done: continue
    t0 = time.time(); P = Problem(kap, 0, X, hol=True, N=N, center=center); lam, t, rd = P.solve()
    rec = dict(r=kap, X=X, N=P.N, t_low=lam, t_lp=t, rounds=rd, secs=round(time.time() - t0, 1), per_form=P.mineigs(P.w),
               primes=P.primes, mean_ap=[float(w @ (2*np.cos(P.th))) for w in P.w])
    with open(out, 'a') as fh: fh.write(json.dumps(rec) + '\n')
