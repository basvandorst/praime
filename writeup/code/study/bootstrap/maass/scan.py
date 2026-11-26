# Resumable scan of t*(r) for round 165.  Usage: scan.py <eps> <X> <r_lo> <r_hi> <step> [kmax] [J] [out]
# Appends one JSON line per point to the output file and skips points already there.  With step >= 0.1 the multiples of 0.5 come first.
import sys, json, os, time, numpy as np
from maass import Problem
eps = int(sys.argv[1]); X = float(sys.argv[2]); lo, hi, step = map(float, sys.argv[3:6])
kmax = int(sys.argv[6]) if len(sys.argv) > 6 else 4; J = int(sys.argv[7]) if len(sys.argv) > 7 else 61
out = sys.argv[8] if len(sys.argv) > 8 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'maass', 'scan_eps%d_X%g_k%d_J%d.jsonl' % (eps, X, kmax, J))
os.makedirs(os.path.dirname(out), exist_ok=True)
done = set()
if os.path.exists(out):
    for line in open(out): done.add(round(json.loads(line)['r'], 4))
rs = [round(lo + i*step, 4) for i in range(int(round((hi - lo)/step)) + 1)]
if step >= 0.0999: rs = [r for r in rs if abs(r*2 - round(r*2)) < 1e-9] + [r for r in rs if abs(r*2 - round(r*2)) >= 1e-9]
for r in rs:
    if r in done: continue
    t0 = time.time(); P = Problem(r, eps, X, J=J, kmax=kmax)
    lam, t, rd = P.solve()
    rec = dict(r=r, eps=eps, X=X, kmax=kmax, J=J, N=P.N, t_low=lam, t_lp=t, rounds=rd, secs=round(time.time() - t0, 1),
               per_form=P.mineigs(P.w), primes=P.primes, mean_ap=[float(w @ (2*np.cos(P.th))) for w in P.w],
               atoms=[[(round(float(P.th[j]), 5), round(float(w[j]), 5)) for j in np.nonzero(w > 1e-6)[0]] for w in P.w])
    with open(out, 'a') as fh: fh.write(json.dumps(rec) + '\n')
