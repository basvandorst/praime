# Round 167: scan of t*(l1, l2) on SL(3,Z).  Usage: gl3_scan.py <X> local|coarse     (resumable; appends JSON lines)
import sys, json, os, time, numpy as np
from gl3 import GL3
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'maass')
X = float(sys.argv[1]); mode = sys.argv[2]
out = os.path.join(D, 'gl3_%s_X%g.jsonl' % (mode, X))
key = sorted(json.load(open(os.path.join(D, 'gl3_key.json')))['forms'], key=lambda f: f['laplace'])
pts = []
if mode == 'local':
    for f in key[:5]:
        l1, l2 = f['lambda']; pts.append((l1, l2, 'exact'))
        c1, c2 = round(l1, 1), round(l2, 1)
        for i in range(-4, 5):
            for j in range(-4, 5):
                a, b = round(c1 + 0.1*i, 4), round(c2 + 0.1*j, 4)
                if 0 <= b <= a: pts.append((a, b, 'grid'))
else:
    for i in range(41):
        for j in range(i + 1): pts.append((0.5*i, 0.5*j, 'grid'))
done = set()
if os.path.exists(out):
    for line in open(out): d = json.loads(line); done.add((round(d['l1'], 6), round(d['l2'], 6)))
for l1, l2, kind in pts:
    if (round(l1, 6), round(l2, 6)) in done: continue
    done.add((round(l1, 6), round(l2, 6)))
    t0 = time.time(); P = GL3(l1, l2, X); lam, t, rd = P.solve()
    rec = dict(l1=l1, l2=l2, kind=kind, X=X, N=P.N, t_low=lam, t_lp=t, rounds=rd, secs=round(time.time() - t0, 1), per_form=P.per_form,
               primes=P.primes, a_p=[[z.real, z.imag] for z in P.a_p])
    with open(out, 'a') as fh: fh.write(json.dumps(rec) + '\n')
