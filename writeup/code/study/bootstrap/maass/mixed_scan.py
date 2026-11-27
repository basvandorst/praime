# Round 171: sequential bootstrap.  Usage: mixed_scan.py <eps> <X> <r_lo> <r_hi> <step> [truth]
# For each r: convex problem (L(f) with the central zero if eps = 1, Sym^2) plus L(f x g) for every true form g with
# r_g <= r - 0.05.  With "truth": instead of solving, the true forms at their exact r with their true angles (a check).
import sys, json, os, time, numpy as np
from maass import Problem, types
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'maass')
key = [f for f in json.load(open(os.path.join(D, 'lmfdb_key.json'))) if f['ap']]
def build(r, eps, X):
    T = types(r, eps, 2); used = []
    for g in key:
        if g['r'] > r - 0.05 + 1e-9: continue
        e = 1 if g['parity'] == 'odd' else 0; d = (eps + e) % 2; rg = g['r']; used.append(rg)
        thg = {int(p): np.arccos(a/2) for p, a in g['ap'].items()}
        fn = lambda ph: ph
        fn.pmp = (lambda thg: lambda p, m, th: 2*np.cos(m*th)*2*np.cos(m*thg[p]))(thg)
        T.append(('x%.2f' % rg, [d + 1j*(r + rg), d - 1j*(r + rg), d + 1j*(r - rg), d - 1j*(r - rg)], 2, fn))
    N = max(16, int(np.ceil(np.log(X)/2*(2*r + 12 + (max(used) if used else 0))/np.pi)))
    return T, used, N
eps = int(sys.argv[1]); X = float(sys.argv[2])
if 'truth' in sys.argv:
    for f in key:
        if (f['parity'] == 'odd') != (eps == 1): continue
        T, used, N = build(f['r'], eps, X); P = Problem(f['r'], eps, X, T=T, center=(eps == 1), N=N, J=3)
        th = np.array([np.arccos(f['ap'][str(p)]/2) for p in P.primes]); w = np.zeros((len(P.primes), 3)); out = []
        for A, x, c in zip(P.A, P.T, P.c):
            cf = x[3]; coef = np.array([(cf.pmp(p, m, th[P.primes.index(p)]) if hasattr(cf, 'pmp') else cf(m*th[P.primes.index(p)]))*np.log(p)*p**(-m/2) for p, m in P.pp])
            out.append('%s %+.1e' % (x[0], np.linalg.eigvalsh(A - np.einsum('i,ijk->jk', coef, P.Spm))[0]/(2*P.B.a*x[2])))
        print('true form r = %.4f (%s), %d known forms below: %s' % (f['r'], f['parity'], len(used), '  '.join(out)), flush=True)
    sys.exit()
lo, hi, step = map(float, sys.argv[3:6])
out = os.path.join(D, 'mixed_eps%d_X%g.jsonl' % (eps, X)); done = set()
if os.path.exists(out):
    for line in open(out): done.add(round(json.loads(line)['r'], 4))
for r in [round(lo + i*step, 4) for i in range(int(round((hi - lo)/step)) + 1)]:
    if r in done: continue
    t0 = time.time(); T, used, N = build(r, eps, X); P = Problem(r, eps, X, T=T, center=(eps == 1), N=N); lam, t, rd = P.solve()
    rec = dict(r=r, eps=eps, X=X, N=N, known=used, t_low=lam, t_lp=t, rounds=rd, secs=round(time.time() - t0, 1), per_form=P.mineigs(P.w),
               names=[x[0] for x in P.T], primes=P.primes, mean_ap=[float(w @ (2*np.cos(P.th))) for w in P.w])
    with open(out, 'a') as fh: fh.write(json.dumps(rec) + '\n')
