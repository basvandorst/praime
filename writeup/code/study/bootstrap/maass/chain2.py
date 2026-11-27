# Round 173: the chain of round 172 continued above r = 20 (blind prediction; writeup/data/study/preregister/round173.txt).
# Usage: chain2.py <X> <r_hi>
import sys, json, os, time, numpy as np
from maass import Problem, types
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'maass')
X = float(sys.argv[1]); RHI = float(sys.argv[2])
found = json.load(open(os.path.join(D, 'chain_X%g_found.json' % X)))
for g in found: g['ap'] = {int(p): a for p, a in g['ap'].items()}
n0 = len(found)
def solve(r, eps, J=61):
    T = types(r, eps, 2); used = [g for g in found if g['r'] <= r - 0.05 + 1e-9]
    for g in used:
        d = (eps + g['eps']) % 2; rg = g['r']; thg = {p: np.arccos(np.clip(a/2, -1, 1)) for p, a in g['ap'].items()}
        fn = lambda ph: ph
        fn.pmp = (lambda thg: lambda p, m, th: 2*np.cos(m*th)*2*np.cos(m*thg[p]))(thg)
        T.append(('x%.2f' % rg, [d + 1j*(r + rg), d - 1j*(r + rg), d + 1j*(r - rg), d - 1j*(r - rg)], 2, fn))
    N = max(16, int(np.ceil(np.log(X)/2*(2*r + 12 + (max(g['r'] for g in used) if used else 0))/np.pi)))
    P = Problem(r, eps, X, T=T, center=(eps == 1), N=N, J=J); lam, t, rd = P.solve()
    return t, {p: float(w @ (2*np.cos(P.th))) for p, w in zip(P.primes, P.w)}, len(used)
def record(r0, eps):
    """r0: first feasible grid point of a piece (the grid point before it is infeasible)."""
    cand = None
    for k in range(3):
        if solve(r0 + 0.01*k, eps, J=121)[0] >= -0.0005: cand = r0 + 0.01*k; break
    if cand is not None:
        lo, hi = r0 - 0.01, cand
        for _ in range(11):
            mid = (lo + hi)/2
            if solve(mid, eps, J=121)[0] >= -0.0005: hi = mid
            else: lo = mid
        rh = hi; how = 'left crossing'
    else:
        a, b = r0 - 0.01, r0 + 0.01
        for _ in range(12):
            m1, m2 = a + (b - a)/3, b - (b - a)/3
            if solve(m1, eps, J=121)[0] < solve(m2, eps, J=121)[0]: a = m1
            else: b = m2
        rh = (a + b)/2; how = 'maximum'
    t, ap, nk = solve(rh, eps, J=121)
    g = dict(r=rh, eps=eps, ap=ap, t=t, how=how, r0=r0, known=nk); found.append(g)
    print('FOUND %s form: r^ = %.5f (%s; t* = %+.6f; %d fed-back forms); a^_p: %s' % ('odd' if eps else 'even', rh, how, t, nk, ' '.join('%d:%+.4f' % kv for kv in ap.items())), flush=True)
rng = lambda a, b: [round(a + 0.01*i, 2) for i in range(int(round((b - a)/0.01)) + 1)]
pts = [(r, e) for r in rng(20.01, RHI) for e in (1, 0)]
pts.sort(); prev = {0: False, 1: False}; out = os.path.join(D, 'chain2_X%g.jsonl' % X); open(out, 'w').close(); t00 = time.time()
for r, eps in pts:
    t, ap, nk = solve(r, eps); feas = t >= -0.002
    with open(out, 'a') as fh: fh.write(json.dumps(dict(r=r, eps=eps, t_lp=t, known=nk, ap=ap)) + '\n')
    if feas and not prev[eps]: record(r, eps)
    prev[eps] = feas
json.dump(found[n0:], open(os.path.join(D, 'chain2_X%g_found.json' % X), 'w'), indent=1)
print('--- %d forms predicted in (20, %g]; %.0f s' % (len(found) - n0, RHI, time.time() - t00))
