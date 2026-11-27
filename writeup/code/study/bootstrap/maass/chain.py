# Round 172: the self-contained chain (see writeup/data/study/preregister/round172.txt).  Usage: chain.py <X>
import sys, json, os, time, numpy as np
from maass import Problem, types
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'maass')
X = float(sys.argv[1]); found = []          # dicts: r, eps, ap {p: a_p}
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
    t0, ap0, _ = solve(r0, eps, J=121)
    if t0 >= -0.0005:
        lo, hi = r0 - 0.01, r0
        for _ in range(10):
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
pts = [(r, 1) for a, b in ((9.50, 9.58), (12.13, 12.23), (14.31, 14.41), (15.90, 20.00)) for r in rng(a, b)] + [(r, 0) for a, b in ((13.74, 13.82), (17.65, 17.98), (19.28, 19.47)) for r in rng(a, b)]
pts.sort(); prev = {0: False, 1: False}; out = os.path.join(D, 'chain_X%g.jsonl' % X); open(out, 'w').close(); t00 = time.time()
for r, eps in pts:
    t, ap, nk = solve(r, eps); feas = t >= -0.002
    with open(out, 'a') as fh: fh.write(json.dumps(dict(r=r, eps=eps, t_lp=t, known=nk, ap=ap)) + '\n')
    if feas and not prev[eps]: record(r, eps)
    prev[eps] = feas
json.dump(found, open(os.path.join(D, 'chain_X%g_found.json' % X), 'w'), indent=1)
key = [f for f in json.load(open(os.path.join(D, 'lmfdb_key.json'))) if f['r'] < 20]
print('--- comparison with the key (%d forms recorded, %d true forms below 20; %.0f s)' % (len(found), len(key), time.time() - t00))
for g in found:
    c = [f for f in key if (f['parity'] == 'odd') == (g['eps'] == 1)]; f = min(c, key=lambda f: abs(f['r'] - g['r']))
    print('r^ = %.5f %s: nearest true %.5f (difference %+.5f); a_2 %+.4f (true %+.4f), a_3 %+.4f (true %+.4f), a_5 %+.4f (%+.4f), a_7 %+.4f (%+.4f), a_11 %+.4f (%+.4f), a_13 %+.4f (%+.4f)' % (
        g['r'], 'odd ' if g['eps'] else 'even', f['r'], g['r'] - f['r'], g['ap'][2], f['ap']['2'], g['ap'][3], f['ap']['3'], g['ap'][5], f['ap']['5'], g['ap'][7], f['ap']['7'], g['ap'][11], f['ap']['11'], g['ap'][13], f['ap']['13']))
