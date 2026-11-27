# Iterative refinement of the LP of c3opt.py: variables centred at the current point beta and scaled by eps (the current
# defect), so that float64 HiGHS resolves improvements far below 1e-10.  All data (rows, right-hand sides, the defect) are
# computed in mpmath at 250 digits.  Usage: c3refine.py BASISFILE START.json D_from D_to   (degrees 2*D)
# Writes c3ref_d<2D>.json (same format as c3opt).
import sys; sys.dont_write_bytecode = True
import json, time
import numpy as np, mpmath as mp
from scipy.optimize import linprog
mp.mp.dps = 250
basis = json.load(open(sys.argv[1])); NMAX = basis['NMAX']
mom = [mp.mpf(a) for a, r in basis['moments']]
vals = {int(n): [mp.mpf(a) for a, r in v] for n, v in basis['values'].items()}
def chi(n): return (0, 1, -1, -1, 1)[n % 5]
primes = [p for p in range(2, NMAX + 1) if all(p % d for d in range(2, int(p**0.5) + 1))]
pw = {p: [p**m for m in range(1, 40) if p**m <= NMAX] for p in primes}
OUT = '/home/bas/praime/writeup/code/study/bootstrap/maass/theory/'

def setup(D):
    M = mp.matrix(D + 1, D + 1)
    for j in range(D + 1):
        for k in range(D + 1): M[j, k] = mom[j + k]
    Ci = mp.inverse(mp.cholesky(M))
    Cm = [[Ci[i, k] for k in range(D + 1)] for i in range(D + 1)]
    G = {n: [mp.fsum(Cm[i][k]*vals[n][k] for k in range(D + 1)) for i in range(D + 1)] for n in vals}
    G0 = [mp.fsum(Cm[i][k]*mom[k] for k in range(D + 1)) for i in range(D + 1)]
    cf = {p: [[mp.log(p)*mp.power(p, -mp.mpf(m)/2)*G[p**m][i] for i in range(D + 1)] for m in range(1, len(pw[p]) + 1)] for p in primes}
    return Cm, G, G0, cf

def polyx(cfp, beta): return [mp.fsum(r[i]*beta[i] for i in range(len(beta))) for r in cfp]
def peval(a, x): return mp.fsum(a[m - 1]*x**m for m in range(1, len(a) + 1))
def minP(a):
    """min over [-1, 1] of sum_m a[m-1] x^m: endpoints and every real critical point, found from sign changes of P'
    on a 4001-point grid (float scan, then bisection and Newton in mpmath) and from polyroots (loose filter, Newton-polished)."""
    M = len(a)
    d1 = lambda x: mp.fsum(m*a[m - 1]*x**(m - 1) for m in range(1, M + 1))
    d2 = lambda x: mp.fsum(m*(m - 1)*a[m - 1]*x**(m - 2) for m in range(2, M + 1))
    cand = [mp.mpf(-1), mp.mpf(1)]
    af = np.array([float(v) for v in a]); xg = np.linspace(-1, 1, 4001)
    dg = sum((m + 1)*af[m]*xg**m for m in range(M))
    for i in np.nonzero(np.sign(dg[:-1])*np.sign(dg[1:]) <= 0)[0]:
        lo, hi = mp.mpf(xg[i]), mp.mpf(xg[i + 1]); flo = d1(lo)
        for _ in range(80):
            mid = (lo + hi)/2; fm = d1(mid)
            if (fm > 0) == (flo > 0): lo, flo = mid, fm
            else: hi = mid
        cand.append((lo + hi)/2)
    if M > 1:
        try:
            der = [m*a[m - 1] for m in range(M, 0, -1)]
            while len(der) > 1 and der[0] == 0: der = der[1:]
            for r in mp.polyroots(der, maxsteps=400, extraprec=400):
                if abs(mp.im(r)) < 1e-6 and -1.01 < mp.re(r) < 1.01:
                    x = mp.re(r)
                    for _ in range(30):
                        h = d2(x)
                        if h == 0: break
                        x = x - d1(x)/h
                    if -1 < x < 1: cand.append(x)
        except Exception:
            pass
    return min((peval(a, x), x) for x in cand)

def defect(beta, cf):
    tot = mp.mpf(0); per = {}
    for p in primes:
        a = polyx(cf[p], beta); Pc = peval(a, mp.mpf(chi(p))); mn, x = minP(a)
        per[p] = (Pc - mn, x, a); tot += Pc - mn
    return tot, per

PLP = 150   # primes and prime powers used in the LP rows (all n <= NMAX enter the exact defect)
def refine(D, beta, log=print, rounds=30):
    """LP in SVD coordinates: W maps beta to (g(0), coefficients log p p^{-m/2} g(log p^m), p^m <= PLP); W = U S V^T;
    the step is eps * V S^{-1} y, so every LP coefficient is an entry of U (size <= 1)."""
    t0 = time.time(); Cm, G, G0, cf = setup(D); nb = D + 1
    lp_pr = [p for p in primes if p <= PLP]; nP = len(lp_pr)
    keys = [(p, m) for p in lp_pr for m in range(1, len(pw[p]) + 1) if p**m <= PLP]
    W = mp.matrix(1 + len(keys), nb)
    for i in range(nb): W[0, i] = G0[i]
    for r, (p, m) in enumerate(keys):
        for i in range(nb): W[r + 1, i] = cf[p][m - 1][i]
    U, S, Vt = mp.svd_r(W)
    rank = sum(1 for j in range(nb) if S[j] > S[0]*mp.mpf(10)**-200)
    T = [[Vt[j, i]/S[j] for j in range(rank)] for i in range(nb)]          # beta-step = T y
    Uf = {k: [float(U[r + 1, j]) for j in range(rank)] for r, k in enumerate(keys)}
    U0 = [float(U[0, j]) for j in range(rank)]
    log('deg %d: singular values %s .. %s, rank %d' % (2*D, mp.nstr(S[0], 3), mp.nstr(S[rank - 1], 3), rank))
    beta = [mp.mpf(b) for b in beta] + [mp.mpf(0)]*(nb - len(beta))
    g0 = mp.fsum(G0[i]*beta[i] for i in range(nb)); beta = [b/g0 for b in beta]
    Def, per = defect(beta, cf)
    log('deg %d start: Defect/g(0) = %s' % (2*D, mp.nstr(Def, 6)))
    xs = {p: set([mp.mpf(i)/50 - 1 for i in range(101)] + [per[p][1]]) for p in lp_pr}
    # B >= 0: sampled rows B(t; beta + eps T y) >= 0 and the leading coefficient >= 0, in the y coordinates
    tg = [mp.mpf(i)/20 for i in range(1201)] + [60*mp.mpf(1.03)**i for i in range(1, 260)]
    Bc = [[mp.fsum(Cm[i][k]*t**(2*k) for k in range(nb)) for i in range(nb)] for t in tg] + [[Cm[i][D] for i in range(nb)]]
    BY = []
    for row in Bc:
        ry = [mp.fsum(row[i]*T[i][j] for i in range(nb)) for j in range(rank)]
        w = 1/max(abs(v) for v in ry); BY.append((row, [float(v*w) for v in ry], w))
    TR = 1e2; eps = Def
    for rd in range(rounds):
        A, b = [], []
        for ip, p in enumerate(lp_pr):
            a = per[p][2]; Pc = peval(a, mp.mpf(chi(p))); ms = [m for m in range(1, len(pw[p]) + 1) if p**m <= PLP]
            for x in xs[p]:
                rhs = (peval(a, x) - Pc)/eps
                if rhs > 1e7: continue
                xf = float(x); row = np.zeros(rank + nP)
                row[:rank] = -sum(np.array(Uf[(p, m)])*xf**m for m in ms)
                row[rank + ip] = 1; A.append(row); b.append(float(rhs))
        for row, ryw, w in BY:
            Bv = mp.fsum(row[i]*beta[i] for i in range(nb))
            r_ = np.zeros(rank + nP); r_[:rank] = [-v for v in ryw]; A.append(r_); b.append(float(Bv*w/eps))
        obj = np.zeros(rank + nP)
        for ip, p in enumerate(lp_pr):
            x = chi(p); ms = [m for m in range(1, len(pw[p]) + 1) if p**m <= PLP]
            obj[:rank] += sum(np.array(Uf[(p, m)])*float(x)**m for m in ms); obj[rank + ip] = -1
        Aeq = np.concatenate([U0, np.zeros(nP)])[None, :]
        res = linprog(obj, A_ub=np.array(A), b_ub=np.array(b), A_eq=Aeq, b_eq=[0.0],
                      bounds=[(-TR, TR)]*rank + [(None, None)]*nP, method='highs',
                      options={'primal_feasibility_tolerance': 1e-10, 'dual_feasibility_tolerance': 1e-10})
        if res.status != 0:
            log('  round %d: LP status %d' % (rd, res.status)); TR /= 10; continue
        y = [mp.mpf(float(v)) for v in res.x[:rank]]
        hit = max(abs(float(v)) for v in res.x[:rank]) > 0.98*TR
        nbeta = [beta[i] + eps*mp.fsum(T[i][j]*y[j] for j in range(rank)) for i in range(nb)]
        g0 = mp.fsum(G0[i]*nbeta[i] for i in range(nb)); nbeta = [v/g0 for v in nbeta]
        nDef, nper = defect(nbeta, cf)
        added = 0
        for p in lp_pr:
            x = nper[p][1]
            if x not in xs[p]: xs[p].add(x); added += 1
        if nDef < Def:
            beta, Def, per = nbeta, nDef, nper
            log('  round %2d: Defect/g(0) = %s (LP model %.3e, cuts +%d, TR %.0e%s)  [%.0f s]' % (rd, mp.nstr(Def, 8), res.fun*float(eps), added, TR, ' hit' if hit else '', time.time() - t0))
            eps = max(Def, mp.mpf(10)**-200)
            if hit: TR *= 10
        else:
            log('  round %2d: no gain (%s, LP model %.3e), cuts +%d, TR %.0e' % (rd, mp.nstr(nDef, 4), res.fun*float(eps), added, TR))
            if added == 0: TR /= 10
            if TR < 1e-3: break
    bcoef = [mp.fsum(beta[i]*Cm[i][k] for i in range(nb)) for k in range(nb)]
    top = sorted(per.items(), key=lambda kv: -kv[1][0])[:5]
    log('deg %d: Defect/g(0) = %s, bound %s; carriers %s' % (2*D, mp.nstr(Def, 10), mp.nstr(5*mp.exp(-2*Def), 25),
        ', '.join('%d:%s@%s' % (p, mp.nstr(v[0], 3), mp.nstr(v[1], 5)) for p, v in top)))
    gl = ', '.join('%d:%s' % (n, mp.nstr(mp.fsum(G[n][i]*beta[i] for i in range(nb)), 3)) for n in sorted(G)[:22])
    log('   g(log n)/g(0): ' + gl)
    return beta, Def, dict(D=D, beta=[mp.nstr(v, 60) for v in beta], b=[mp.nstr(v, 80) for v in bcoef], defect=mp.nstr(Def, 20),
                           per={str(p): [mp.nstr(v[0], 6), mp.nstr(v[1], 10)] for p, v in per.items() if v[0] > mp.mpf(10)**-60})

if __name__ == '__main__':
    st = json.load(open(sys.argv[2])); beta = st['beta']
    for D in range(int(sys.argv[3]), int(sys.argv[4]) + 1):
        beta, Def, r = refine(D, beta)
        json.dump(r, open(OUT + 'c5ref2_d%d.json' % (2*D), 'w'), indent=1)
        beta = r['beta']
