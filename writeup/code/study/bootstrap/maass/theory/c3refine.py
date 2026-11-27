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
def chi(n): r = n % 3; return 0 if r == 0 else (1 if r == 1 else -1)
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
    cand = [mp.mpf(-1), mp.mpf(1)]
    der = [m*a[m - 1] for m in range(len(a), 0, -1)]
    while len(der) > 1 and der[0] == 0: der = der[1:]
    if len(der) > 1:
        try:
            for r in mp.polyroots(der, maxsteps=300, extraprec=300):
                if abs(mp.im(r)) < mp.mpf(10)**-60 and -1 < mp.re(r) < 1: cand.append(mp.re(r))
        except mp.libmp.NoConvergence:
            cand += [mp.mpf(i)/4000 - 1 for i in range(8001)]
    return min((peval(a, x), x) for x in cand)

def defect(beta, cf):
    tot = mp.mpf(0); per = {}
    for p in primes:
        a = polyx(cf[p], beta); Pc = peval(a, mp.mpf(chi(p))); mn, x = minP(a)
        per[p] = (Pc - mn, x, a); tot += Pc - mn
    return tot, per

def refine(D, beta, log=print, rounds=40):
    t0 = time.time(); Cm, G, G0, cf = setup(D); nb = D + 1; nP = len(primes)
    beta = [mp.mpf(b) for b in beta] + [mp.mpf(0)]*(nb - len(beta))
    g0 = mp.fsum(G0[i]*beta[i] for i in range(nb)); beta = [b/g0 for b in beta]
    Def, per = defect(beta, cf)
    log('deg %d start: Defect/g(0) = %s' % (2*D, mp.nstr(Def, 6)))
    xs = {p: set([mp.mpf(i)/50 - 1 for i in range(101)] + [per[p][1]]) for p in primes}
    tg = [mp.mpf(i)/20 for i in range(1201)] + [60*mp.mpf(1.03)**i for i in range(1, 200)]
    Bt = [[mp.fsum(Cm[i][k]*t**(2*k) for k in range(nb)) for i in range(nb)] for t in tg]
    Bw = [1/max(abs(v) for v in row) for row in Bt]
    TR = 1e4; eps = Def
    for rd in range(rounds):
        # rows: sig_p - P_p(x; d) <= r_p(x)/eps,   r_p(x) = P_p(x; beta) - P_p(chi_p; beta)
        A, b = [], []
        for ip, p in enumerate(primes):
            a = per[p][2]; Pc = peval(a, mp.mpf(chi(p)))
            for x in xs[p]:
                rhs = (peval(a, x) - Pc)/eps
                if rhs > 1e6: continue
                row = np.zeros(nb + nP)
                row[:nb] = [-float(mp.fsum(cf[p][m - 1][i]*x**m for m in range(1, len(cf[p]) + 1))) for i in range(nb)]
                row[nb + ip] = 1; A.append(row); b.append(float(rhs))
        for t, row, w in zip(tg, Bt, Bw):
            Bv = mp.fsum(row[i]*beta[i] for i in range(nb))
            r_ = np.zeros(nb + nP); r_[:nb] = [-float(v*w) for v in row]; A.append(r_); b.append(float(Bv*w/eps))
        obj = np.zeros(nb + nP)
        for ip, p in enumerate(primes):
            x = chi(p); obj[:nb] += [float(mp.fsum(cf[p][m - 1][i]*mp.mpf(x)**m for m in range(1, len(cf[p]) + 1))) for i in range(nb)]
            obj[nb + ip] = -1
        Aeq = np.concatenate([[float(v) for v in G0], np.zeros(nP)])[None, :]
        res = linprog(obj, A_ub=np.array(A), b_ub=np.array(b), A_eq=Aeq, b_eq=[0.0],
                      bounds=[(-TR, TR)]*nb + [(None, None)]*nP, method='highs',
                      options={'primal_feasibility_tolerance': 1e-9, 'dual_feasibility_tolerance': 1e-9})
        if res.status != 0:
            log('  round %d: LP status %d' % (rd, res.status)); eps *= 10; continue
        d = [mp.mpf(float(v)) for v in res.x[:nb]]
        hit = max(abs(float(v)) for v in res.x[:nb]) > 0.98*TR
        nbeta = [beta[i] + eps*d[i] for i in range(nb)]
        g0 = mp.fsum(G0[i]*nbeta[i] for i in range(nb)); nbeta = [v/g0 for v in nbeta]
        nDef, nper = defect(nbeta, cf)
        added = 0
        for p in primes:
            x = nper[p][1]
            if x not in xs[p]: xs[p].add(x); added += 1
        if nDef < Def:
            beta, Def, per = nbeta, nDef, nper
            log('  round %2d: Defect/g(0) = %s (LP model %.3e, cuts +%d, TR %.0e%s)  [%.0f s]' % (rd, mp.nstr(Def, 8), res.fun*float(eps), added, TR, ' hit' if hit else '', time.time() - t0))
            eps = max(Def, mp.mpf(10)**-200)
            if hit: TR *= 30
        else:
            log('  round %2d: no gain (%s), cuts +%d, TR %.0e' % (rd, mp.nstr(nDef, 4), added, TR))
            if added == 0 and TR < 1e2: break
            if added == 0: TR /= 30
    bcoef = [mp.fsum(beta[i]*Cm[i][k] for i in range(nb)) for k in range(nb)]
    top = sorted(per.items(), key=lambda kv: -kv[1][0])[:5]
    log('deg %d: Defect/g(0) = %s, bound %s; carriers %s' % (2*D, mp.nstr(Def, 10), mp.nstr(3*mp.exp(-2*Def), 25),
        ', '.join('%d:%s@%s' % (p, mp.nstr(v[0], 3), mp.nstr(v[1], 5)) for p, v in top)))
    gl = ', '.join('%d:%s' % (n, mp.nstr(mp.fsum(G[n][i]*beta[i] for i in range(nb)), 3)) for n in sorted(G)[:22])
    log('   g(log n)/g(0): ' + gl)
    return beta, Def, dict(D=D, beta=[mp.nstr(v, 60) for v in beta], b=[mp.nstr(v, 80) for v in bcoef], defect=mp.nstr(Def, 20),
                           per={str(p): [mp.nstr(v[0], 6), mp.nstr(v[1], 10)] for p, v in per.items() if v[0] > mp.mpf(10)**-60})

if __name__ == '__main__':
    st = json.load(open(sys.argv[2])); beta = st['beta']
    for D in range(int(sys.argv[3]), int(sys.argv[4]) + 1):
        beta, Def, r = refine(D, beta)
        json.dump(r, open(OUT + 'c3ref_d%d.json' % (2*D), 'w'), indent=1)
        beta = r['beta']
