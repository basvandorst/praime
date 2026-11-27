# Optimisation of Defect(g)/g(0) over B (even, degree 2D in t, B >= 0 on R) for g^ = Xi_chi^2 B, chi = chi_{-3}.
import sys; sys.dont_write_bytecode = True
# Defect(g) = sum_p [P_p(chi(p)) - min_{[-1,1]} P_p],  P_p(x) = sum_m x^m log p p^{-m/2} g(m log p).
# LP (HiGHS) in the orthonormal basis of Xi^2 dt, cutting planes in x; then an exact evaluation of the Defect in mpmath.
# Usage: c3opt.py BASISFILE D [D2 ...]      Output: c3opt_d<2D>.json (B coefficients as decimal strings) and a log on stdout.
import sys, json, time, math
import numpy as np, mpmath as mp
from scipy.optimize import linprog
mp.mp.dps = 120
basis = json.load(open(sys.argv[1]))
NMAX = basis['NMAX']
mom = [mp.mpf(a) for a, r in basis['moments']]
vals = {int(n): [mp.mpf(a) for a, r in v] for n, v in basis['values'].items()}
def chi(n): return (0, 1, -1, -1, 1)[n % 5]
primes = [p for p in range(2, NMAX + 1) if all(p % d for d in range(2, int(p**0.5) + 1))]
pw = {p: [p**m for m in range(1, 40) if p**m <= NMAX] for p in primes}

def setup(D):
    M = mp.matrix(D + 1, D + 1)
    for j in range(D + 1):
        for k in range(D + 1): M[j, k] = mom[j + k]
    L = mp.cholesky(M); Ci = mp.inverse(L)          # rows of Ci: orthonormal polynomials in t^2
    Cm = [[Ci[i, k] for k in range(D + 1)] for i in range(D + 1)]
    G = {n: [mp.fsum(Cm[i][k]*vals[n][k] for k in range(D + 1)) for i in range(D + 1)] for n in vals}
    G0 = [mp.fsum(Cm[i][k]*mom[k] for k in range(D + 1)) for i in range(D + 1)]
    return Cm, G, G0

def coeffs(p, G, D):
    """c[m-1][i] = log p p^{-m/2} G_i(p^m)."""
    return [[mp.log(p)*mp.power(p, -mp.mpf(m)/2)*G[p**m][i] for i in range(D + 1)] for m in range(1, len(pw[p]) + 1)]

def Ppoly(c, beta):
    return [mp.fsum(ci[i]*beta[i] for i in range(len(beta))) for ci in c]   # coefficients of x^1..x^M

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

def defect(beta, G, D):
    tot = mp.mpf(0); per = {}
    for p in primes:
        a = Ppoly(coeffs(p, G, D), beta)
        Pchi = mp.fsum(a[m - 1]*mp.mpf(chi(p))**m for m in range(1, len(a) + 1))
        mn, x = minP(a)
        per[p] = (Pchi - mn, x); tot += Pchi - mn
    return tot, per

def solve(D, log=print):
    t0 = time.time()
    Cm, G, G0 = setup(D)
    nb = D + 1; nP = len(primes)
    cf = {p: np.array([[float(v) for v in row] for row in coeffs(p, G, D)]) for p in primes}   # (M_p, nb)
    obj = np.zeros(nb + nP)
    for ip, p in enumerate(primes):
        x = chi(p); obj[:nb] += sum(cf[p][m - 1]*x**m for m in range(1, cf[p].shape[0] + 1)); obj[nb + ip] = -1
    xs = {p: list(np.linspace(-1, 1, 41)) for p in primes}
    tg = np.concatenate([np.linspace(0, 60, 1201), 60*1.03**np.arange(1, 200)])
    Bt = np.array([[float(mp.fsum(Cm[i][k]*mp.mpf(t)**(2*k) for k in range(nb))) for i in range(nb)] for t in tg])
    Bt = Bt/np.abs(Bt).max(axis=1, keepdims=True)
    lead = np.array([float(Cm[i][D]) for i in range(nb)]); lead = lead/np.abs(lead).max()
    for rd in range(60):
        rows = []
        for ip, p in enumerate(primes):
            for x in xs[p]:
                r = np.zeros(nb + nP); r[:nb] = -sum(cf[p][m - 1]*x**m for m in range(1, cf[p].shape[0] + 1)); r[nb + ip] = 1
                rows.append(r)
        A = np.array(rows); b = np.zeros(len(rows))
        AB = np.hstack([-Bt, np.zeros((len(tg), nP))]); A = np.vstack([A, AB, np.concatenate([-lead, np.zeros(nP)])[None, :]])
        b = np.concatenate([b, np.zeros(len(tg) + 1)])
        Aeq = np.concatenate([[float(v) for v in G0], np.zeros(nP)])[None, :]
        res = linprog(obj, A_ub=A, b_ub=b, A_eq=Aeq, b_eq=[1.0], bounds=[(None, None)]*(nb + nP), method='highs',
                      options={'primal_feasibility_tolerance': 1e-10, 'dual_feasibility_tolerance': 1e-10})
        if res.status != 0:
            log('  LP status %d: %s' % (res.status, res.message)); break
        beta = res.x[:nb]; s = res.x[nb:]; added = 0
        fine = np.linspace(-1, 1, 4001)
        for ip, p in enumerate(primes):
            Pv = sum((cf[p][m - 1] @ beta)*fine**m for m in range(1, cf[p].shape[0] + 1))
            j = int(np.argmin(Pv))
            if Pv[j] < s[ip] - 1e-13*max(1, abs(s[ip])):
                xs[p].append(float(fine[j])); added += 1
        if not added: break
    bm = [mp.mpf(float(v)) for v in beta]
    dtot, per = defect(bm, G, D)
    g0 = mp.fsum(G0[i]*bm[i] for i in range(nb))
    bcoef = [mp.fsum(bm[i]*Cm[i][k] for i in range(nb)) for k in range(nb)]     # B(t) = sum_k bcoef[k] t^{2k}
    log('deg %d: LP %.3e (rounds %d), exact Defect/g(0) = %s, bound 5 exp(-2 D/g0) = %s  [%.0f s]' % (
        2*D, res.fun if res.status == 0 else float('nan'), rd + 1, mp.nstr(dtot/g0, 8), mp.nstr(5*mp.exp(-2*dtot/g0), 15), time.time() - t0))
    top = sorted(per.items(), key=lambda kv: -kv[1][0])[:6]
    log('   largest per-prime defects: ' + ', '.join('%d:%s@x=%s' % (p, mp.nstr(d/g0, 3), mp.nstr(x, 4)) for p, (d, x) in top))
    gl = ', '.join('%d:%s' % (n, mp.nstr(mp.fsum(G[n][i]*bm[i] for i in range(nb))/g0, 4)) for n in sorted(G)[:14])
    log('   g(log n)/g(0): ' + gl)
    return dict(D=D, beta=[mp.nstr(v, 30) for v in bm], b=[mp.nstr(v, 40) for v in bcoef], defect=mp.nstr(dtot/g0, 20),
                per={str(p): [mp.nstr(d/g0, 6), mp.nstr(x, 8)] for p, (d, x) in per.items() if d/g0 > mp.mpf(10)**-30})

if __name__ == '__main__':
    for D in map(int, sys.argv[2:]):
        r = solve(D)
        json.dump(r, open('/home/bas/praime/writeup/code/study/bootstrap/maass/theory/c5opt_d%d.json' % (2*D), 'w'), indent=1)
