# Post-registration tools for round 165: t*(r) with a refined angle grid, the peak of t*(r) and the zero crossings.
# Two forms (L(f), Sym^2) by default; Sym^3 and Sym^4 are slack near the edges (checked in the scans).
# Usage: refine.py peak <eps> <X> <r_lo> <r_hi>   |   refine.py edge <eps> <X> <r_infeasible> <r_feasible>   |   refine.py at <eps> <X> <r>
import sys, numpy as np
from maass import Problem

HOL = False
def tstar(r, eps, X, levels=2, kmax=2, tol=1e-8, N=None, info=False):
    th = np.linspace(0, np.pi, 61); d = np.pi/60
    for lev in range(levels + 1):
        P = Problem(r, eps, X, kmax=kmax, th=th, N=N, hol=HOL); lam, t, rd = P.solve(tol=tol, ncut=6)
        if lev == levels: break
        atoms = sorted({float(P.th[j]) for w in P.w for j in np.nonzero(w > 1e-7)[0]})
        d /= 4
        new = {round(min(np.pi, max(0.0, a + k*d)), 12) for a in atoms for k in range(-4, 5)}
        th = np.array(sorted(set(np.round(np.linspace(0, np.pi, 61), 12)) | new))
    return (t, P) if info else t

def peak(eps, X, lo, hi, iters=34):
    g = (np.sqrt(5) - 1)/2; a, b = lo, hi
    c = b - g*(b - a); d = a + g*(b - a); fc = tstar(c, eps, X); fd = tstar(d, eps, X)
    for _ in range(iters):
        if fc > fd: b, d, fd = d, c, fc; c = b - g*(b - a); fc = tstar(c, eps, X)
        else: a, c, fc = c, d, fd; d = a + g*(b - a); fd = tstar(d, eps, X)
    return (a + b)/2, max(fc, fd)

def edge(eps, X, bad, good, iters=30):
    for _ in range(iters):
        mid = (bad + good)/2
        if tstar(mid, eps, X) >= 0: good = mid
        else: bad = mid
    return bad, good

if __name__ == '__main__':
    mode = sys.argv[1]; eps = int(sys.argv[2]); X = float(sys.argv[3])
    if eps == 2: HOL = True; eps = 0          # eps = 2: holomorphic weight as the parameter
    if mode == 'peak':
        r, t = peak(eps, X, float(sys.argv[4]), float(sys.argv[5]))
        t2, P = tstar(r, eps, X, info=True)
        print('peak eps %d X %g: r = %.8f, t* = %.3e; a_p: %s' % (eps, X, r, t2, P.summary()), flush=True)
    elif mode == 'edge':
        bad, good = edge(eps, X, float(sys.argv[4]), float(sys.argv[5]))
        t2, P = tstar(good, eps, X, info=True)
        print('edge eps %d X %g: t* changes sign between %.8f and %.8f; a_p at the feasible side: %s' % (eps, X, bad, good, P.summary()), flush=True)
    else:
        for r in sys.argv[4].split(','):
            t2, P = tstar(float(r), eps, X, info=True)
            print('at eps %d X %g r %s: t* = %.4e; per form %s; a_p: %s' % (eps, X, r, t2, ' '.join('%.2e' % v for v in P.mineigs(P.w)), P.summary()), flush=True)
