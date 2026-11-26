# Sequential bootstrap, second version.  Input: the Gamma factor and pole of zeta, the Beurling form of
# the Euler product, positivity.  No primes are put in.
# Per support X: (1) float64: if positivity fails, add one point by a grid scan and polish all points
# jointly (Nelder-Mead on the smallest eigenvalue); (2) high precision: Gauss-Newton on the block of
# eigenvalues below 2|lam_min|, with a step guard, until the form is positive.
import sys, time, numpy as np
from scipy.optimize import minimize
from mpmath import mp, mpf, log, eigsy, matrix, lu_solve, sqrt, nstr
from weil import build_A, comb, sectors
import scan as boot2
DPS = int(sys.argv[1]); mp.dps = DPS
def lam_min(a, N, A, pos):
    E, O = sectors(a, N, A - comb(a, N, pos))
    return min(min(eigsy(E, eigvals_only=True)), min(eigsy(O, eigvals_only=True)))
def gn_step(a, N, A, pos):
    Bs = sectors(a, N, A - comb(a, N, pos))
    dec = [eigsy(B) for B in Bs]
    lmin = min(min(d[0]) for d in dec)
    if lmin >= 0: return None, lmin
    h = mpf(10)**(-DPS//3); L = abs(lmin)
    Ms = []
    for k in range(len(pos)):
        pp = list(pos); pp[k] = pos[k] + h; pm = list(pos); pm[k] = pos[k] - h
        Ms.append(sectors(a, N, (comb(a, N, pp) - comb(a, N, pm))/(2*h)))
    # block: eigenvalues below 1e4 |lam_min| (and at least 3 per sector).  Target: remove the negative
    # eigenvalues, leave the rest of the block unchanged; rows weighted by the scale of the eigenvalues.
    J = []; b = []
    for bi, (ev, V) in enumerate(dec):
        order = sorted(range(len(ev)), key=lambda i: ev[i])
        idx = [i for i in order if ev[i] < 10**4*L]
        if len(idx) < 3: idx = order[:3]
        cols = matrix(V.rows, len(idx))
        for c, i in enumerate(idx):
            for r in range(V.rows): cols[r, c] = V[r, i]
        red = [cols.T*Ms[k][bi]*cols for k in range(len(pos))]
        sc = [sqrt(max(abs(ev[i]), L)) for i in idx]
        for i in range(len(idx)):
            for j in range(i, len(idx)):
                wgt = (1 if i == j else sqrt(mpf(2)))/(sc[i]*sc[j])
                J.append([-wgt*red[k][i, j] for k in range(len(pos))])
                b.append(wgt*(-ev[idx[i]]) if (i == j and ev[idx[i]] < 0) else mpf(0))
    d = lu_solve(matrix(J), matrix(b))
    return [mp.re(v) for v in d], lmin
def refine(a, N, A, pos, iters=25):
    pos = [mpf(x) for x in pos]; lmin = lam_min(a, N, A, pos)
    for it in range(iters):
        if lmin >= 0: break
        d, _ = gn_step(a, N, A, pos)
        if d is None: break
        s = mpf(1); ok = False
        while s > mpf(1)/64:
            new = [pos[k] + s*d[k] for k in range(len(pos))]
            if min(new) <= 1 or max(new) > mp.e**(2*a): s /= 2; continue
            l2 = lam_min(a, N, A, new)
            if l2 > lmin: pos, lmin, ok = new, l2, True; break
            s /= 2
        if not ok: break
    return pos, lmin
def polish(F, pos):
    if not pos: return pos, F.mineig([])
    f = lambda v: 1e3 if (min(v) <= 1.001 or max(v) > F.X) else -F.mineig(list(v))
    r = minimize(f, np.array(pos), method='Nelder-Mead', options=dict(xatol=1e-13, fatol=1e-16, maxiter=4000, initial_simplex=None))
    r = minimize(f, r.x, method='Nelder-Mead', options=dict(xatol=1e-14, fatol=1e-17, maxiter=4000))
    return sorted(r.x), -r.fun
if __name__ == '__main__':
    pos = []; t0 = time.time()
    for Xs in sys.argv[2].split(','):
        X = mpf(Xs); a = log(X)/2
        N = max(12, int(8 + 4*log(X)/log(2)))
        A = build_A(a, N); F = boot2.Form(float(X), N)
        fp = [float(x) for x in pos]; added = []; dropped = []
        fp, e = polish(F, fp); fp = list(fp)
        for trial in range(3):
            if e > -1e-10: break
            grid = np.arange(1.03, float(X), 0.01)
            vals = [F.mineig(fp + [y]) for y in grid]
            y = grid[int(np.argmax(vals))]
            cand, e2 = polish(F, fp + [y])
            if e2 <= e + 1e-13: break
            added.append(round(float(y), 2)); fp, e = list(cand), e2
        for x in sorted(fp, reverse=True):          # prune: keep only points that positivity needs
            rest = [v for v in fp if v != x]
            if F.mineig(rest) > -1e-10:
                dropped.append(round(x, 3)); fp = rest
        if added or dropped: pos = [mpf(x) for x in fp]
        else: pos = [p if abs(float(p) - q) < 1e-9 else mpf(q) for p, q in zip(sorted(pos), sorted(fp))]
        pos, lmin = refine(a, N, A, sorted(pos))
        print('support X=%s N=%d  float lam_min %.2e  added near %s dropped %s  final lam_min=%s  (%.0fs)' % (Xs, N, e, added, dropped, nstr(lmin, 3), time.time() - t0))
        for x in pos:
            r = mp.nint(x); print('     %s   |x - nearest integer| = %s' % (nstr(x, 28), nstr(abs(x - r), 3)))
        sys.stdout.flush()
