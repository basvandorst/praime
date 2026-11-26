# Robust LP for certificates: F = sum_a c_a f_a (triangles, nested width grid), support [-A, A].
# Constraints: Fhat >= 0 on [1, 40] (grid), Fhat(0) = 0, moments int x^{2k} F = 0 for k in S, int_1^40 Fhat = 1.
# Minimise F(0) = sum c_a.  Then verify Fhat >= 0 on [1, 400] on a fine grid.
import numpy as np, sys
from scipy.optimize import linprog
def run(A, S, h=0.02, Xi=40.0, dxi=0.004):
    a = np.arange(h, A + 1e-9, h)
    xi = np.arange(1.0, Xi, dxi)
    M = a[None,:]*np.sinc(a[None,:]*xi[:,None])**2
    eq = [a] + [2*a**(2*k+1)/((2*k+1)*(2*k+2)) for k in S]
    norm = M.sum(axis=0)*dxi
    A_eq = np.vstack(eq + [norm]); b_eq = np.zeros(len(eq)+1); b_eq[-1] = 1.0
    r = linprog(np.ones(len(a)), A_ub=-M, b_ub=np.zeros(len(xi)), A_eq=A_eq, b_eq=b_eq,
                bounds=[(-50, 50)]*len(a), method='highs')
    if r.status != 0: return None
    c = r.x
    xv = np.arange(1.0, 400.0, 0.001)
    Fv = np.zeros_like(xv)
    for ai, ci in zip(a, c):
        if abs(ci) > 1e-12: Fv += ci*ai*np.sinc(ai*xv)**2
    mom = [float(np.dot(c, e)) for e in eq]
    return r.fun, Fv.min(), mom, np.abs(c).max()
for S in [[], [1], [2]]:
    for A in [1.0, 1.1, 1.2, 1.5, 2.0, 2.2]:
        out = run(A, S)
        if out is None: print("S=%s A=%.2f LP failed" % (S, A)); continue
        f0, mn, mom, cmax = out
        print("S=%-4s A=%.2f  min F(0) = %+.5f  | check: min Fhat on [1,400] = %+.2e, moments = %s, max|c| = %.1f"
              % (S, A, f0, mn, ["%.1e" % m for m in mom], cmax)); sys.stdout.flush()
