# Minimal-conductor principle, float64.  Datum: Gamma factor prod_i Gamma_R(s + kappa_i), simple pole,
# conductor 1.  For a positive Beurling system P let lam(P) = smallest eigenvalue of the truncated form at
# support X.  P is feasible for conductor q iff lam(P) + log q >= 0, so the smallest conductor that admits
# P is exp(-lam(P)).  Greedy search for max_P lam(P): add the best point, polish all, repeat.
import sys, numpy as np
from mpmath import mp, mpf, log
import datum as boot6
from weil import build_A
class FormG(boot6.Form):
    def __init__(self, X, N, kappas):
        boot6.Form.__init__(self, X, N, True, kappas[0], 1.0)
        n = 2*N + 1
        for k in kappas[1:]:
            B = build_A(log(mpf(X))/2, N, pole=False, kappa=k, logq=0)
            self.A = self.A + np.array([[float(B[i, j]) for j in range(n)] for i in range(n)])
def greedy(F, start=(), maxpts=14, verbose=False):
    pos = list(start); one = lambda p: [1]*len(p)
    pos, e = boot6.polish(F, pos, one(pos))
    while len(pos) < maxpts:
        grid = np.arange(1.03, F.X, 0.01)
        vals = [F.mineig(pos + [y], one(pos) + [1]) for y in grid]
        y = grid[int(np.argmax(vals))]
        cand, e2 = boot6.polish(F, pos + [y], one(pos) + [1])
        if e2 <= e + 1e-9: break
        pos, e = cand, e2
        if verbose: print('     ', len(pos), 'points, lam', e, flush=True)
    return sorted(pos), e
if __name__ == '__main__':
    kappas = [int(c) for c in sys.argv[1]]
    pos = []
    for X in [float(v) for v in sys.argv[2].split(',')]:
        N = max(12, int(8 + 4*np.log(X)/np.log(2)))
        F = FormG(X, N, kappas)
        pos, e = greedy(F, start=pos)                 # continuation: start from the previous support
        pos = [x for x in pos if x < X - 0.02]        # points at the edge carry no weight
        e = F.mineig(pos, [1]*len(pos))
        print('kappas %s X=%-4g max lam found %.8f  -> smallest conductor found %.6f   points: %s' % (sys.argv[1], X, e, np.exp(-e), ' '.join('%.5f' % x for x in pos)), flush=True)
