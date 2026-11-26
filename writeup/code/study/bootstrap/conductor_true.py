# The smallest conductor that admits (a) the true prime-ideal norms of a field, (b) the local maximiser
# reached from them after a perturbation, at support X.  Float64.
import numpy as np
import datum as boot6
from conductor import FormG
rng = np.random.default_rng(1)
cases = [('Q(sqrt-3)', [0, 1], 3, [3, 4, 7, 7, 13, 13]),
         ('Q(i)     ', [0, 1], 4, [2, 5, 5, 9, 13, 13]),
         ('Q(sqrt 5)', [0, 0], 5, [4, 5, 9, 11, 11]),
         ('Q(sqrt 2)', [0, 0], 8, [2, 7, 7, 9])]
for name, kappas, d, norms in cases:
    for X in (4.0, 6.0, 8.0, 10.0, 12.0):
        N = max(12, int(8 + 4*np.log(X)/np.log(2)))
        F = FormG(X, N, kappas)
        p = [float(x) for x in norms if x < X]
        lt = F.mineig(p, [1]*len(p))
        start = [x*(1 + 0.01*rng.standard_normal()) for x in p]
        q, lq = boot6.polish(F, start, [1]*len(p))
        print('%s d=%d X=%-4g true norms: conductor %.6f | perturbed 1%% then maximised: conductor %.6f at %s' % (name, d, X, np.exp(-lt), np.exp(-lq), ' '.join('%.5f' % v for v in sorted(q))), flush=True)
