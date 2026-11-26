# How tightly is each prime pinned? For support e^L = X, move one prime p -> p +- 10^-k (others true),
# find the largest shift exponent k at which the truncated form still has a negative eigenvalue.
import sys, time
from mpmath import mp, mpf, log, eigsy
from sympy import primerange
from weil import build_A, comb, sectors
X = mpf(sys.argv[1]); N = int(sys.argv[2]); mp.dps = int(sys.argv[3])
a = log(X)/2
A = build_A(a, N)
P = [mpf(p) for p in primerange(2, int(X) + 1) if p < X]
def mineig(pos):
    E, O = sectors(a, N, A - comb(a, N, pos))
    return min(min(eigsy(E, eigvals_only=True)), min(eigsy(O, eigvals_only=True)))
E, O = sectors(a, N, A - comb(a, N, P))
ev = sorted(list(eigsy(E, eigvals_only=True)) + list(eigsy(O, eigvals_only=True)))
print('X', X, 'N', N, 'dps', mp.dps, 'lowest eigenvalues:', [mp.nstr(v, 3) for v in ev[:8]], flush=True)
kmax = mp.dps - 8
for i, p in enumerate(P):
    out = []
    for sgn in (+1, -1):
        lo, hi = 1, kmax          # negative at 10^-lo (checked), not negative at 10^-hi (checked)
        def neg(k):
            Q = list(P); Q[i] = p + sgn*mpf(10)**(-k)
            return mineig(Q) < 0
        if not neg(lo): out.append('feasible at 0.1'); continue
        if neg(hi): out.append('>%d' % hi); continue
        while hi - lo > 1:
            mid = (lo + hi)//2
            if neg(mid): lo = mid
            else: hi = mid
        out.append('1e-%d' % lo)
    print('  prime %d: infeasible down to shift +%s, -%s' % (int(p), out[0], out[1]), flush=True)
