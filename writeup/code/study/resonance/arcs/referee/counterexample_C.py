# Stucky Lemma 4.2, case T <= M ("Ricci Thm C"): claims, for ARBITRARY coefficients c(a) on ideals with
# N a in [N, 2N]:  sum_{|m|<=M} int_{-T}^{T} |sum_a c(a) lambda^m(a) Na^{-it}|^2 dt << (MT + N min(M,T)) sum|c|^2.
# Counterexample: c = 1 on the ideals (k), k rational integer (arg = 0, so lambda^m((k)) = 1 for all m).
# Exact LHS = (2M+1) * int_{-T}^{T} |sum_k k^{-2it}|^2 dt = (2M+1) * sum_{k,l} 2 sin(2T log(k/l)) / (2 log(k/l)) (k != l), 2T (k = l).
import numpy as np
from mpmath import mp
for N in [10**4, 10**6, 10**8]:
    ks = np.arange(int(np.ceil(np.sqrt(N))), int(np.floor(np.sqrt(2*N))) + 1, dtype=np.float64)
    K = len(ks); T = 1.0
    lk = np.log(ks)
    # int_{-T}^{T} e^{-2it(log k - log l)} dt = 2 sin(2T d)/(2d), d = log k - log l
    tot = 0.0
    for i in range(K):
        d = lk[i] - lk
        v = np.where(d == 0, 2*T, 2*np.sin(2*T*d)/np.where(d == 0, 1, 2*d))
        tot += v.sum()
    for M in [N, int(N**0.75)]:
        lhs = (2*M + 1)*tot
        sumc2 = K
        D_claim = M*T + N*min(M, T)
        D_elem = N + M*T + (M + T)*np.sqrt(N)
        print(f"N={N:.0e} M={M:.2e} T=1 K={K}: LHS/sum|c|^2 = {lhs/sumc2:.3e}; claimed D = MT+N min(M,T) = {D_claim:.3e}; "
              f"ratio = {lhs/sumc2/D_claim:.2f}; elementary D = N+MT+(M+T)sqrt N = {D_elem:.3e}, ratio = {lhs/sumc2/D_elem:.2f}")
