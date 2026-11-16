# Item E (race): log-time share of x in [1000, 1e8] with #{chi(p)=+1} > #{chi(p)=-1}, p <= x, p != P.
import numpy as np, math
X = 10**8
# odd-only sieve: index i <-> 2i+1
h = X // 2 + 1
s = np.ones(h, bool); s[0] = False
for i in range(1, (int(X**0.5) - 1) // 2 + 1):
    if s[i]:
        p = 2 * i + 1
        s[p * p // 2::p] = False
primes = np.concatenate(([2], 2 * np.nonzero(s)[0] + 1)).astype(np.int64)
primes = primes[primes <= X]
del s
print('pi(1e8) =', len(primes))
lp = np.log(primes.astype(np.float64))
lo, hi = math.log(1000), math.log(X)
# interval k: [p_k, p_{k+1}) clipped to [1000, 1e8]
left = np.maximum(lp, lo)
right = np.minimum(np.append(lp[1:], hi), hi)
w = np.clip(right - left, 0, None)
assert abs(w.sum() - (hi - lo)) < 1e-9
for P in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]:
    qr = np.zeros(P, np.int8) - 1
    for a in range(1, P):
        qr[a * a % P] = 1
    qr[0] = 0
    chi = qr[primes % P].astype(np.int32)
    diff = np.cumsum(chi)
    share = w[diff > 0].sum() / (hi - lo)
    tie = w[diff == 0].sum() / (hi - lo)
    chi_odd = chi.copy(); chi_odd[0] = 0          # drop p = 2
    d2 = np.cumsum(chi_odd)
    share_odd = w[d2 > 0].sum() / (hi - lo)
    print('P = %2d: chi(2)=%+d  share incl. p=2: %.5f  share odd primes only: %.5f' % (P, chi[0], share, share_odd))
    print('P = %2d: share(+1 ahead) = %.4f   (ties %.4f)   first x>=1000 with +1 ahead: %s' % (
        P, share, tie, primes[np.nonzero((diff > 0) & (primes >= 1000))[0][0]] if np.any((diff > 0) & (primes >= 1000)) else '-'))
