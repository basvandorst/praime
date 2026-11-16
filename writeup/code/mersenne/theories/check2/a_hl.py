# Item A: Hardy-Littlewood predictions for wheel-24 bytes, chi-square, entropies, MI.
import sys, numpy as np, math
from sympy import primerange

OFFS = [-1, 1, 5, 7, 11, 13, 17, 19]
R = 10**7
primes_big = np.array(list(primerange(47, R)), dtype=np.float64)

def singular_series(offs):
    """Return S(U) for all subsets U of offs, mask bit (n-1-i) <-> position i."""
    n = len(offs)
    N = 1 << n
    pc = np.array([bin(m).count('1') for m in range(N)])
    logS = np.zeros(N)
    # r = 2, 3: nu = 0 -> factor 1/(1-1/r)^k
    for r in (2, 3):
        logS += -pc * math.log(1 - 1 / r)
    # r = 5..43: compute nu via residue bitmasks
    for r in primerange(5, 47):
        resmask = np.zeros(N, np.int64)
        for i, o in enumerate(offs):
            b = n - 1 - i
            v = 1 << (o % r)
            lo = 1 << b
            # subsets containing bit b: those m with bit b set
            idx = np.nonzero(np.arange(N) & lo)[0]
            resmask[idx] |= v
        nu = np.array([bin(int(x)).count('1') for x in resmask])
        with np.errstate(divide='ignore'):
            logS += np.log(np.maximum(1 - nu / r, 0.0)) - pc * math.log(1 - 1 / r)
    # r >= 47: nu = k (all offsets distinct mod r since max diff 44 < 47)
    fk = np.zeros(n + 1)
    tail_sum = 1 / (R * math.log(R))   # sum_{r>R} 1/r^2 approx
    for k in range(n + 1):
        fk[k] = np.sum(np.log1p(-k / primes_big) - k * np.log1p(-1 / primes_big)) - k * (k - 1) / 2 * tail_sum
    logS += fk[pc]
    return np.exp(logS), pc

def lam_sums(J, kmax, shift):
    lam = np.zeros(kmax + 1)
    for a in range(1, J + 1, 5_000_000):
        b = min(J + 1, a + 5_000_000)
        L = 1 / np.log(24.0 * np.arange(a, b) + shift)
        p = np.ones_like(L)
        for k in range(kmax + 1):
            lam[k] += p.sum()
            p *= L
    return lam

def exact_from_all(f, n):
    """g(T) = sum_{U superset T} (-1)^{|U-T|} f(U)."""
    g = f.copy()
    N = 1 << n
    ar = np.arange(N)
    for b in range(n):
        lo = (ar & (1 << b)) == 0
        g[lo] -= g[ar[lo] | (1 << b)]
    return g

def H(p):
    p = p[p > 0]
    return -np.sum(p * np.log2(p))

S8, pc8 = singular_series(OFFS)
print('S({-1,+1}) / (12*C2) =', S8[0b11000000] / (12 * 0.6601618158468696))
imp = np.nonzero(S8 == 0)[0]
print('impossible byte values:', len(imp))

for J in (416666, 41666666):
    obs = np.load('A_counts_%d.npy' % J).astype(float)
    lam = lam_sums(J, 8, 9.0)
    E = exact_from_all(S8 * lam[pc8], 8)
    print('\nJ =', J, ' sum E =', E.sum(), ' min E (possible) =', E[S8 > 0].min(), ' max |E| impossible =', np.abs(E[S8 == 0]).max())
    print('  observed in impossible values:', obs[S8 == 0].sum())
    poss = S8 > 0
    chi_all = np.sum((obs[poss] - E[poss]) ** 2 / E[poss])
    ge5 = E >= 5
    chi_5 = np.sum((obs[ge5] - E[ge5]) ** 2 / E[ge5])
    print('  possible values: %d  chi2 = %.2f' % (poss.sum(), chi_all))
    print('  values with E>=5: %d  chi2 = %.2f' % (ge5.sum(), chi_5))
    print('  entropy observed %.5f  predicted %.5f bits' % (H(obs / obs.sum()), H(E / E.sum())))

# mutual information of consecutive bytes, J = 41,666,666
J = 41666666
pair = np.load('A_pairs_%d.npy' % J).astype(float)
P = pair.reshape(256, 256) / pair.sum()
px, py = P.sum(1), P.sum(0)
mi_obs = H(px) + H(py) - H(P.ravel())
print('\nMI observed (plug-in) = %.5f bits  (pairs %d)' % (mi_obs, pair.sum()))
OFF16 = OFFS + [24 + o for o in OFFS]
S16, pc16 = singular_series(OFF16)
for shift in (9.0, 21.0):
    lam16 = lam_sums(J - 1, 16, shift)
    E16 = exact_from_all(S16 * lam16[pc16], 16)
    neg = E16[E16 < 0]
    Q = np.clip(E16, 0, None).reshape(256, 256); Q /= Q.sum()
    mi_pred = H(Q.sum(1)) + H(Q.sum(0)) - H(Q.ravel())
    print('MI predicted (16-pos HL, L = 1/ln(24j+%g)) = %.5f bits; negative cells %d (min %.3g)' % (shift, mi_pred, len(neg), neg.min() if len(neg) else 0))
    # consistency: marginal of the predicted joint vs 8-position prediction
# MI of predicted joint at a single j (no mixing over j), at j = J/2
Lmid = 1 / math.log(24 * J / 2 + 21)
E1 = exact_from_all(S16 * Lmid ** pc16.astype(float), 16)
Q = np.clip(E1, 0, None).reshape(256, 256); Q /= Q.sum()
print('MI predicted at single j=J/2 (no mixture) = %.5f' % (H(Q.sum(1)) + H(Q.sum(0)) - H(Q.ravel())))
