"""Independent twin k-th neighbour law at a few d (nothing imported from the study). Own S_{>13}: exact residue
counts for 13 < q <= 1300 (spans <= 1202), generic factor for 1300 < q < 1e7 (all points distinct mod q).
P(k-th next twin start at d) = rho0^2 S(B)/S({0,2}) E_a[ [z^{k-1}] prod_h (1 - p_h + p_h z)^{M_ah} exp(Q_a(z)) ],
B = {0,2,d,d+2}, p_h = rho0^2 S(B u {h,h+2})/S(B), Q = 1/2 sum_{h!=k} R_hk w_h w_k, as a power series in z."""
import sys, math, itertools
import numpy as np
X = 10 ** 7
s = np.ones(X + 1, bool); s[:2] = False
for i in range(2, math.isqrt(X) + 1):
    if s[i]: s[i * i::i] = False
PR = np.nonzero(s)[0]
MIDQ = [int(q) for q in PR if 13 < q <= 1300]
FARQ = PR[PR > 1300].astype(float)
FAR = {k: float(np.sum(np.log1p(-k / FARQ) - k * np.log1p(-1 / FARQ))) for k in (2, 4, 6, 8)}
PP = 30030; DEN = math.prod(1 - 1 / q for q in (2, 3, 5, 7, 11, 13))

def lnS(sets):
    sets = np.asarray(sets, np.int64); N, k = sets.shape
    out = np.full(N, FAR[k])
    for q in MIDQ:
        r = np.sort(sets % q, axis=1); nu = 1 + np.count_nonzero(np.diff(r, axis=1), axis=1)
        out += np.log(1 - nu / q) - k * math.log(1 - 1 / q)
    return out

def law(d, lam, K=4):
    rho0 = 1 / (lam * DEN)
    allc = np.array([a for a in range(PP) if math.gcd(a, PP) == 1 and math.gcd(a + 2, PP) == 1], np.int64)
    cls = allc[(np.gcd(allc + d, PP) == 1) & (np.gcd(allc + d + 2, PP) == 1)]
    h = np.arange(6, d, 6, dtype=np.int64); n = len(h)
    B = np.array([0, 2, d, d + 2])
    lB = lnS([B])[0]; l02 = lnS([[0, 2]])[0]
    base = rho0 ** 2 * math.exp(lB - l02)
    M = ((np.gcd(cls[:, None] + h, PP) == 1) & (np.gcd(cls[:, None] + h + 2, PP) == 1)).astype(float)
    l1 = lnS(np.column_stack([np.tile(B, (n, 1)), h, h + 2])) - lB
    p = rho0 ** 2 * np.exp(l1)
    R = np.zeros((n, n))
    ij = np.array(list(itertools.combinations(range(n), 2)))
    if len(ij):
        hi, hj = h[ij[:, 0]], h[ij[:, 1]]
        l2 = lnS(np.column_stack([np.tile(B, (len(ij), 1)), hi, hi + 2, hj, hj + 2])) - lB
        rr = np.exp(l2 - l1[ij[:, 0]] - l1[ij[:, 1]]) - 1
        R[ij[:, 0], ij[:, 1]] = rr; R[ij[:, 1], ij[:, 0]] = rr
    poly = np.zeros((len(cls), K)); poly[:, 0] = 1
    for i in range(n):
        m = M[:, i] * p[i]; new = poly * (1 - m)[:, None]; new[:, 1:] += poly[:, :-1] * m[:, None]; poly = new
    u = p / (1 - p)
    a = [-u] + [u * (-u) ** (t - 1) * (1 + u) for t in range(1, K)]
    A = [M * a[t] for t in range(K)]; AR = [A[t] @ R for t in range(K)]
    Qs = np.zeros((len(cls), K))
    for t in range(K):
        for m_ in range(t + 1):
            Qs[:, t] += 0.5 * np.sum(AR[m_] * A[t - m_], axis=1)
    E = np.zeros_like(Qs); E[:, 0] = np.exp(Qs[:, 0])
    for t in range(1, K):
        E[:, t] = sum(mm * Qs[:, mm] * E[:, t - mm] for mm in range(1, t + 1)) / t
    G = np.array([[sum(poly[c, m_] * E[c, t - m_] for m_ in range(t + 1)) for t in range(K)] for c in range(len(cls))])
    return base * G.sum(0) / len(allc)

if __name__ == '__main__':
    lam = math.log(1.5e9)
    ref = np.load('rlaw.npz'); dl = list(ref['d'])
    for d in map(int, sys.argv[1].split(',')):
        mine = law(d, lam)
        theirs = ref['out'][:, 1, dl.index(d)]
        print(d, ' '.join('k%d %.8e/%.8e rel %.1e' % (k + 1, mine[k], theirs[k], mine[k] / theirs[k] - 1 if theirs[k] else 0) for k in range(4)), flush=True)
