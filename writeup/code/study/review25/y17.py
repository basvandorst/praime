"""Spot check of round 88: the k-th neighbour law per D at y = 13 and y = 17, own code (power series, own S_{>y},
nothing from the study) against r47/r88's construction at y = 17 for the same D."""
import sys, os, math, time, itertools
import numpy as np
LAM = math.log(1.05e9); K = 4
X = 10 ** 7
s = np.ones(X + 1, bool); s[:2] = False
for i in range(2, math.isqrt(X) + 1):
    if s[i]: s[i * i::i] = False
PR = np.nonzero(s)[0]

def law(D, y):
    SMALL = [q for q in (2, 3, 5, 7, 11, 13, 17) if q <= y]; PP = math.prod(SMALL)
    rho0 = 1 / (LAM * math.prod(1 - 1 / q for q in SMALL))
    MIDQ = [int(q) for q in PR if y < q <= D + 1]; FARQ = PR[PR > D + 1].astype(float)
    FAR = {k: float(np.sum(np.log1p(-k / FARQ) - k * np.log1p(-1 / FARQ))) for k in (2, 3, 4)}
    def lnS(sets):
        sets = np.asarray(sets, np.int64); N, k = sets.shape; out = np.full(N, FAR[k])
        for q in MIDQ:
            r = np.sort(sets % q, axis=1); nu = 1 + np.count_nonzero(np.diff(r, axis=1), axis=1)
            out += np.log(1 - nu / q) - k * math.log(1 - 1 / q)
        return out
    allc = np.array([a for a in range(PP) if math.gcd(a, PP) == 1], np.int64)
    cls = allc[np.gcd(allc + D, PP) == 1]
    pts = np.arange(2, D, 2, dtype=np.int64); n = len(pts)
    l2 = lnS([[0, D]])[0]
    l1 = lnS(np.column_stack([np.zeros(n, np.int64), np.full(n, D), pts])) - l2; p = rho0 * np.exp(l1)
    ij = np.array(list(itertools.combinations(range(n), 2)))
    R = np.zeros((n, n))
    ll = lnS(np.column_stack([np.zeros(len(ij), np.int64), np.full(len(ij), D), pts[ij[:, 0]], pts[ij[:, 1]]])) - l2
    rr = np.exp(ll - l1[ij[:, 0]] - l1[ij[:, 1]]) - 1; R[ij[:, 0], ij[:, 1]] = rr; R[ij[:, 1], ij[:, 0]] = rr
    M = (np.gcd(cls[:, None] + pts, PP) == 1).astype(np.float64)
    poly = np.zeros((len(cls), K)); poly[:, 0] = 1
    for i in range(n):
        m = M[:, i] * p[i]; new = poly * (1 - m)[:, None]; new[:, 1:] += poly[:, :-1] * m[:, None]; poly = new
    u = p / (1 - p); a = [-u] + [u * (-u) ** (t - 1) * (1 + u) for t in range(1, K)]
    Qs = np.zeros((len(cls), K))
    for m_ in range(K):
        AR = (M * a[m_]) @ R
        for t in range(m_, K):
            Qs[:, t] += 0.5 * np.sum(AR * (M * a[t - m_]), axis=1)
        del AR
    E = np.zeros_like(Qs); E[:, 0] = np.exp(Qs[:, 0])
    for t in range(1, K):
        E[:, t] = sum(mm * Qs[:, mm] * E[:, t - mm] for mm in range(1, t + 1)) / t
    G = np.zeros_like(Qs)
    for t in range(K):
        G[:, t] = sum(poly[:, mm] * E[:, t - mm] for mm in range(t + 1))
    return rho0 * math.exp(l2) * G.sum(0) / len(allc)

def r88_style(D):
    os.environ['R47_Y'] = '17'; sys.path.insert(0, '/home/bas/praime/writeup/code/study')
    import r47
    NZ = 48; Z = np.exp(2j * np.pi * np.arange(NZ) / NZ); NCLS = 92160
    RHO0 = 1 / (LAM * r47.RHO_DEN)
    st = r47.setup(D); M = st['M']; n = M.shape[1]; p = RHO0 * st['s1']
    logG = np.zeros((M.shape[0], NZ), complex)
    logG += M @ np.log1p(p[:, None] * (Z[None, :] - 1))
    W = p[:, None] * (Z[None, :] - 1) / (1 + p[:, None] * (Z[None, :] - 1))
    for m in range(NZ):
        Xm = M * W[:, m][None, :]
        logG[:, m] += 0.5 * np.sum((Xm.real @ st['R'] + 1j * (Xm.imag @ st['R'])) * Xm, axis=1)
    coef = np.fft.fft(np.exp(logG), axis=1).real / NZ
    return RHO0 * st['S2'] * coef[:, :K].sum(0) / NCLS

if __name__ == '__main__':
    for D in map(int, sys.argv[1].split(',')):
        t0 = time.time(); a13 = law(D, 13); a17 = law(D, 17); b17 = r88_style(D)
        print(D, 'y13', np.round(a13, 9), '\n   y17 own', np.round(a17, 9), '\n   y17 r88', np.round(b17, 9),
              '\n   own17/r88-1', np.round(a17 / np.where(b17 != 0, b17, 1) - 1, 8), ' y17/y13-1', np.round(a17 / np.where(a13 != 0, a13, 1) - 1, 6), round(time.time() - t0), 's', flush=True)
