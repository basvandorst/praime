"""High-precision Weil form D_P on [-L/2, L/2] in the monomial basis (2x/L)^i.
Exact polynomial correlations (python-flint fmpq_poly), mpmath for the rest.
Same conventions as wform.py."""
import mpmath as mp
import flint, pickle, os

CACHE = '/tmp/claude-1000/bootstrap_theory/cij_%d.pkl'


def corr_polys(N):
    """S[i][j] = coefficient list (as Fractions strings) of c_ij(v)+c_ji(v), and c0[i][j] = c_ij(0).
    c_ij(v) = int_{v-1}^{1} s^i (s-v)^j ds, 0<=v<=2."""
    fn = CACHE % N
    if os.path.exists(fn):
        return pickle.load(open(fn, 'rb'))
    v = flint.fmpq_poly([0, 1])
    one_minus_v = flint.fmpq_poly([1, -1])
    pw = [flint.fmpq_poly([1])]
    for n in range(2 * N + 2):
        pw.append(pw[-1] * one_minus_v)
    vp = [flint.fmpq_poly([1])]
    for n in range(N + 1):
        vp.append(vp[-1] * v)
    from math import comb
    c = [[None] * N for _ in range(N)]
    for i in range(N):
        for j in range(N):
            acc = flint.fmpq_poly([0])
            for k in range(i + 1):
                n = j + k + 1
                term = (pw[n] - flint.fmpq_poly([(-1) ** n])) * vp[i - k] * flint.fmpq(comb(i, k), n)
                acc += term
            c[i][j] = acc
    S = [[[str(x) for x in (c[i][j] + c[j][i]).coeffs()] for j in range(N)] for i in range(N)]
    c0 = [[str(c[i][j].coeffs()[0]) if c[i][j].coeffs() else '0' for j in range(N)] for i in range(N)]
    pickle.dump((S, c0), open(fn, 'wb'))
    return S, c0


def frac(s):
    if '/' in s:
        a, b = s.split('/')
        return mp.mpf(int(a)) / int(b)
    return mp.mpf(int(s))


def kfun(u):
    return mp.exp(-u / 2) / (1 - mp.exp(-2 * u))


def form(N, L, P, dps=80):
    mp.mp.dps = dps
    L = mp.mpf(L)
    S, c0 = corr_polys(N)
    Sm = [[[frac(x) for x in S[i][j]] for j in range(N)] for i in range(N)]
    C0 = [[frac(c0[i][j]) for j in range(N)] for i in range(N)]
    a0 = mp.digamma(mp.mpf(1) / 4) - mp.log(mp.pi)
    K = [None] + [mp.quad(lambda u: kfun(u) * u ** m, [0, L / 2, L]) for m in range(1, 2 * N + 2)]
    ktail = mp.quad(kfun, [L, mp.inf])
    al = [L / 2 * mp.quad(lambda s: s ** i * mp.exp(L * s / 4), [-1, 1]) for i in range(N)]
    be = [(-1) ** i * al[i] for i in range(N)]
    atoms = []
    for x in P:
        lx = mp.log(x)
        k = 1
        while k * lx < L:
            atoms.append((2 * k * lx / L, lx * mp.mpf(x) ** (-mp.mpf(k) / 2)))
            k += 1
    Q = mp.matrix(N, N)
    G = mp.matrix(N, N)
    r = 2 / L
    for i in range(N):
        for j in range(i, N):
            s = Sm[i][j]
            arch_int = mp.fsum(s[m] * r ** m * K[m] for m in range(1, len(s)))
            val = al[i] * be[j] + be[i] * al[j]
            val += L / 2 * (a0 * C0[i][j] - arch_int + 2 * C0[i][j] * ktail)
            if atoms:
                val -= L / 2 * mp.fsum(w * mp.polyval(s[::-1], vv) for vv, w in atoms)
            Q[i, j] = Q[j, i] = val
            G[i, j] = G[j, i] = L / 2 * C0[i][j]
    return Q, G


def eig_gen(Q, G, k=3):
    R = mp.cholesky(G)          # G = R R^T, R lower
    Ri = mp.inverse(R)
    M = Ri * Q * Ri.T
    M = (M + M.T) / 2
    E = mp.eigsy(M, eigvals_only=True)
    E = sorted([E[i] for i in range(len(E))])
    return E[:k]
