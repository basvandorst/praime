# Truncated Weil quadratic form on functions supported in [-a, a] (u = log x), a = L/2.
# Basis phi_k(u) = exp(i pi k u / a), k = -N..N; the matrix is real symmetric.
# Q = A - C(P): A = pole + archimedean part (fixed), C(P) = comb of a set P of generalized primes.
from mpmath import mp, mpf, matrix, exp, log, sqrt, pi, euler, sin, cos, cosh, eigsy
from mpmath.calculus.quadrature import GaussLegendre

def gl_nodes(lo, hi, panels, degree=6):
    base = GaussLegendre(mp).calc_nodes(degree, mp.prec)
    out = []
    h = (hi - lo) / panels
    for i in range(panels):
        c = lo + (i + mpf(1)/2) * h
        for x, w in base:
            out.append((c + x*h/2, w*h/2))
    return out

def S_matrix(a, N, w):
    """S_jk(w) = g_jk(w) + g_jk(-w) for 0 <= w <= 2a, j,k = -N..N."""
    n = 2*N + 1
    sn = [sin(pi*m*w/a) for m in range(-N, N+1)]
    cs = [cos(pi*m*w/a) for m in range(-N, N+1)]
    S = matrix(n, n)
    for j in range(n):
        S[j, j] = 2*(2*a - w)*cs[j]
        for k in range(j):
            d = j - k
            v = (-1)**d * (2*a/(pi*d)) * (sn[k] - sn[j])
            S[j, k] = v; S[k, j] = v
    return S

def build_A(a, N, pole=True, kappa=0, logq=0):
    n = 2*N + 1
    nodes = gl_nodes(mpf(0), 2*a, max(8, N//2 + 2))
    Ks = [mpf(0)]*(N+1); Ps = [mpf(0)]*(N+1); Dp = [mpf(0)]*(N+1); Da = [mpf(0)]*(N+1)
    for u, w in nodes:
        e2 = exp(-2*u); ek = exp(-(mpf(1)/2 + kappa)*u); k1 = ek/(1 - e2); ch = 2*cosh(u/2) if pole else mpf(0)
        for m in range(N+1):
            c = pi*m/a
            sm = sin(c*u); cm = cos(c*u)
            Ks[m] += w*k1*sm
            Ps[m] += w*ch*sm
            Dp[m] += w*ch*2*(2*a - u)*cm
            Da[m] += w*(4*a*e2 - 2*(2*a - u)*ek*cm)/(1 - e2)
    const = (logq - euler - log(pi))*2*a - 2*a*log(1 - exp(-4*a))
    A = matrix(n, n)
    sg = lambda m: (1 if m >= 0 else -1)
    for j in range(-N, N+1):
        A[j+N, j+N] = Dp[abs(j)] + Da[abs(j)] + const
        for k in range(-N, j):
            d = j - k
            pk = sg(k)*Ps[abs(k)]; pj = sg(j)*Ps[abs(j)]
            kk = sg(k)*Ks[abs(k)]; kj = sg(j)*Ks[abs(j)]
            v = (-1)**d * (2*a/(pi*d)) * ((pk - pj) - (kk - kj))
            A[j+N, k+N] = v; A[k+N, j+N] = v
    return A

def comb(a, N, positions, signs=None):
    n = 2*N + 1
    C = matrix(n, n)
    for i, x in enumerate(positions):
        x = mpf(x); lx = log(x); k = 1
        sg = 1 if signs is None else signs[i]
        while k*lx <= 2*a:
            C += sg**k * lx * x**(-mpf(k)/2) * S_matrix(a, N, k*lx)
            k += 1
    return C

def sectors(a, N, H):
    """Return (even, odd) matrices in orthonormal cos / sin bases."""
    E = matrix(N+1, N+1); O = matrix(N, N)
    nrm = lambda k: sqrt(2*a) if k == 0 else sqrt(4*a)
    for k in range(N+1):
        for l in range(N+1):
            ks = [k] if k == 0 else [k, -k]; ls = [l] if l == 0 else [l, -l]
            E[k, l] = sum(H[x+N, y+N] for x in ks for y in ls)/(nrm(k)*nrm(l))
    for k in range(1, N+1):
        for l in range(1, N+1):
            O[k-1, l-1] = (H[k+N, l+N] - H[k+N, -l+N] - H[-k+N, l+N] + H[-k+N, -l+N])/(4*a)
    return E, O

if __name__ == '__main__':
    import sys, time
    from sympy import primerange
    mp.dps = int(sys.argv[3]) if len(sys.argv) > 3 else 60
    X = mpf(sys.argv[1]); N = int(sys.argv[2])   # support: e^L = X
    a = log(X)/2
    t = time.time()
    A = build_A(a, N)
    P = list(primerange(2, int(X) + 1))
    H = A - comb(a, N, P)
    E, O = sectors(a, N, H)
    ev = sorted(eigsy(E, eigvals_only=True)); od = sorted(eigsy(O, eigvals_only=True))
    print('X', X, 'N', N, 'time %.0fs' % (time.time() - t))
    print(' even:', [mp.nstr(v, 5) for v in ev[:6]])
    print(' odd :', [mp.nstr(v, 5) for v in od[:6]])
    E0, O0 = sectors(a, N, A)
    print(' no primes, even:', [mp.nstr(v, 5) for v in sorted(eigsy(E0, eigvals_only=True))[:4]])
