"""Referee check S1: trace formula of T_p o W_N on the new space, N squarefree composite.
LHS from PARI eigenforms (mfeigenbasis, mfatkineigenvalues), RHS from Hurwitz class numbers (qfbhclassno).
Also: sign convention of W_N vs root number by the functional equation (lfuncheckfeq with both signs),
old forms: trace of T_p W_N on the full cusp space (matrices) vs on the new space."""
import cypari2
from fractions import Fraction as Fr
pari = cypari2.Pari()
pari.allocatemem(400 * 10 ** 6)
gp = pari

def H(n):
    return Fr(str(gp(f'qfbhclassno({n})')))

def U(m, x2):
    # U_m(x) for even m as polynomial in x^2 (x2 = x^2), exact: recurrence U0=1, U1=2x, U_{n+1}=2x U_n - U_{n-1}
    # work with x symbolic via pairs (even part coefficient list in x)
    # simple: build polynomial coefficients
    c = [[1], [0, 2]]
    for n in range(1, m):
        a = [0] + [2 * t for t in c[n]]
        b = c[n - 1] + [0] * (len(a) - len(c[n - 1]))
        c.append([a[i] - b[i] for i in range(len(a))])
    poly = c[m]
    s = Fr(0)
    for i, t in enumerate(poly):
        if t:
            assert i % 2 == 0
            s += t * x2 ** (i // 2)
    return s

def rhs(N, k, p):
    s = Fr(1, 2) * H(4 * p * N)
    sg = (-1) ** (k // 2 - 1)
    r = 1
    while r * r * N * N < 4 * p * N:
        s += sg * U(k - 2, Fr(r * r * N, 4 * p)) * H(4 * p * N - r * r * N * N)
        r += 1
    if k == 2:
        s -= p + 1
    return s

def omega(N):
    return len(gp(f'factor({N})[,1]'))

def lhs_data(N, k):
    gp(f'mf = mfinit([{N},{k}],0)')
    d = int(gp('mfdim(mf)'))
    if d == 0:
        return d, None
    gp('L = mfeigenbasis(mf); W = mfatkineigenvalues(mf, %d); F = mffields(mf)' % N)
    norb = int(gp('#L'))
    orbs = []
    for i in range(1, norb + 1):
        wv = [int(gp(f'W[{i}][{j}]')) for j in range(1, int(gp(f'#W[{i}]')) + 1)]
        assert len(set(wv)) == 1, (N, k, wv)
        deg = int(gp(f'poldegree(F[{i}])'))
        orbs.append((i, wv[0], deg))
    return d, orbs

def lhs(N, k, p, orbs):
    if orbs is None:
        return Fr(0)
    tot = Fr(0)
    for i, w, deg in orbs:
        a = gp(f'my(c=mfcoef(L[{i}],{p})); if(type(c)=="t_POLMOD", trace(c), {deg}*c)')
        tot += w * Fr(str(a))
    # eps_f = (-1)^{k/2} w_N(f); sum eps_f a_f(p) / p^{(k-2)/2}
    return (-1) ** (k // 2) * tot / Fr(p) ** ((k - 2) // 2)

# ---- 1. sign convention: root number from the functional equation for rational newforms
print('== root number from functional equation (lfuncheckfeq, both signs) vs (-1)^(k/2) w_N from mfatkineigenvalues')
nfe = 0
for (N, k) in [(14, 2), (15, 2), (21, 2), (30, 2), (42, 2), (11, 2), (37, 2), (6, 4), (10, 4), (14, 4), (15, 4), (21, 4), (30, 4)]:
    d, orbs = lhs_data(N, k)
    if orbs is None:
        continue
    for i, w, deg in orbs:
        if deg != 1:
            continue
        gp(f'co = mfcoefs(L[{i}], 3000)[^1]')
        res = []
        for eps in (1, -1):
            v = float(gp(f'lfuncheckfeq(lfuncreate([co, 0, [0,1], {k}, {N}, {eps}]))'))
            res.append(v)
        eps_fe = 1 if res[0] < res[1] else -1
        eps_al = (-1) ** (k // 2) * w
        nfe += 1
        print(f'  N={N} k={k} orb {i}: feq bits (+1,-1) = {res[0]:.0f},{res[1]:.0f} -> eps={eps_fe}; (-1)^(k/2) w_N = {eps_al}',
              'OK' if eps_fe == eps_al else 'MISMATCH')

# ---- 2. trace formula on the new space
print('== trace formula, new space')
levels = [6, 10, 14, 15, 21, 30, 42, 66, 70, 78, 105, 110, 22, 33, 26, 55, 65, 77]
P = [int(x) for x in gp('primes([2,200])')]
bad = 0
cnt = 0
for k in (2, 4):
    for N in levels:
        if k == 4 and N > 110:
            continue
        d, orbs = lhs_data(N, k)
        nb = 0
        for p in P:
            if N % p == 0:
                continue
            a, b = lhs(N, k, p, orbs), rhs(N, k, p)
            cnt += 1
            if a != b:
                nb += 1
                if nb <= 3:
                    print('   MISMATCH', N, k, p, a, b)
        bad += nb
        print(f'  k={k} N={N:4d} dim_new={d:3d} omega={omega(N)}  mismatches {nb}')
print(f'trace formula: {cnt} pairs, {bad} mismatches')

# ---- 3. old forms: Tr(T_p W_N) on full cusp space vs new space (matrices)
print('== old forms: trace of T_p W_N, full cusp space vs new space')
for (N, k) in [(22, 2), (33, 2), (55, 2), (66, 2), (77, 2), (110, 2), (30, 4), (42, 4), (26, 4), (22, 4)]:
    out = []
    for flag in (1, 0):
        gp(f'mfx = mfinit([{N},{k}],{flag}); A = mfatkininit(mfx,{N}); M = A[2]')
        dim = int(gp('mfdim(mfx)'))
        if dim == 0:
            out.append(None); continue
        # normalise: M^2 must be a scalar c*Id; use M/sqrt(c) with c a rational square
        c = gp('(M^2)[1,1]')
        okscal = bool(gp('M^2 == (M^2)[1,1]*matid(#M)'))
        out.append((dim, okscal, str(c)))
    tr = []
    for p in [3, 5, 7, 13, 17, 19, 23, 29, 31, 37, 41, 43]:
        if N % p == 0:
            continue
        vals = []
        for flag in (1, 0):
            gp(f'mfx = mfinit([{N},{k}],{flag})')
            if int(gp('mfdim(mfx)')) == 0:
                vals.append(Fr(0)); continue
            gp(f'A = mfatkininit(mfx,{N}); M = A[2]; c=(M^2)[1,1]; s=sqrtint(numerator(c))/sqrtint(denominator(c))')
            assert bool(gp('issquare(c)')), 'W^2 scalar not a square'
            v = gp(f'trace(M*mfheckemat(mfx,{p}))/s')
            vals.append(Fr(str(v)))
        tr.append((p, vals[0], vals[1]))
    eq = all(a == b for _, a, b in tr)
    print(f'  N={N} k={k} spaces(full,new)={out}  Tr full == Tr new for {len(tr)} primes: {eq}  sample {tr[:3]}')
