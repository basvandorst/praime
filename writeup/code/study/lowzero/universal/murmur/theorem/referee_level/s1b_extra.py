"""Extra: trace formula at k = 12 (old forms from level 1 and from proper divisors present) and k = 6;
dimension deviation for omega(N) = 4 with all q = 11 mod 12."""
import sys
sys.path.insert(0, '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/referee_level')
from fractions import Fraction as Fr
import cypari2
gp = cypari2.Pari(); gp.allocatemem(400 * 10 ** 6)
def H(n): return Fr(str(gp(f'qfbhclassno({n})')))
def U(m, x2):
    c = [[1], [0, 2]]
    for n in range(1, m):
        a = [0] + [2 * t for t in c[n]]; b = c[n - 1] + [0] * (len(a) - len(c[n - 1]))
        c.append([a[i] - b[i] for i in range(len(a))])
    return sum(t * x2 ** (i // 2) for i, t in enumerate(c[m]) if t)
def rhs(N, k, p):
    s = Fr(1, 2) * H(4 * p * N); r = 1
    while r * r * N * N < 4 * p * N:
        s += (-1) ** (k // 2 - 1) * U(k - 2, Fr(r * r * N, 4 * p)) * H(4 * p * N - r * r * N * N); r += 1
    return s - (p + 1 if k == 2 else 0)
tot = bad = 0
for N, k in [(6, 12), (10, 12), (15, 12), (30, 6), (14, 10)]:
    gp(f'mf = mfinit([{N},{k}],0); L = mfeigenbasis(mf); W = mfatkineigenvalues(mf,{N}); F = mffields(mf)')
    gp(f'mfc = mfinit([{N},{k}],1)')
    nold = int(gp('mfdim(mfc)')) - int(gp('mfdim(mf)'))
    for p in gp('primes([2,60])'):
        p = int(p)
        if N % p == 0: continue
        s = Fr(0)
        for i in range(1, int(gp('#L')) + 1):
            w = int(gp(f'W[{i}][1]')); deg = int(gp(f'poldegree(F[{i}])'))
            a = gp(f'my(c=mfcoef(L[{i}],{p})); if(type(c)=="t_POLMOD", trace(c), {deg}*c)')
            s += w * Fr(str(a))
        lhs = (-1) ** (k // 2) * s / Fr(p) ** ((k - 2) // 2)
        tot += 1; bad += (lhs != rhs(N, k, p))
    print(f'N={N} k={k}: dim new {int(gp("mfdim(mf)"))}, old {nold}; cumulative pairs {tot}, mismatches {bad}', flush=True)
for N in [11 * 23 * 47 * 59, 2 * 11 * 23 * 47]:
    d = int(gp(f'mfdim([{N},2],0)')); phi = int(gp(f'eulerphi({N})')); om = int(gp(f'omega({N})'))
    dev = Fr(d) - Fr(phi, 12)
    print(f'N={N}: dim S_2^new = {d}, phi/12 = {float(Fr(phi,12)):.3f}, dev = {dev} ({float(dev):.3f}), 2^(omega-1) = {2**(om-1)}')
