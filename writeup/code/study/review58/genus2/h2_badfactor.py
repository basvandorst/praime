"""H2, bad primes: is the functional-equation check sensitive to the Euler factor at p = 2? For one curve with 2 | N, build the
L-function with the factor at 2 replaced (removed; X-coefficient sign flipped) and compare lfuncheckfeq with the original.
Also the share of c carried by the bad primes, for the 51 curves of H1."""
import sys, os, math, time, glob; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, cypari2
from g2data import load
pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6)
rows = {r['key']: r for r in load()}
key = sys.argv[1] if len(sys.argv) > 1 else 'G2_N2862_17'; r = rows[key]
L = pari('my(L = lfungenus2([%s,%s])); L[6] = lfunrootres(L)[3]; L' % (r['spec'][1], r['spec'][2]))
N, k, V = pari.lfunparams(L); eps = pari.lfunrootres(L)[2]
P2 = 1 / pari.lfuneuler(L, 2); print(key, 'N', N, 'eps', eps, 'Euler factor at 2: 1/(%s)' % P2)
t = time.time(); print('  original: lfuncheckfeq %s (%.0f s)' % (pari.lfuncheckfeq(L), time.time() - t), flush=True)
nmax = 40000
a = pari.lfunan(L, nmax)
for lab, Pn in (('factor at 2 removed', pari('1 + 0*x')), ('top coefficient flipped', pari.subst(P2, 'x', pari('-x')) if False else None)):
    if Pn is None:
        co = pari.Vecrev(P2); co[len(co) - 1] = -co[len(co) - 1]; Pn = pari.Polrev(co)
    # B(s) = A(s) * P2(2^-s) / Pn(2^-s): coefficients via the 2-part
    vo = [0] * (nmax + 1)
    Pc = [complex(c).real for c in pari.Vecrev(P2)]; Q = [complex(c).real for c in pari.Vecrev(Pn)]
    J = int(math.log2(nmax)) + 1; inv = [1.0] + [0.0] * J
    for j in range(1, J + 1): inv[j] = -sum(Q[i] * inv[j - i] for i in range(1, min(j, len(Q) - 1) + 1))
    corr = [sum(Pc[i] * inv[j - i] for i in range(0, min(j, len(Pc) - 1) + 1)) for j in range(J + 1)]   # coefficients in 2^-js
    b = [0.0] * (nmax + 1)
    for n in range(1, nmax + 1):
        m, e = n, 0
        while m % 2 == 0: m //= 2; e += 1
        # b(n) = sum_{j <= e} corr_j a(n / 2^j)
        b[n] = sum(corr[j] * float(a[n // 2 ** j - 1]) for j in range(e + 1))
    vec = pari('Vec')([round(x) for x in b[1:]]) if all(abs(x - round(x)) < 1e-9 for x in b[1:]) else pari('Vec')(b[1:])
    Lb = pari.lfuncreate([vec, 0, V, k, N, eps])
    t = time.time(); print('  %s: lfuncheckfeq %s (%.0f s)' % (lab, pari.lfuncheckfeq(Lb), time.time() - t), flush=True)
# share of c from bad primes
sh = []
for fn in glob.glob('cown/*.npz'):
    z = np.load(fn); p = z['p']; N = int(z['N']); bad = np.array([N % int(q) == 0 for q in p])
    sh.append((os.path.basename(fn)[:-4], float(np.exp(z['logG'][bad].sum())), [int(q) for q in p[bad]]))
v = np.array([s[1] for s in sh]); print('product of the bad-prime factors of c over 51 curves: median %.3f, range %.3f..%.3f' % (np.median(v), v.min(), v.max()))
