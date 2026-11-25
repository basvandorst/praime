"""basic checks: own theta vs scipy; own Conrey characters -> Gauss-sum root numbers vs frozen eps; L(1/2+i gamma) = 0 at
stored zeros with own character tables (mpmath); Kronecker tables; fundamental discriminants."""
import sys, math, time
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import numpy as np
from rlib import *
from scipy.special import loggamma
import mpmath as mp

t = np.linspace(95, 1000, 5000)
for a in (0, 1):
    d = np.abs(theta(t, 7, a) - (t / 2 * math.log(7 / math.pi) + np.imag(loggamma(0.25 + a / 2 + 0.5j * t))))
    print('theta own vs scipy, a=%d: max diff %.2e' % (a, d.max()))

C = np.load(FZ + 'pred_C.npz')
labs = [str(x) for x in C['labels']]
Z = load_group(CPLXF)
maxd = 0; bad = []
tabs = {}
for k, lab in enumerate(labs):
    if not lab.startswith('N') or lab.endswith('_conj'): continue
    N, m = int(lab.split('_')[0][1:]), int(lab.split('_')[1][1:])
    tab = conrey_table(N, m); tabs[lab] = tab
    a = 0 if abs(tab[N - 1] - 1) < 1e-9 else 1
    # primitivity: no proper divisor d of N with chi trivial on n = 1 mod d
    prim = True
    for d_ in range(1, N):
        if N % d_ == 0:
            if all(abs(tab[n] - 1) < 1e-9 for n in range(1, N) if math.gcd(n, N) == 1 and n % d_ == 1):
                prim = False
    e = gauss_eps(tab, a)
    d = abs(e - complex(C['eps'][k])); maxd = max(maxd, d)
    if d > 1e-8 or a != int(C['a'][k]) or not prim: bad.append((lab, a, int(C['a'][k]), prim, d))
    ek = complex(C['eps'][k + 1])
    if abs(ek - np.conj(e)) > 1e-8: bad.append((lab, 'conj side eps', ek))
print('complex: max |eps(own Gauss sum, own Conrey table) - frozen eps| = %.2e; problems: %s' % (maxd, bad))

# L(1/2 + i gamma, chi) at stored zeros, own tables (mpmath dirichlet)
mp.mp.dps = 20
t0 = time.time()
rng = np.random.default_rng(5)
for lab in ['N101_m2', 'N124_m3', 'N388_m3', 'N235_m2', 'N173_m2']:
    N = int(lab.split('_')[0][1:]); tab = tabs[lab]
    z = Z[lab]
    chi = [mp.mpc(x.real, x.imag) for x in tab]
    out = []
    for g in (z[z > 0][rng.integers(0, 20)], z[z < 0][-1 - rng.integers(0, 20)], z[(z > 300)][3]):
        v0 = abs(mp.dirichlet(mp.mpf(0.5) + 1j * mp.mpf(g), chi))
        v1 = abs(mp.dirichlet(mp.mpf(0.5) + 1j * mp.mpf(g + 0.05), chi))
        out.append('%.1f: %.1e / %.1e' % (g, v0, v1))
    print(lab, '|L(1/2+ig)| at a stored zero / at g+0.05:', '; '.join(out), '(%.0fs)' % (time.time() - t0), flush=True)

# real: tables are characters mod |D| (periodic, primitive) and D fundamental
P = np.load(FZ + 'pred_OO.npz'); Ds = [int(x) for x in P['Ds']]
nf = [D for D in Ds if not ((D % 4 == 1 and all(e == 1 for e in factor(abs(D)).values())) or (D % 4 == 0 and (D // 4) % 4 in (2, 3) and all(e == 1 for p, e in factor(abs(D // 4)).items())))]
print('OO: %d discriminants; non-fundamental: %s' % (len(Ds), nf))
D = -399; tab = real_table(D)
print('periodicity check kron(-399, n + 399) == kron(-399, n):', all(kron(D, n + 399) == kron(D, n) for n in range(1, 2000)))
z = load_group(R120F)['D-399']
chi = [mp.mpf(int(x)) for x in tab]
print('D-399 |L| at zero / off:', abs(mp.dirichlet(0.5 + 1j * mp.mpf(z[200]), chi)), abs(mp.dirichlet(0.5 + 1j * mp.mpf(z[200] + 0.05), chi)))
