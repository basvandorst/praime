"""Round 104 (preregister/round104.txt): for every fundamental discriminant 3 <= |D| <= DMAX: first zero gamma_1 of
L(s, chi_D) (PARI lfunzeros), L(1, chi_D), the first split prime P(D); Spearman correlations of the scaled height
g = gamma_1 ln|D| / (2 pi) with ln L(1) and ln P. usage: r104.py [DMAX]"""
import sys, json, math, os, time
import cypari2
from scipy.stats import spearmanr
pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6)
DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
t0 = time.time(); rows = []
for D in list(range(-DMAX, -2)) + list(range(5, DMAX + 1)):
    if not int(pari.isfundamental(D)) or D == 1: continue
    L = pari.lfuncreate(D)
    T = 4.0; z = []
    while not z and T < 200:
        z = [float(t) for t in pari.lfunzeros(L, T)]; T *= 2
    g1 = min(z)
    L1 = float(pari.lfun(L, 1))
    P = next(p for p in map(int, pari.primes(400)) if int(pari.kronecker(D, p)) == 1) if any(int(pari.kronecker(D, p)) == 1 for p in map(int, pari.primes(400))) else 3000
    rows.append({'D': D, 'gamma1': g1, 'L1': L1, 'P': P, 'g': g1 * math.log(abs(D)) / (2 * math.pi)})
out = {'DMAX': DMAX, 'n': len(rows), 'seconds': time.time() - t0, 'rows': rows, 'tests': {}}
for sgn, name in ((-1, 'negative'), (1, 'positive')):
    R = [r for r in rows if (r['D'] < 0) == (sgn < 0)]
    s1 = spearmanr([math.log(r['L1']) for r in R], [r['g'] for r in R])[0]
    s2 = spearmanr([math.log(r['P']) for r in R], [r['g'] for r in R])[0]
    low = sorted(R, key=lambda r: r['g'])[:10]
    out['tests'][name] = {'n': len(R), 'spearman_lnL1_g': s1, 'spearman_lnP_g': s2, 'lowest10': [(r['D'], round(r['gamma1'], 4), round(r['g'], 4), r['P'], round(r['L1'], 4)) for r in low]}
    print(name, len(R), 'rho(lnL1, g) = %.3f  rho(lnP, g) = %.3f' % (s1, s2))
    print('  lowest scaled first zeros (D, gamma1, g, P, L1):', out['tests'][name]['lowest10'])
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'r104.json'), 'w'), indent=1)
print('%d discriminants, %.0f s' % (len(rows), time.time() - t0))
