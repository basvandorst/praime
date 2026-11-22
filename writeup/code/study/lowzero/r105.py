"""Round 105 (preregister/round105.txt): race bias between chi_D(p) = -1 and +1 primes up to 1e5 against the first zero
and the exact sum S of round 104. usage: r105.py"""
import json, math, os
import numpy as np
from sympy import primerange
from scipy.stats import spearmanr
import cypari2
pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)
H = os.path.dirname(os.path.abspath(__file__))
rows = json.load(open(os.path.join(H, 'r104.json')))['rows']
P = np.array(list(primerange(3, 10 ** 5)), dtype=np.int64)
def jacobi_vec(D, P):
    # Kronecker (D/p) for odd primes p via Euler's criterion, vectorised with Python ints in chunks
    Dm = D % P
    out = np.array([pow(int(a), (int(p) - 1) // 2, int(p)) for a, p in zip(Dm, P)], dtype=np.int64)
    return np.where(out == 0, 0, np.where(out == 1, 1, -1))
for r in rows:
    D = r['D']
    chi = jacobi_vec(D, P)
    chi = np.concatenate([[int(pari.kronecker(D, 2))], chi])            # p = 2 first
    lead = np.cumsum(-chi)                                              # (# -1) - (# +1)
    m = chi != 0
    r['bias'] = float(np.mean(lead[m] > 0))
    L = pari.lfuncreate(D); a = 1 if D < 0 else 0
    r['S'] = 2 * float(pari.lfun(L, 1, 1)) / r['L1'] + math.log(abs(D) / math.pi) + float(pari.psi((1 + a) / 2.0))
out = {}
for name, sel in (('negative', lambda r: r['D'] < 0), ('positive', lambda r: r['D'] > 0)):
    R = [r for r in rows if sel(r)]
    s1 = spearmanr([r['bias'] for r in R], [r['g'] for r in R])[0]
    s2 = spearmanr([r['bias'] for r in R], [r['S'] for r in R])[0]
    low = sorted(R, key=lambda r: r['g'])[:100]
    med = float(np.median([r['bias'] for r in low])); medall = float(np.median([r['bias'] for r in R]))
    e163 = [r['bias'] for r in R if r['D'] == -163]
    out[name] = {'spearman_bias_g': s1, 'spearman_bias_S': s2, 'median_bias_lowest100': med, 'median_bias_all': medall, 'bias_-163': e163}
    print(name, out[name])
json.dump({'tests': out, 'rows': [{'D': r['D'], 'bias': r['bias'], 'S': r['S'], 'g': r['g']} for r in rows]}, open(os.path.join(H, 'r105.json'), 'w'), indent=0)
