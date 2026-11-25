"""consec_null.py: the nulls N1 and N2 as numbers, for a range (X0, X1].

N1: labels of p_n and p_(n+k) independent with the Chebotarev (Sato-Tate) marginals.
N2: the abelian projections (chi(p_n), chi(p_(n+k))) follow Lemke Oliver-Soundararajan (LOS, arXiv 1603.03720v4,
    Main Conjecture p. 2-3, Conjecture 1.3 and Conjecture 1.4); the non-abelian part of each label is drawn independently
    inside its fibre with its Chebotarev probability.  For a real character chi of prime conductor q (and for q = 4):
      sum_{p_n <= x} chi(p_n) chi(p_(n+1)) = -li(x) log(2 pi log x / q) / (2 log x) + O(x / (log x)^(11/4))     (k = 1)
      sum_{p_n <= x} chi(p_n) chi(p_(n+k)) = -li(x) / (2 (k - 1) log x)       + O(x / (log x)^(11/4))     (k >= 2)
    (k = 1: projection of c1 = (phi/2)(1/phi - delta_ab) and c2 (1.1)-(1.2) on chi x chi, = Conjecture 1.3; k >= 2: Conj. 1.4.)
    Over a range the cumulative forms are differenced: b_k = [G(X1) - G(X0)] / [li(X1) - li(X0)].
    P(chi_n = s, chi_(n+k) = t) = (1 + s t b_k)/4 (marginals 1/2; the one-point Chebyshev bias is O(x^(-1/2) log x)).
usage: consec_null.py X0 X1  -> prints b_k and the N1 / N2 tables per scheme (JSON on stdout with --json)."""
import sys, json, math
import numpy as np
from mpmath import li as mli

MARG = {'C23': [1 / 2, 1 / 6, 1 / 3], 'C31': [1 / 2, 1 / 6, 1 / 3], 'C47': [1 / 2, 1 / 10, 1 / 5, 1 / 5],
        'G': [1 / 2] + [1 / 8] * 4, 'E37m2': [1 / 2, 1 / 6, 1 / 3], 'E37m3': [1 / 8, 3 / 16, 3 / 16, 1 / 4, 1 / 8, 1 / 8],
        'C39': [1 / 2, 1 / 8, 1 / 8, 1 / 4], 'E37s': [1 / 2, 1 / 2, 0.0]}
CHI = {'C23': ([-1, 1, 1], 23), 'C31': ([-1, 1, 1], 31), 'C47': ([-1, 1, 1, 1], 47), 'G': ([-1, 1, 1, 1, 1], 4),
       'E37m2': ([-1, 1, 1], 37), 'E37m3': ([1, 1, 1, -1, -1, -1], 3)}

def li(x): return float(mli(max(x, 3.0)))
def G(x, k, q):
    x = max(x, 3.0); L = math.log(x)
    return -li(x) * math.log(2 * math.pi * L / q) / (2 * L) if k == 1 else -li(x) / (2 * (k - 1) * L)
def bk(X0, X1, k, q): return (G(X1, k, q) - G(X0, k, q)) / (li(X1) - li(X0))

def tables(X0, X1):
    out = {}
    for nm, P in MARG.items():
        P = np.array(P); ent = {'N1': {str(k): np.outer(P, P).tolist() for k in (1, 2, 3)}}
        if nm in CHI:
            ch = np.array(CHI[nm][0]); q = CHI[nm][1]; ent['b_LOS'] = {}; ent['N2'] = {}
            for k in (1, 2, 3):
                b = bk(X0, X1, k, q); ent['b_LOS'][str(k)] = b
                # P(a, b) = P(chi_a, chi_b) * P(a | chi_a) * P(b | chi_b),  P(a | chi) = 2 P(a) (marginal of chi = 1/2)
                T = np.array([[(1 + ch[i] * ch[j] * b) / 4 * (2 * P[i]) * (2 * P[j]) for j in range(len(P))] for i in range(len(P))])
                ent['N2'][str(k)] = T.tolist()
        out[nm] = ent
    return out

if __name__ == '__main__':
    X0, X1 = float(sys.argv[1]), float(sys.argv[2])
    t = tables(X0, X1)
    if '--json' in sys.argv: print(json.dumps(t)); sys.exit()
    print('range (%.3g, %.3g]: lambda at midpoint (log) %.3f' % (X0, X1, math.log(math.sqrt(max(X0, 3) * X1))))
    for nm, e in t.items():
        if 'b_LOS' in e: print('%-6s q = %2d  b_1 = %+.5f  b_2 = %+.5f  b_3 = %+.5f' % (nm, CHI[nm][1], e['b_LOS']['1'], e['b_LOS']['2'], e['b_LOS']['3']))
