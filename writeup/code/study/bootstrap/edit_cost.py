# The conductor cost of a finite edit of the primes.  Claim (L = infinity, under RH): for a system P that
# differs from the primes in finitely many points,  lambda_min(P) = 2 inf_t h_P(t),
# h_P(t) = Re sum_n (Lambda(n) - Lambda_P(n)) n^(-1/2 - it),  so the smallest conductor is exp(-2 inf h_P).
# Here: finite supports X, float64, compared with the limit.
import numpy as np
from sympy import primerange
from conductor import FormG
def term(x, t):                      # Re log x / (x^(1/2+it) - 1): a point x with all its powers
    return (np.log(x)/(x**(0.5 + 1j*t) - 1)).real
t = np.linspace(0, 4000, 4000001)
edits = [('remove 2', [2], []), ('remove 3', [3], []), ('remove 7', [7], []), ('add 6', [], [6.0]),
         ('move 5 to 5.5', [5], [5.5]), ('move 2 to 2.1', [2], [2.1])]
for name, rem, add in edits:
    h = sum(term(x, t) for x in rem) - sum(term(x, t) for x in add)
    lim = 2*h.min()
    row = []
    for X in (8.0, 16.0, 32.0, 64.0):
        N = max(12, int(8 + 4*np.log(X)/np.log(2)))
        F = FormG(X, N, [0])
        P = [float(p) for p in primerange(2, int(X) + 1) if p not in rem] + [a for a in add if a <= X]
        row.append('%.4f' % F.mineig(P, [1]*len(P)))
    print('%-14s lambda_min at X = 8, 16, 32, 64: %s   limit 2 inf h = %.4f   conductor exp(-limit) = %.3f' % (name, ' '.join(row), lim, np.exp(-lim)), flush=True)
print('closed forms: remove p -> -2 log p/(sqrt p + 1): p=2 %.4f, p=3 %.4f, p=7 %.4f; add x -> -2 log x/(sqrt x - 1): x=6 %.4f' % (-2*np.log(2)/(2**.5 + 1), -2*np.log(3)/(3**.5 + 1), -2*np.log(7)/(7**.5 + 1), -2*np.log(6)/(6**.5 - 1)))
