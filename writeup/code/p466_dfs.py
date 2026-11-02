"""Puzzle 466: all odd B <= X with beta*B - alpha*phi(B) = beta*K (Lehmer-style DFS: a skeleton S of
small primes, then the last prime r from a division, or the last two r < s from a factorisation of
alpha^2 u^2 - D (alpha u + beta K), D = alpha u - beta S, u = phi(S)).
Modes:
  p466_dfs.py validate            DFS vs a phi sieve for 120 equations below 1e7
  p466_dfs.py n7                  pair shapes of n = 7; rediscovers a(7) = 30057431145 (Resta)
  p466_dfs.py n3n9 X              every pair shape of n = 3 and n = 9, B <= X
  p466_dfs.py reduced X           odd 3-free B' <= X with B' - 2 phi(B') = K', K' in {1,3,5,7,9,15,17,21,51}
"""
import sys, math
from fractions import Fraction
import cypari2
pari = cypari2.Pari()
pari.allocatemem(2 * 10**8, silent=True)
from sympy import primerange, isprime

def solve(alpha, beta, K, X, minp, forbid3, need3, maxprime=10**6):
    primes = [p for p in primerange(minp, maxprime) if not (forbid3 and p == 3)]
    out = set()
    nodes = [0]

    def finish(S, u, last):
        D = alpha * u - beta * S
        if D <= 0:
            return
        num = alpha * u - beta * K
        if num > 0 and num % D == 0:
            r = num // D
            if r > last and S * r <= X and isprime(r) and S % r:
                out.add(S * r)
        M = alpha * alpha * u * u - D * (alpha * u + beta * K)
        if M == 0:
            return
        sgn = 1 if M > 0 else -1
        A = abs(M)
        if A > 10**60:
            return
        for d in pari.divisors(A):
            d = int(d)
            e = A // d
            for x, y in ((d, sgn * e), (-d, -sgn * e)):
                if (x + alpha * u) % D or (y + alpha * u) % D:
                    continue
                r, s = (x + alpha * u) // D, (y + alpha * u) // D
                if not (last < r < s) or S * r * s > X:
                    continue
                if isprime(r) and isprime(s):
                    out.add(S * r * s)

    def reachable(S, u, j):
        num, den, prod = beta * S, alpha * u, S
        k = j
        while k < len(primes) and prod * primes[k] <= X:
            q = primes[k]
            num *= q; den *= q - 1; prod *= q
            if num > den:
                return True
            k += 1
        return num > den

    def dfs(S, u, i):
        nodes[0] += 1
        if beta * S <= alpha * u and not reachable(S, u, i):
            return
        last = primes[i - 1] if i else 1
        if not need3 or S % 3 == 0:
            finish(S, u, last)
            if S > 1 and beta * S - alpha * u == beta * K:
                out.add(S)
        for j in range(i, len(primes)):
            p = primes[j]
            if S * p * p > X:
                break
            if need3 and S % 3 and p > 3:
                break
            pe, ph = p, p - 1
            while S * pe * p <= X or S * pe <= X // p:
                S2, u2 = S * pe, u * ph
                if beta * S2 >= alpha * u2:
                    if S2 * (beta * S2 - alpha * u2) > beta * K * S2:
                        break
                dfs(S2, u2, j + 1)
                pe *= p; ph *= p
                if S * pe > X:
                    break

    dfs(1, 1, 0)
    return sorted(b for b in out if beta * b - alpha * (int(pari.eulerphi(b))) == beta * K), nodes[0]


def pair_hits(n, a, c, X):
    pa, pc = int(pari.eulerphi(a)), int(pari.eulerphi(c))
    rho = Fraction(a, pa)
    sol, nodes = solve(rho.numerator, rho.denominator, a + n, X, 3 if a % 3 == 0 or n % 3 else 5,
                       a % 3 != 0 and n % 3 == 0, a % 3 == 0 and n % 3 == 0)
    hits = []
    for B in sol:
        N = int(pari.eulerphi(B))
        if N % pa or N % pc:
            continue
        p, q = N // pa + 1, N // pc + 1
        if isprime(p) and isprime(q) and a * p == B - n and c * q == B + n:
            hits.append(B)
    return sol, nodes, hits


def validate():
    import numpy as np
    N = 10**7
    phi = np.arange(N + 1, dtype=np.int64)
    sv = np.ones(N + 1, dtype=np.bool_); sv[:2] = False
    for i in range(2, int(N**0.5) + 1):
        if sv[i]: sv[i*i::i] = False
    for p in np.nonzero(sv)[0]:
        p = int(p); v = phi[p::p]; v -= v // p
    B = np.arange(N + 1, dtype=np.int64)
    bad = tot = 0
    for alpha, need3 in ((2, 0), (3, 1)):
        for K in range(1, 120, 2):
            mask = (B - alpha * phi == K) & (B % 2 == 1) & (B > 1)
            if need3: mask &= (B % 3 == 0)
            sieve = set(np.nonzero(mask)[0].tolist())
            sol, _ = solve(alpha, 1, K, N, 3, False, bool(need3))
            tot += len(sieve)
            bad += set(sol) != sieve
    print('120 equations below 1e7: %d sieve solutions, %d mismatches' % (tot, bad))


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'validate':
        validate()
    elif mode == 'n7':
        for a, c in ((2, 16), (14, 28), (98, 112)):
            sol, nodes, hits = pair_hits(7, a, c, 10**11)
            print('n=7 (%d,%d): %d B-solutions <= 1e11, hits %s' % (a, c, len(sol), hits))
    elif mode == 'n3n9':
        X = int(float(sys.argv[2]))
        for n, a, c in ((3, 2, 8), (3, 6, 12), (3, 12, 18), (3, 18, 24), (3, 48, 54), (9, 6, 24),
                        (9, 18, 36), (9, 36, 54), (9, 54, 72), (9, 144, 162)):
            sol, nodes, hits = pair_hits(n, a, c, X)
            print('n=%d (%d,%d): %d B-solutions <= %.0e, hits %s' % (n, a, c, len(sol), X, hits), flush=True)
    elif mode == 'reduced':
        X = int(float(sys.argv[2]))
        for K in (1, 3, 5, 7, 9, 15, 17, 21, 51):
            sol, nodes = solve(2, 1, K, X, 5, True, False)
            print("B' - 2 phi(B') = %d, B' odd 3-free <= %.0e: %s" % (K, X, sol), flush=True)
