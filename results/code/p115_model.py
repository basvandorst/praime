# Expected number of k-digit squares N = A^2 (resp. N = p^2, p prime) whose decimal digits use at most two values.
# Digit DP over N mod M (M = 83160 = 2^3*3^3*5*7*11) gives the exact count of digit strings in every class r;
# a number in class r is a square with probability rho(r)/(2 sqrt N) (rho = number of square roots of r mod M),
# and a prime square with probability [r in R] * M/|R| / (2 sqrt N ln sqrt N), R = squares of units mod M.
import numpy as np, math, itertools, sys
M = 83160
u = np.arange(M)
rho = np.bincount(u * u % M, minlength=M).astype(float)
units = u[np.gcd(u, M) == 1]
Rmask = np.zeros(M, bool); Rmask[units * units % M] = True; nR = Rmask.sum()
pow10 = [pow(10, i, M) for i in range(80)]
def dp_counts(digits, L):                       # count strings of length L over 'digits' by value mod M (no leading rule)
    c = np.zeros(M); c[0] = 1
    for i in range(L):
        n = np.zeros(M)
        for d in digits: n += np.roll(c, d * pow10[i] % M)
        c = n
    return c
def expected(k, prime):
    tot = 0.0
    for a, b in itertools.combinations(range(10), 2):
        tail_ab = dp_counts((a, b), k - 1); tail_a = dp_counts((a,), k - 1); tail_b = dp_counts((b,), k - 1)
        for lead in (a, b):
            if lead == 0: continue
            other = b if lead == a else a
            # strings with lead and tail over {a,b}, minus the repdigit lead...lead (counted below once)
            c = np.roll(tail_ab, lead * pow10[k - 1] % M) - np.roll(tail_a if lead == a else tail_b, lead * pow10[k - 1] % M)
            N = (lead + 0.5) * 10 ** (k - 1); s = math.sqrt(N)
            tot += (c[Rmask].sum() * M / nR / (2 * s * math.log(s))) if prime else (c * rho).sum() / (2 * s)
    for d in range(1, 10):                       # repdigits
        r = sum(d * pow10[i] for i in range(k)) % M; N = d * (10 ** k - 1) / 9; s = math.sqrt(N)
        tot += (Rmask[r] * M / nR / (2 * s * math.log(s))) if prime else rho[r] / (2 * s)
    return tot
if __name__ == "__main__":
    K = int(sys.argv[1])
    for k in range(2, K + 1):
        print(k, f"all squares: {expected(k, False):.3f}", f"prime squares: {expected(k, True):.4f}", flush=True)
