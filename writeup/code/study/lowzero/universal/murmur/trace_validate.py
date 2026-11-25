"""Validate the k = 2 class-number formula for sum_f eps_f a_f(p) (N squarefree, p prime, p not dividing N)
    T(p,N) = (1/2) H(4pN) + sum_{r >= 1, r^2 N < 4p} H(4pN - r^2 N^2) - (p+1)
(Skoruppa-Zagier as quoted by Zubrilina Sec. 2 and Martin 2505.01093 eq. (2.1)) against PARI:
sum_f eps_f a_f(p) = -trace(T_p * W_N) on S_2^new(N) (eps_f = -(W_N eigenvalue) for k = 2), with T_p from mfheckemat and
W_N from mfatkininit (constant C = 1 checked)."""
import cypari2
pari = cypari2.Pari()
pari.allocatemem(600 * 10 ** 6)


def T_formula(p, N):
    s = pari('qfbhclassno(%d)' % (4 * p * N)) / 2
    r = 1
    while r * r * N < 4 * p:
        s += pari('qfbhclassno(%d)' % (4 * p * N - r * r * N * N))
        r += 1
    return s - (p + 1)


tot = bad = 0
for N, pmax in ((389, 1700), (390, 1700), (431, 1800), (101, 1300)):
    pari('mf = mfinit([%d, 2], 0); A = mfatkininit(mf, %d);' % (N, N))
    C = pari('A[3]')
    nb = nn = 0
    withr = 0
    for p in pari('primes([2, %d])' % pmax):
        p = int(p)
        if N % p == 0:
            continue
        val = -pari('trace(mfheckemat(mf, %d) * A[2])/A[3]' % p)
        tf = T_formula(p, N)
        nn += 1
        nb += (val != tf)
        withr += (4 * p > N)
    tot += nn
    bad += nb
    print(f"N={N}: dim S_2^new = {int(pari('mfdim(mf)'))}, C = {C}; primes p <= {pmax}: {nn} (with r-terms: {withr}); mismatches: {nb}")
print("TOTAL checked:", tot, " mismatches:", bad)
