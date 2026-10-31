from sympy import primerange, n_order
def revdigits(p, L):
    # 10-adic digits of -1/p, least significant first = period of 1/p read backwards
    m = 10**L; Y = (-pow(p, -1, m)) % m
    return [(Y // 10**k) % 10 for k in range(L)], Y
def simulate(b, N, J, L):
    d = b[:N]; c = 0
    for n in range(N, L):
        v = d[n-J] + d[n-N] + c; d.append(v % 10); c = v // 10
    return d
if __name__=="__main__":
    bad = 0; tested = 0
    for p in primerange(7, 1500):
        T = n_order(10, p); L = 3*T + 130
        b, Y = revdigits(p, L)
        for N in range(2, 41):
            for J in range(1, N):
                ok_sim = simulate(b, N, J, L) == b
                div = (10**N + 10**J - 1) % p == 0
                crit = div and ((10**N + 10**J - 1)//p + 10**J * (Y % 10**(N-J))) < 10**N
                tested += 1
                if ok_sim != crit: bad += 1; print("mismatch", p, N, J, ok_sim, crit)
    print("tested", tested, "mismatches", bad)
    b,_ = revdigits(7, 12); print(b, simulate(b,3,2,12)==b)
