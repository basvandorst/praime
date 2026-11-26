# Check omega_f = 12/((k-1) N Y'_f), Y'_f = L(1,Sym^2 f)/zeta(2), at levels with one rational newform (k=2),
# against the Petersson formula: omega_f = 1 - 2 pi sum_{c = N m} S(1,1;c)/c J_1(4 pi/c)  (single form, i^{-2} = -1).
import cypari2, math, cmath
pari = cypari2.Pari(); pari.allocatemem(300*10**6); pari.set_real_precision(30)
def kloost(a, b, c):
    s = 0.0
    for x in range(1, c):
        if math.gcd(x, c) == 1:
            xi = pow(x, -1, c)
            s += math.cos(2 * math.pi * (a * x + b * xi) / c)
    return s
for N, lab in ((11, "[0,-1,1,-10,-20]"), (19, "[0,1,1,-9,-15]"), (17, "[1,-1,1,-1,-14]")):  # dim S_2(Gamma_0(N)) = 1
    E = pari(f"ellinit({lab})"); assert int(pari.ellglobalred(E)[0]) == N
    L = pari.lfunsympow(E, 2)
    Lval = float(pari.lfun(L, 2))
    Yp = Lval / (math.pi**2 / 6)
    om_claim = 12 / (N * Yp)
    tot = 0.0
    for m in range(1, 1500):
        c = N * m
        tot += kloost(1, 1, c) / c * float(pari.besselj(1, 4 * math.pi / c))
    om_pet = 1 - 2 * math.pi * tot
    print(N, lab, "L(1,Sym2)=", round(Lval, 8), "12/(N Y') =", round(om_claim, 6), " Petersson omega =", round(om_pet, 6))
