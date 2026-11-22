"""Independent check of r42.rho_twin for Q = 2310: S({0,2,d,d+2})/S({0,2})^2 from theory2.S1 (general code), divided
by the factors p | 2310, which are (1 - nu/p)/(1 - nu/p)^2 with nu = |{0, 2} mod p|."""
import theory2 as T2, r42
Q = 2310
fac = 1.0
for p in (2, 3, 5, 7, 11):
    nu = len({0 % p, 2 % p}); fac *= (1 - nu / p) / (1 - nu / p) ** 2
R = r42.rho_twin(3000, Q)
worst = 0.0
for k in (1, 2, 13, 17, 169, 221, 1001, 2431, 2999):
    d = Q * k
    v = T2.S1([0, 2, d, d + 2]) / T2.S1([0, 2]) ** 2 / fac
    worst = max(worst, abs(v / R[k] - 1))
print('largest relative difference', worst)
