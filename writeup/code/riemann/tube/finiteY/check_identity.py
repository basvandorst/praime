# Check A_k(s) = int_0^inf u^k/k! h(u) e^{-su} du  ==  sum_rho r_rho (s-rho)^{-(k+1)}, r = 1/(rho(1-rho)).
import math, numpy as np, mpmath as mp
from sympy import primerange
gam = float(mp.euler); L2P = math.log(2*math.pi); N = 10**7
pp = []
for p in primerange(2, N+1):
    q = p
    while q <= N: pp.append((q, math.log(p))); q *= p
pp.sort(); n = np.array([1.0]+[a for a,_ in pp]); w = np.array([0.0]+[b for _,b in pp])
A = np.cumsum(w/n); psi = np.cumsum(w)
def H(y):
    i = np.searchsorted(n, y, side='right') - 1
    x = np.minimum(1/y, 1-1e-15); T = -0.5*np.log1p(-x*x) - (np.arctanh(x)/x - 1)
    return y*A[i] - psi[i] - y*np.log(y) + (1+gam)*y - L2P + T
du = 1e-5; u = np.arange(0, math.log(N), du) + du/2; h = H(np.exp(u))
print("H(1+) =", float(H(np.array([1.0000001]))[0]), " C =", 2+gam-math.log(4*math.pi))
zs = [complex(mp.zetazero(j)) for j in range(1, 301)]
for k, s in [(20, 3+0j), (20, 3+14.134725j), (30, 3.5+21.02204j), (12, 2+0j)]:
    lw = k*np.log(np.maximum(u, 1e-300)) - math.lgamma(k+1) - s.real*u
    data = np.sum(np.exp(lw)*h*np.exp(-1j*s.imag*u))*du
    # missing tail u > ln N: bound with |h| <= 0.05 e^{u/2}
    zsum = sum(1/(r*(1-r))*(s-r)**(-(k+1)) + 1/(r.conjugate()*(1-r.conjugate()))*(s-r.conjugate())**(-(k+1)) for r in zs)
    print(f"k={k} s={s}: data {data:.6e}  zeros {zsum:.6e}  rel.diff {abs(data-zsum)/abs(zsum):.2e}")
