import mpmath as mp
mp.mp.dps = 40
chi = [0, 1, -1, -1, 1]; q = 5
def L(s): return mp.power(q, -s) * mp.fsum(chi[a] * mp.zeta(s, mp.mpf(a) / q) for a in range(1, q))
print('L(0)=', L(0), ' mp.dirichlet(0)=', mp.dirichlet(0, chi))
L1 = mp.diff(L, 0, 1); L2 = mp.diff(L, 0, 2)
B0 = L2 / (2 * L1)
print('B0 = L''(0)/(2L''(0)) =', B0)
# series check: L'/L(s) - 1/s at small s
for s in [mp.mpf('1e-8'), mp.mpf('-1e-8')]:
    print(' check', mp.diff(L, s) / L(s) - 1 / s)
d1 = mp.dirichlet(1, chi, 1) / mp.dirichlet(1, chi); L = (lambda s: mp.dirichlet(s, chi)) if False else L
print("L(1)", mp.dirichlet(1,chi), '2 ln phi/sqrt5', 2*mp.log((1+mp.sqrt(5))/2)/mp.sqrt(5))
print("L'/L(1)=", d1)
Ca = d1 - B0 - mp.log(2)
Cb = 2 * d1 + mp.log(mp.mpf(q) / mp.pi) + mp.digamma(mp.mpf(1) / 2)
print('C (y=1) =', Ca); print('C (Hadamard) =', Cb, ' diff', Ca - Cb)
# T_e(1)
print('T_e(1) series', mp.nsum(lambda k: 1 / ((2*k) * (2*k+1)), [1, mp.inf]), ' 1-ln2', 1 - mp.log(2))
# Functional-equation derivation of B0: for even chi, Lambda(s)=(q/pi)^{s/2}Gamma(s/2)L(s) = eps Lambda(1-s), real chi eps=1
# => L'/L(s) = -L'/L(1-s) - ln(q/pi) - 1/2 psi(s/2) - 1/2 psi((1-s)/2) ; psi(s/2)/2 = -1/s - gamma/2 + O(s)
B0fe = -d1 - mp.log(q/mp.pi) + mp.euler/2 - mp.digamma(mp.mpf(1)/2)/2
print('B0 from functional eq =', B0fe)
# independent L'(1) by numeric differentiation of dirichlet away from 1 (symmetric)
h = mp.mpf('1e-10')
Lp = (mp.dirichlet(1+h, chi) - mp.dirichlet(1-h, chi)) / (2*h)
print("L'/L(1) numeric diff", Lp / mp.dirichlet(1, chi))
