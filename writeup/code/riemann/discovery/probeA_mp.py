"""Independent check of the near-central zero: L(1/2, chi_590) mod 395,137 with mpmath's own Dirichlet L-function."""
import mpmath as mp
from sympy import primitive_root
mp.mp.dps = 20
q = 395137; N = q - 1; k = 590
g = int(primitive_root(q)); ind = [0] * q; v = 1
for j in range(N):
    ind[v] = j; v = v * g % q
w = mp.exp(2j * mp.pi / N)
chi = [0] + [w ** ((k * ind[a]) % N) for a in range(1, q)]
for s in (mp.mpf(0.5), mp.mpc(0.5, 1.655759586864e-6)):
    print(f'mpmath: L({mp.nstr(s, 15)}, chi_590) = {mp.nstr(mp.dirichlet(s, chi), 12)}', flush=True)
