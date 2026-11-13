import math, numpy as np
from cert import sample
from kbm import K0
G = np.load('hp2419.npy')
for (x, a, beta, m, c) in [(1000.0, 12, 12, 2, 1.2), (5000.0, 60, 14, 2, 1.2), (9900.0, 130, 16, 2, 1.2), (150.0, 1.9, 6, 2, 1.2)]:
    T = c*beta*x/a
    z, d, info = sample(x, G, T, a, beta, m, 1e-12)
    k0 = K0(beta, m)
    print(f'x={x} a={a} beta={beta} T={T:.1f} NT={info["NT"]}: delta/K(0) = {d/k0:.3e}; parts/K(0): ' + ', '.join(f'{k} {v/k0:.2e}' for k, v in info.items() if k != 'NT'))
