import numpy as np
from scipy.stats import unitary_group
rng = np.random.default_rng(1)
N = 812; ls = np.array([0.5, 1, 2, 4, 6, 7, 7.5, 8])
res = []
for rep in range(24):
    U = unitary_group.rvs(N, random_state=rng)
    th = np.sort(np.angle(np.linalg.eigvals(U)) % (2*np.pi)) * N / (2*np.pi)
    u = np.concatenate([th - N, th, th + N])  # periodic
    out = []
    for l in ls:
        x = np.arange(0, N, 0.01)
        c = np.searchsorted(u, x + l) - np.searchsorted(u, x)
        out.append(c.var())
    res.append(out)
res = np.array(res)
print('CUE N=812 mean Sigma2', np.round(res.mean(0), 4))
print('CUE N=812 std  Sigma2', np.round(res.std(0, ddof=1), 4))
