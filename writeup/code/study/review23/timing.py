import sys, time, math
sys.path.insert(0, '/home/bas/praime/writeup/code/study')
import numpy as np, r47
RHO0 = 1 / (math.log(1.001e12) * r47.RHO_DEN); NZ = 48; Z = np.exp(2j * np.pi * np.arange(NZ) / NZ)
tot = 0
for D in (80, 160, 240, 320):
    t0 = time.time()
    st = r47.setup(D); M = st['M']; p = RHO0 * st['s1']
    logG = M @ np.log1p(p[:, None] * (Z[None, :] - 1))
    W = p[:, None] * (Z[None, :] - 1) / (1 + p[:, None] * (Z[None, :] - 1))
    for m in range(NZ):
        X = M * W[:, m][None, :]
        V = X.real @ st['R'] + 1j * (X.imag @ st['R']); logG[:, m] += 0.5 * np.sum(V * X, axis=1)
    print(D, round(time.time() - t0, 2), 's', flush=True)
