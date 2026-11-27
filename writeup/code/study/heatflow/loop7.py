"""Loop 7: backward flow (tau < 0) of Haar matrices: share of zeros that left the unit circle, against the two-body law."""
import sys, numpy as np
from scipy.stats import unitary_group
N = 100; n_mat = int(sys.argv[1]); taus = [float(a) for a in sys.argv[2:]]
M = 128; G = 512 * M; xg = np.exp(2j * np.pi * np.arange(M) / M); j = np.arange(N + 1); th = 2 * np.pi * np.arange(G) / G
lost = {T: [] for T in taus}
rng = np.random.default_rng(31)
for b in range(n_mat // 500):
  for U in unitary_group.rvs(N, size=500, random_state=rng):
    z = np.linalg.eigvals(U)
    e = np.fft.fft(np.prod(1 + xg[:, None] * z[None, :], axis=1))[: N + 1] / M
    for T in taus:
        a = e * np.exp(-(T / N ** 2) * j * (N - j))
        r = (np.fft.ifft(a, n=G) * G * np.exp(-0.5j * N * th) / np.sqrt(e[N])).real
        lost[T].append(N - np.count_nonzero(r * np.roll(r, -1) < 0))
print("   tau     share off the circle     two-body, small-gap law   ratio")
for T in taus:
    v = np.array(lost[T]) / N; m, s = v.mean(), v.std(ddof=1) / np.sqrt(n_mat)
    f = 2 / (9 * np.pi) * (2 * abs(T)) ** 1.5
    print(f" {T:6.2f}   {m:9.6f} +- {s:8.6f}     {f:9.6f}              {m/f:6.3f} +- {s/f:5.3f}")
