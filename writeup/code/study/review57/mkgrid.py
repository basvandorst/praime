import numpy as np, math
k0 = int(math.floor(2048 * math.log(1e8)))
ks = np.arange(k0, 43862)
X = np.array([math.floor(math.exp(k / 2048)) for k in ks], dtype=np.int64)
np.savez('grid_disc.npz', X=X, ks=ks)
print(len(X), X[0], X[-1])
