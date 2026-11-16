# B: heuristic survival factor F = prod_{k<=K} (1 - P_k/k), averaged over p mod 5,7,11,13
import numpy as np
from math import log
from sympy import primerange, factorint

# tail constant prod_{r>13} r(r-2)/(r-1)^2, computed directly to 2e7 plus tail estimate
rs = np.array(list(primerange(17, 20_000_000)), dtype=np.float64)
logt = np.sum(np.log1p(-1.0 / (rs - 1.0) ** 2))
# tail r > 2e7: sum log(1 - 1/(r-1)^2) ~ -sum 1/r^2 ~ -1/(R ln R)
R = 2e7
logt += -1.0 / (R * log(R))
TAIL = float(np.exp(logt))
# compare with twin prime constant
C2 = 0.66016181584686957393
small = 1.0
for r in (3, 5, 7, 11, 13):
    small *= r * (r - 2) / (r - 1) ** 2
print("tail computed", TAIL, " via C2", C2 / small)

SMALL = (5, 7, 11, 13)


def model(p_size, K, cls):
    p4 = cls % 4
    p3 = cls % 3
    ks = np.arange(1, K + 1)
    # mod 8 condition on q = 2kp+1 depends on kp mod 4
    q8 = (2 * ((ks * p4) % 4) + 1) % 8
    allowed = (q8 == 1) | (q8 == 7)
    # r = 3 local factor
    f3 = np.where(((2 * ks * p3 + 1) % 3) == 0, 0.0, 1.5)
    # r > 13: TAIL * prod_{r | k, r > 13} (r-1)/(r-2)
    big = np.ones(K)
    for i, k in enumerate(ks):
        for r in factorint(int(k)):
            if r > 13:
                big[i] *= (r - 1) / (r - 2)
    base = np.where(allowed, 2.0 * f3 * TAIL * big / np.log(2.0 * ks * p_size + 1.0), 0.0)
    # residues grid
    grids = np.meshgrid(*[np.arange(1, r) for r in SMALL], indexing="ij")
    grids = [g.ravel() for g in grids]
    logF = np.zeros(grids[0].shape)
    F = np.ones(grids[0].shape)
    for i, k in enumerate(ks):
        if base[i] == 0.0:
            continue
        loc = np.full(grids[0].shape, base[i])
        for r, g in zip(SMALL, grids):
            killed = ((2 * int(k) * g + 1) % r) == 0
            loc = loc * np.where(killed, 0.0, r / (r - 1))
        F = F * (1.0 - loc / k)
    return F.mean(), F.size


for p_size, K in ((1e6, 3000), (2.6e5, 100)):
    vals = {}
    for cls in (1, 5, 7, 11):
        vals[cls], n = model(p_size, K, cls)
    m = np.mean(list(vals.values()))
    print(f"p={p_size:g} K={K} combos={n}")
    print("  absolute", {c: round(v, 4) for c, v in vals.items()})
    print("  relative", {c: round(v / m, 4) for c, v in vals.items()}, "mean", round(m, 5))
