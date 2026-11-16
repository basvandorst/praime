# D + E: own segmented sieve of Eratosthenes over all integers, segments aligned to 24.
# D: byte j = bits (MSB first) of primality of 24j-1, 24j+1, 24j+5, 24j+7, 24j+11, 24j+13, 24j+17, 24j+19
# E: twin pairs (n, n+2), n = 11 mod 12 ("up") and n = 5 mod 12 ("down")
import numpy as np
from math import isqrt

JMAX = 41_666_666
LIMIT = 10**9
TOP = 24 * (JMAX + 1) + 24 + 64
small = np.ones(isqrt(TOP) + 2, dtype=bool)
small[:2] = False
for i in range(2, isqrt(len(small)) + 1):
    if small[i]:
        small[i * i::i] = False
sp = np.nonzero(small)[0]
sp = sp[sp >= 5]  # 2 and 3 never divide the numbers we read (all coprime to 6)

B = 1 << 20
hist = np.zeros(256, dtype=np.int64)
pos = {223: [], 251: []}
up = down = 0
up_all = down_all = 0  # including small pairs
extra = []  # twin pairs with n+2 > 1e9
total_twins_below = 0  # all twin pairs with p+2 < 1e9 (incl. (3,5))
cols = np.array([1, 5, 7, 11, 13, 17, 19])
J0 = 1
while J0 <= JMAX:
    nb = min(B, JMAX - J0 + 1)
    base = 24 * (J0 - 1)
    L = 24 * (nb + 1) + 24  # one extra row of slack for n+2
    seg = np.ones(L, dtype=bool)
    for r in sp:
        r = int(r)
        if r * r >= base + L:
            break
        st = max(r * r, ((base + r - 1) // r) * r)
        seg[st - base::r] = False
    if base == 0:
        seg[0] = seg[1] = False
        # 2 and 3 are prime but irrelevant: leave as is (never read for bytes; twins handled below)
    rows = seg[: 24 * (nb + 1)].reshape(nb + 1, 24)
    b = rows[:nb, 23].astype(np.int64) << 7
    for sh, c in zip(range(6, -1, -1), cols):
        b |= rows[1:nb + 1, c].astype(np.int64) << sh
    hist += np.bincount(b, minlength=256)
    for v in (223, 251):
        idx = np.nonzero(b == v)[0]
        pos[v].extend((J0 + idx).tolist())
    # twins: owned n in [base, base + 24*nb)
    own = 24 * nb
    n = base + np.arange(own)
    tw = seg[:own] & seg[2:own + 2]
    m12 = n % 12
    sel_u = tw & (m12 == 11)
    sel_d = tw & (m12 == 5)
    nu = n[sel_u]
    nd = n[sel_d]
    up += int(np.count_nonzero((nu >= 23) & (nu + 2 <= LIMIT)))
    down += int(np.count_nonzero((nd >= 23) & (nd + 2 <= LIMIT)))
    up_all += int(np.count_nonzero(nu + 2 <= LIMIT))
    down_all += int(np.count_nonzero(nd + 2 <= LIMIT))
    extra += [int(x) for x in np.concatenate([nu, nd]) if x + 2 > LIMIT]
    J0 += nb

print("D hist[223]", hist[223], "hist[251]", hist[251], "hist[255]", hist[255])
print("D 223 positions", pos[223])
print("D 251 positions", pos[251])
print("D first 251 j", pos[251][0] if pos[251] else None, "-> p =", 24 * pos[251][0] - 1 if pos[251] else None)
print("D up to j=416666: 223:", sum(1 for j in pos[223] if j <= 416666), " 251:", sum(1 for j in pos[251] if j <= 416666))
print("E up", up, "down", down, "sum", up + down)
print("E up_all (incl 11,13)", up_all, "down_all (incl 5,7 and 17,19)", down_all, "+ (3,5) =", up_all + down_all + 1)
print("E twin pairs with n+2 > 1e9 within sieve range (n < 24*(JMAX+1)):", extra)
np.save("hist.npy", hist)
