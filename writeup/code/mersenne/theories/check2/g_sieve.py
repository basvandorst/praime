# Item G, stage 1: segmented odd-only sieve below 1e9.
# Saves the primes in each class 1, 5, 7, 11 mod 12 (uint32, sorted) and a packed odd-only primality bit array.
import numpy as np, math, time, os
X = 10**9                      # sieve x < X
HALF = X // 2                  # odd x = 2i+1, i in [0, HALF)
B = 1 << 24                    # odd numbers per segment (multiple of 8)
t0 = time.time()
lim = math.isqrt(X) + 1
bs = np.ones(lim + 1, bool); bs[:2] = False
for i in range(2, math.isqrt(lim) + 1):
    if bs[i]: bs[i*i::i] = False
base = [int(p) for p in np.nonzero(bs)[0] if p > 2]
cls = {a: [] for a in (1, 5, 7, 11)}
packed = []
total = 0
for i0 in range(0, HALF, B):
    i1 = min(HALF, i0 + B)
    seg = np.ones(i1 - i0, bool)
    if i0 == 0: seg[0] = False            # x = 1
    xlo = 2 * i0 + 1; xhi = 2 * (i1 - 1) + 1
    for p in base:
        if p * p > xhi: break
        x0 = max(p * p, ((xlo + p - 1) // p) * p)
        if x0 % 2 == 0: x0 += p
        if x0 > xhi: continue
        seg[(x0 - 1) // 2 - i0::p] = False
    packed.append(np.packbits(seg))
    x = (2 * (np.nonzero(seg)[0] + i0) + 1).astype(np.uint32)
    total += len(x)
    r = x % 12
    for a in cls:
        cls[a].append(x[r == a])
out = os.path.dirname(os.path.abspath(__file__))
np.save(os.path.join(out, 'G_bits.npy'), np.concatenate(packed))
for a in cls:
    arr = np.concatenate(cls[a]); cls[a] = None
    np.save(os.path.join(out, 'G_p%d.npy' % a), arr)
    print('class', a, len(arr))
print('odd primes < 1e9:', total, ' (pi(1e9) - 1 = 50847533)', ' time %.1f s' % (time.time() - t0))
