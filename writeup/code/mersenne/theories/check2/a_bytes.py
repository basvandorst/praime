# Item A: wheel-24 bytes. Independent recomputation.
# Byte j = bits (MSB first) of primality of 24j+o, o in OFFS.
import sys, numpy as np, math
OFFS = [-1, 1, 5, 7, 11, 13, 17, 19]
J = int(sys.argv[1]) if len(sys.argv) > 1 else 41666666

def small_primes(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]

top = 24 * (J + 1) + 19
base = small_primes(int(math.isqrt(top)) + 2)
base5 = [int(p) for p in base if p >= 5]

# wheel sieve: for each offset o, array over j in a block; mark j with p | 24j+o
B = 2_000_000
counts = np.zeros(256, np.int64)
pair = np.zeros(65536, np.int64)
prev_last = None
nprimes = 0
for j0 in range(1, J + 1, B):
    j1 = min(J + 1, j0 + B)   # j in [j0, j1)
    n = j1 - j0
    byte = np.zeros(n, np.uint16)
    for bit, o in enumerate(OFFS):
        isp = np.ones(n, bool)
        for p in base5:
            if p * p > 24 * (j1 - 1) + o:
                break
            inv = pow(24, -1, p)
            r = (-o * inv) % p              # j = r mod p  <=> p | 24j+o
            # start: smallest j >= max(j0, jmin) with j = r mod p, where 24j+o >= p*p
            jmin = max(j0, -(-(p * p - o) // 24))
            js = jmin + ((r - jmin) % p)
            if js < j1:
                isp[js - j0::p] = False
        nprimes += int(isp.sum())
        byte |= (isp.astype(np.uint16) << (7 - bit))
    counts += np.bincount(byte, minlength=256)
    # consecutive pairs (byte_j, byte_{j+1})
    if prev_last is not None:
        pair[int(prev_last) * 256 + int(byte[0])] += 1
    pair += np.bincount(byte[:-1].astype(np.int64) * 256 + byte[1:], minlength=65536)
    prev_last = byte[-1]
    print('block', j0, file=sys.stderr, flush=True)

np.save('A_counts_%d.npy' % J, counts)
np.save('A_pairs_%d.npy' % J, pair)
print('J', J, 'top', 24 * J + 19, 'primes among wheel values', nprimes)
