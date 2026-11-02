import numpy as np, gmpy2, sys
N = int(sys.argv[1])
ok = np.ones(N + 1, dtype=np.bool_)
ok[0] = False
ok[3::2] = False
lim = N
sv = np.ones(lim + 1, dtype=np.bool_); sv[:2] = False
for i in range(2, int(lim**0.5) + 1):
    if sv[i]: sv[i*i::i] = False
for q in np.nonzero(sv)[0]:
    q = int(q)
    if q % 4 != 1:
        continue
    c = 2
    while gmpy2.powmod(c, (q - 1) // 2, q) != q - 1:
        c += 1
    s = int(gmpy2.powmod(c, (q - 1) // 4, q))
    for t in (s, q - s):
        start = t if t * t + 1 != q else t + q
        ok[start::q] = False
idx = np.nonzero(ok)[0]
for e in range(1, 2 * len(str(N))):
    x = 10 ** e
    n = int(np.searchsorted(idx, int(gmpy2.isqrt(x - 1)), side='right'))
    if 10 ** e <= (N * N + 1):
        print(e, n)
