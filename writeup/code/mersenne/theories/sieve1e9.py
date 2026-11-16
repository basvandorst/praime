"""Segmented sieve to X: wheel-24 bytes (bits 24j-1, 24j+1, 24j+5, 24j+7, 24j+11, 24j+13, 24j+17, 24j+19, MSB first),
byte counts, consecutive-byte pair counts, the j of bytes 223 and 251, twin counts per block of 1000 bytes
(up both: 12i-1, 12i+1; down both: 12i+5, 12i+7). usage: sieve1e9.py X out.npz"""
import sys, time, numpy as np
X = int(float(sys.argv[1])); out = sys.argv[2]
J = X // 24
r = int((24 * J + 20) ** 0.5) + 2
s = np.ones(r + 1, bool); s[:2] = False
for i in range(2, int(r ** 0.5) + 1):
    if s[i]: s[i * i::i] = False
base = np.nonzero(s)[0].astype(np.int64)
base = base[base > 3]
M = 500_000
OFF = np.array([-1, 1, 5, 7, 11, 13, 17, 19])
W = (1 << np.arange(7, -1, -1)).astype(np.int64)
cnt = np.zeros(256, np.int64); pair = np.zeros(65536, np.int64)
rare = []; twu = []; twd = []
prev = None; t0 = time.time()
for j0 in range(1, J + 1, M):
    j1 = min(j0 + M, J + 1)
    lo = 24 * j0 - 1; hi = 24 * (j1 - 1) + 20
    a = np.ones(hi - lo, bool)
    for p in base:
        if p * p >= hi: break
        st = max(p * p, ((lo + p - 1) // p) * p)
        a[st - lo::p] = False
    idx = (24 * np.arange(j0, j1, dtype=np.int64))[:, None] + OFF[None, :] - lo
    bits = a[idx]
    b = (bits.astype(np.int64) * W).sum(1)
    cnt += np.bincount(b, minlength=256)
    if prev is not None:
        pair[prev * 256 + b[0]] += 1
    pair += np.bincount(b[:-1] * 256 + b[1:], minlength=65536)
    prev = b[-1]
    for v in (223, 251):
        for k in np.nonzero(b == v)[0]:
            rare.append((v, j0 + int(k)))
    up = ((b >> 6) & 3 == 3).astype(np.int64) + ((b >> 2) & 3 == 3)
    dn = ((b >> 4) & 3 == 3).astype(np.int64) + (b & 3 == 3)
    pad = (-len(up)) % 1000
    twu.append(np.pad(up, (0, pad)).reshape(-1, 1000).sum(1)); twd.append(np.pad(dn, (0, pad)).reshape(-1, 1000).sum(1))
    if (j0 - 1) % (20 * M) == 0:
        print(f'j {j0:,} ({24 * j0:.3g}) {time.time() - t0:.0f} s', flush=True)
bu = np.concatenate(twu); bd = np.concatenate(twd); twu = bu; twd = bd
np.savez(out, cnt=cnt, pair=pair, rare=np.array(rare, np.int64).reshape(-1, 2), bu=bu, bd=bd, J=J,
         tu=twu.sum() + 0, td=twd.sum() + 0)
print('done', J, 'bytes', cnt.sum(), 'rare', rare[:40], 'twins up', twu.sum(), 'down', twd.sum(), f'{time.time() - t0:.0f} s')
