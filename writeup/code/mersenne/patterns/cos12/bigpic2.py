"""Big pictures of the cosine rows, period 12, drawn so that the structure survives zooming out.
n sits in row (n-1)//w, column (n-1)%w, lifted by -cos(2 pi n / 12) * amp.
Primes: blue in the up lane (n = 1 or 11 mod 12), orange in the down lane (n = 5 or 7 mod 12), green for twins
(mode twins), pink for Mersenne exponents. Composites: faint (fade 0..1, 0 = hidden).
usage: bigpic2.py N w sx rowh mode fade out.png"""
import sys, zlib, struct, numpy as np
N, w, sx, rowh, mode, fade, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), sys.argv[5], float(sys.argv[6]), sys.argv[7]
P = 12; amp = rowh * 0.34; dot = max(2, min(sx, int(round(rowh))) - 0) if mode != 'thin' else 2
EXP = [2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,216091,756839,859433]
s = np.ones(N + 3, bool); s[:2] = False
for i in range(2, int((N + 2) ** .5) + 1):
    if s[i]: s[i * i::i] = False
n = np.arange(1, N + 1); rows = (N + w - 1) // w
pad = int(rowh) + 4
H, W = int(rows * rowh) + 2 * pad, w * sx + 2 * sx
BG = np.array((11, 18, 32), float)
img = np.zeros((H, W, 3), np.uint8); img[:] = BG.astype(np.uint8)
x = sx + ((n - 1) % w) * sx + sx // 2
y = np.round(pad + ((n - 1) // w) * rowh + rowh / 2 - np.cos(2 * np.pi * n / P) * amp).astype(int)
def paint(mask, col, d):
    xs, ys = x[mask], y[mask]
    for dx in range(d):
        for dy in range(d):
            img[np.clip(ys + dy - d // 2, 0, H - 1), np.clip(xs + dx - d // 2, 0, W - 1)] = col
pr = s[n]; up = np.cos(2 * np.pi * n / P) > 0
twin = pr & (s[n + 2] | s[np.maximum(n - 2, 0)])
if fade > 0:
    c = tuple(int(v) for v in BG + fade * (np.array((120, 140, 170)) - BG)); paint(~pr, c, 2 if rowh >= 10 else 1)
if mode in ('primes', 'thin'):
    paint(pr & up, (56, 189, 248), dot); paint(pr & ~up, (245, 158, 11), dot)
elif mode == 'twins':
    paint(pr & ~twin, (40, 58, 82), dot); paint(twin & up, (56, 189, 248), dot); paint(twin & ~up, (245, 158, 11), dot)
E = np.array([k for k in EXP if k <= N]); m = np.isin(n, E); paint(m, (244, 114, 182), dot + 4)
raw = b''.join(b'\x00' + img[r].tobytes() for r in range(H))
def chunk(t, d): return struct.pack('>I', len(d)) + t + d + struct.pack('>I', zlib.crc32(t + d) & 0xffffffff)
open(out, 'wb').write(b'\x89PNG\r\n\x1a\n' + chunk(b'IHDR', struct.pack('>IIBBBBB', W, H, 8, 2, 0, 0, 0)) + chunk(b'IDAT', zlib.compress(raw, 6)) + chunk(b'IEND', b''))
lit = [c for c in range(w) if np.gcd(c + 1, w) == 1]
print(out, W, 'x', H, 'rows', rows, 'primes up', int((pr & up).sum()), 'down', int((pr & ~up).sum()), 'twins', int(twin.sum()), 'columns coprime to w', len(lit))
