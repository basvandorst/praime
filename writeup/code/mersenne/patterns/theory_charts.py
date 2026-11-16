"""SVG charts for the theories section of prime2-patterns.html, from the outputs in ../theory/."""
import re, math, sys, os, numpy as np
TH = '../theories'
OUT = 'theory'
os.makedirs(OUT, exist_ok=True)
GRID, AX, MUTED, PINK, BLUE, YEL = '#1f2e44', '#5b6b82', '#9fb0c6', '#f472b6', '#38bdf8', '#fde68a'
def frame(w, h, label):
    return [f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{label}" style="width:100%;height:auto;display:block">']
def txt(x, y, s, anchor='start', fill=MUTED, size=12):
    return f'<text x="{x:.1f}" y="{y:.1f}" fill="{fill}" font-size="{size}" text-anchor="{anchor}">{s}</text>'
rows = []
for line in open(f'{TH}/goldrace12.txt'):
    m = re.match(r'\s+([\d,]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s*$', line)
    if m: rows.append((float(m.group(1).replace(',', '')), float(m.group(2)), float(m.group(3)), float(m.group(4))))
rows = np.array(rows)
zrows = []
for line in open(f'{TH}/goldzeros.txt'):
    m = re.match(r'\s+([\d,]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s*$', line)
    if m: zrows.append(float(m.group(3)))
zrows = np.array(zrows)
W, H, L, R, T, B = 1000, 340, 60, 16, 16, 44
x0, x1 = math.log(rows[0, 0]), math.log(4e6)
ymin, ymax = -130, 10
X = lambda v: L + (math.log(v) - x0) / (x1 - x0) * (W - L - R)
Y = lambda v: T + (ymax - v) / (ymax - ymin) * (H - T - B)
s = frame(W, H, "Goldbach's comet, n = 0 mod 12: class-1 deficit against the prime race")
for v in range(-120, 11, 20):
    s.append(f'<line x1="{L}" x2="{W - R}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" stroke="{AX if v == 0 else GRID}" stroke-width="1"/>' + txt(L - 6, Y(v) + 4, str(v), 'end'))
for v in (2.5e5, 5e5, 1e6, 2e6, 4e6):
    s.append(f'<line x1="{X(v):.1f}" x2="{X(v):.1f}" y1="{T}" y2="{H - B}" stroke="{GRID}" stroke-width="1"/>' + txt(X(v), H - B + 18, f'{v:,.0f}', 'end' if v == 4e6 else 'middle'))
mid = np.sqrt(rows[:, 0] * np.append(rows[1:, 0], 4e6))
s.append('<polyline fill="none" stroke="' + MUTED + '" stroke-width="2" stroke-dasharray="6 5" points="' + ' '.join(f'{X(a):.1f},{Y(c):.1f}' for a, c in zip(mid, rows[:, 3])) + '"/>')
s.append('<polyline fill="none" stroke="' + BLUE + '" stroke-width="2" points="' + ' '.join(f'{X(a):.1f},{Y(c):.1f}' for a, c in zip(mid, rows[:, 2])) + '"/>')
s.append('<polyline fill="none" stroke="' + YEL + '" stroke-width="1.5" stroke-dasharray="2 3" points="' + ' '.join(f'{X(a):.1f},{Y(c):.1f}' for a, c in zip(mid, zrows)) + '"/>')
for a, d in zip(mid, rows[:, 1]):
    s.append(f'<circle cx="{X(a):.1f}" cy="{Y(d):.1f}" r="4" fill="{PINK}"><title>n near {a:,.0f}: mean gap {d:.1f}</title></circle>')
lx, ly = L + 12, Y(-122) - 56
s.append(txt(L + 8, H - 6, 'n (log scale)') + txt(lx, ly, '● measured: mean of (1,11)+(11,1) − (5,7)−(7,5) over a bin', 'start', PINK)
         + txt(lx, ly + 16, '— predicted from the prime race at n (no parameter)', 'start', BLUE) + txt(lx, ly + 32, '··· predicted from 359 zeros of L(s, χ₁₂), no prime count', 'start', YEL) + txt(lx, ly + 48, '- - prime squares alone: −S(n)', 'start', MUTED))
s.append('</svg>')
open(f'{OUT}/goldrace.svg', 'w').write('\n'.join(s))
zr = []
for line in open(f'{TH}/zrace.txt'):
    m = re.match(r'\s*(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$', line)
    if m: zr.append([float(v) for v in m.groups()])
zr = np.array(zr)
zeros = {}
for line in open(f'{TH}/zrace.txt'):
    m = re.match(r'P\s+(\d+) .*lowest zero ([\d.]+)', line)
    if m: zeros[int(m.group(1))] = float(m.group(2))
W, H, L, R, T, B = 1000, 320, 56, 16, 16, 58
n = len(zr); bw = (W - L - R) / n
ymax = 0.25
Y = lambda v: T + (ymax - v) / ymax * (H - T - B)
s = frame(W, H, 'Share of time the residues lead, by modulus: measured against the zero prediction')
for v in (0, 0.05, 0.10, 0.15, 0.20, 0.25):
    s.append(f'<line x1="{L}" x2="{W - R}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" stroke="{AX if v == 0 else GRID}" stroke-width="1"/>' + txt(L - 6, Y(v) + 4, f'{v:.2f}', 'end'))
for i, r in enumerate(zr):
    P = int(r[0]); cx = L + (i + 0.5) * bw
    s.append(f'<rect x="{cx - 16:.1f}" y="{Y(r[1]):.1f}" width="14" height="{Y(0) - Y(r[1]):.1f}" rx="2" fill="{PINK}"><title>P = {P}: measured {r[1]:.3f}</title></rect>')
    s.append(f'<rect x="{cx + 2:.1f}" y="{Y(r[5]):.1f}" width="14" height="{Y(0) - Y(r[5]):.1f}" rx="2" fill="{BLUE}"><title>P = {P}: predicted {r[5]:.3f} (lowest zero {zeros[P]:.4f})</title></rect>')
    s.append(txt(cx, H - B + 18, str(P), 'middle', '#d6e2f2') + txt(cx, H - B + 34, f'{zeros[P]:.2f}', 'middle', YEL, 11))
s.append(txt(L - 6, H - B + 18, 'P', 'end') + txt(L - 6, H - B + 34, 'zero', 'end', YEL, 11))
s.append(txt(W - R, T + 14, '■ measured, odd primes to 10⁸ (log-time from 10³)', 'end', PINK) + txt(W - R, T + 30, '■ predicted from L(s, χ_P) alone', 'end', BLUE))
s.append('</svg>')
open(f'{OUT}/zrace.svg', 'w').write('\n'.join(s))
sys.path.insert(0, TH)
from hlbytes import expected
O8 = [-1, 1, 5, 7, 11, 13, 17, 19]
a = np.load(f'{TH}/s1e9.npz')
E, _ = expected(O8, int(a['J']), 9)
obs = a['cnt'].astype(float)
W, H, L, R, T, B = 1000, 420, 64, 16, 16, 44
lo, hi = 0, 7.2
X = lambda v: L + (math.log10(max(v, 1)) - lo) / (hi - lo) * (W - L - R)
Y = lambda v: T + (hi - math.log10(max(v, 1))) / (hi - lo) * (H - T - B)
s = frame(W, H, 'The 256 wheel-24 bytes up to a billion: count against the Hardy-Littlewood prediction')
for e in range(0, 8):
    v = 10 ** e
    s.append(f'<line x1="{L}" x2="{W - R}" y1="{Y(v):.1f}" y2="{Y(v):.1f}" stroke="{GRID}" stroke-width="1"/>' + txt(L - 6, Y(v) + 4, f'10{"⁰¹²³⁴⁵⁶⁷"[e]}', 'end'))
    s.append(f'<line x1="{X(v):.1f}" x2="{X(v):.1f}" y1="{T}" y2="{H - B}" stroke="{GRID}" stroke-width="1"/>' + txt(X(v), H - B + 18, f'10{"⁰¹²³⁴⁵⁶⁷"[e]}', 'middle'))
s.append(f'<line x1="{X(1):.1f}" y1="{Y(1):.1f}" x2="{X(10 ** hi):.1f}" y2="{Y(10 ** hi):.1f}" stroke="{AX}" stroke-width="1"/>')
for v in range(256):
    if E[v] <= 0: continue
    z = (obs[v] - E[v]) / math.sqrt(E[v])
    s.append(f'<circle cx="{X(E[v]):.1f}" cy="{Y(obs[v]):.1f}" r="4" fill="{YEL if v in (223, 251) else BLUE}"><title>byte {v} = {v:08b}: {obs[v]:,.0f} found, {E[v]:,.1f} predicted (z {z:+.1f})</title></circle>')
s.append(txt(L + 8, H - 6, 'predicted count (log)') + txt(W - R, H - B - 10, '27 impossible bytes: predicted 0, found 0 (not drawn)', 'end')
         + txt(L + 12, T + 14, 'found (log)', 'start') + txt(X(13.22) + 10, Y(14) + 4, '223 and 251: the septuplets, 13.2 predicted, 14 and 14 found', 'start', YEL))
s.append('</svg>')
open(f'{OUT}/bytes.svg', 'w').write('\n'.join(s))
print('ok')
