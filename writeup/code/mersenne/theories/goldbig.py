"""Theory 9b/9c on a long range, without building the comet. For a bin [lo, hi) of n and a class pair (a, b) mod 12, the sum of
#{(p, q): p = a, q = b, p + q = n} over the n in the bin is F(hi) - F(lo), F(X) = sum over primes p = a, p < X of #{q = b, q <= X - 1 - p}
(exact integers, by binary search in the sorted prime lists). The square term: the sum of S(n) over the n = c in a sub-bin is
sum over primes s >= 5 of #{q = c - 1 mod 12: q + s^2 in the sub-bin}, also exact. The race weight SIG(n)/Q(n) and the zero weight
change slowly and are taken at the middle of each of 100 sub-bins per bin.
usage: goldbig.py LO HI NBINS [zeros]  (bins equally spaced in log n over [LO, HI))"""
import sys, time, math, numpy as np
LO, HI, NB = int(float(sys.argv[1])), int(float(sys.argv[2])), int(sys.argv[3])
ZE = len(sys.argv) > 4
t0 = time.time()
X = HI + 2
r = int(X ** 0.5) + 2
b0 = np.ones(r + 1, bool); b0[:2] = False
for i in range(2, int(r ** 0.5) + 1):
    if b0[i]: b0[i * i::i] = False
base = np.nonzero(b0)[0]; base = base[base > 2]
parts = []
SEG = 1 << 24
for lo in range(1, X, 2 * SEG):
    hi = min(lo + 2 * SEG, X)
    a = np.ones((hi - lo + 1) // 2, bool)
    for p in base:
        if p * p >= hi: break
        st = max(p * p, ((lo + p - 1) // p) * p)
        if st % 2 == 0: st += p
        a[(st - lo) // 2::p] = False
    v = lo + 2 * np.nonzero(a)[0]
    parts.append(v[(v >= 5) & (v < hi)].astype(np.uint32))
P = np.concatenate(parts); del parts
P = P[P != 1]
C = {c: P[P % 12 == c] for c in (1, 5, 7, 11)}
del P
print(f'primes 5 .. {X:,}: {sum(len(v) for v in C.values()):,} ({time.time() - t0:.0f} s)', flush=True)
sq = np.array([int(s) * int(s) for s in base if s >= 5 and s * s <= X], np.int64)
def cnt(c, y):
    y = np.asarray(y, np.int64)
    out = np.searchsorted(C[c], np.clip(y, 0, 2 ** 32 - 1).astype(np.uint32), side='right')
    return np.where(y < 0, 0, out)
def F(Xe, a, b):
    pa = C[a][C[a] < Xe - 5].astype(np.int64)
    return int(cnt(b, Xe - 1 - pa).sum())
edges = np.unique(np.round(np.exp(np.linspace(math.log(LO), math.log(HI), NB + 1))).astype(np.int64))
NBIN = len(edges) - 1
TYPES = {0: [((1, 11), 2), ((5, 7), -2)], 2: [((1, 1), 1), ((7, 7), -1)], 6: [((1, 5), 2), ((7, 11), -2)], 10: [((5, 5), 1), ((11, 11), -1)]}
Fv = {}
for c, T in TYPES.items():
    for (a, b), _ in T:
        if (a, b) not in Fv:
            Fv[(a, b)] = np.array([F(int(e), a, b) for e in edges], np.int64)
print(f'pair counts at {len(edges)} edges ({time.time() - t0:.0f} s)', flush=True)
def ncount(lo, hi, c): return (hi - 1 - c) // 12 - (lo - 1 - c) // 12
SUB = 100
sub = np.unique(np.concatenate([np.round(np.linspace(edges[i], edges[i + 1], SUB + 1)).astype(np.int64) for i in range(NBIN)]))
mid = (sub[:-1] + sub[1:]) / 2.0
which = np.searchsorted(edges, sub[:-1], side='right') - 1
pim = {c: cnt(c, np.floor(mid)) for c in C}
Qm = np.searchsorted(sq, np.floor(mid), side='right').astype(float)
if ZE:
    Z4 = np.loadtxt('zeros_m4.txt'); Z12 = np.loadtxt('zeros_12.txt')
    def Zsum(x, g):
        l = np.log(x)[:, None]; u = (0.5 + 1j * g)[None, :] * l
        return 2 * (np.exp(u) / u * (1 + 1 / u + 2 / u ** 2)).real.sum(1)
    z4 = Zsum(mid, Z4); z12 = Zsum(mid, Z12)
def Ssum(c):
    b, shift = {0: (11, 0), 2: (1, 0), 6: (5, 0), 10: (7, 2)}[c]
    tot = np.zeros(len(sub), np.int64)
    for i in range(0, len(sub), 400):
        e = sub[i:i + 400]
        y = e[:, None] - 1 - sq[None, :] - shift
        tot[i:i + 400] = cnt(b, y).sum(1)
    return np.diff(tot)
print(f'n in [{LO:,}, {HI:,}), {NBIN} bins of log n, {len(mid)} sub-bins; prediction from the race at n and (if given) from the zeros only')
print(' n mod 12   corr race   corr zeros   corr squares   |err| race   |err| zeros   |err| squares   mean Delta   mean -S')
rows = {}
for c, T in TYPES.items():
    tot = sum(sg * Fv[ab] for ab, sg in T)
    nb = np.array([ncount(edges[i], edges[i + 1], c) for i in range(NBIN)])
    md = np.diff(tot) / nb
    Ss = Ssum(c)
    sig = {0: 2 * (pim[1] + pim[11] - pim[5] - pim[7]), 2: 2 * (pim[1] - pim[7]), 6: 2 * (pim[1] + pim[5] - pim[7] - pim[11]), 10: 2 * (pim[5] - pim[11])}[c]
    pr = np.bincount(which, weights=Ss * sig / Qm, minlength=NBIN) / nb
    ms = -np.bincount(which, weights=Ss, minlength=NBIN) / nb if c != 10 else np.zeros(NBIN)
    if ZE:
        zw = {0: -(1 + 2 * z12 / Qm), 2: -(1 + (z4 + z12 - 0.5) / Qm), 6: -(1 + (2 * z4 - 1) / Qm), 10: (z12 - z4 + 0.5) / Qm}[c]
        pz = np.bincount(which, weights=Ss * zw, minlength=NBIN) / nb
    else:
        pz = np.full(NBIN, np.nan)
    cs = np.corrcoef(md, ms)[0, 1] if c != 10 else float('nan')
    cz = np.corrcoef(md, pz)[0, 1] if ZE else float('nan')
    print(f'{c:>9} {np.corrcoef(md, pr)[0, 1]:>11.3f} {cz:>12.3f} {cs:>14.3f} {np.mean(np.abs(md - pr)):>12.2f} {np.mean(np.abs(md - pz)):>13.2f} {np.mean(np.abs(md - ms)):>15.2f} {md.mean():>12.2f} {ms.mean():>9.2f}', flush=True)
    rows[c] = (md, pr, pz, ms)
for c in (0,):
    md, pr, pz, ms = rows[c]
    print(f'\nper bin (n = {c} mod 12): lower edge, mean Delta, from the race, from the zeros, squares only')
    for i in range(NBIN):
        print(f'  {edges[i]:>14,} {md[i]:>10.2f} {pr[i]:>10.2f} {pz[i]:>10.2f} {ms[i]:>10.2f}')
for c in (2, 6, 10):
    md, pr, pz, ms = rows[c]
    print(f'\nper bin (n = {c} mod 12): mean Delta / race / zeros: ' + '; '.join(f'{x:.1f}/{y:.1f}/{z:.1f}' for x, y, z in zip(md, pr, pz)))
print(f'\n({time.time() - t0:.0f} s)')
