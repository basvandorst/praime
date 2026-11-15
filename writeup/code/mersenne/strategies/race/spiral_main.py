"""Layout 1 (polyspiral): primes on the corner rays, Bateman-Horn predictions.

Ray (P, j, s): f(m) = P*m*(m-1)/2 + j*m + 1 + s  (= (P/2)m^2 + (j - P/2)m + 1 + s),
P = 3..40, j = 1..P, s = -6..6, m >= 1 with 2 <= f(m) <= N.

Bateman-Horn constant C(f) = prod_p (1 - w(p)/p)/(1 - 1/p), w(p)/p = share of m with p | f(m).
  p = 2: share counted directly over m mod 8 (f has period 4 mod 2 when P is odd).
  p odd, p | P: f = (2jm + 2(1+s))/2 mod p: w = 1 if p !| j; else w = p if p | 1+s else 0.
  p odd, p !| P: w = 1 + (D/p) with D = (2j - P)^2 - 8P(1+s) (discriminant of 2f).
Euler product over p < 10^4 (main), and p < 10^6 for the reported rays.
Output: spiral_main.json (every ray) and spiral_main.txt (tables).
"""
import json
import math
import numpy as np

N = 500_000
PB = 10_000
PB_BIG = 1_000_000


def sieve(n):
    s = np.ones(n + 1, bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return s


isp = sieve(max(N, PB_BIG))
odd_primes = np.nonzero(isp[:PB])[0][1:]          # 3 .. 9973
odd_primes_big = np.nonzero(isp[:PB_BIG])[0][1:]

rays = [(P, j, s) for P in range(3, 41) for j in range(1, P + 1) for s in range(-6, 7)]
RP = np.array([r[0] for r in rays], np.int64)
RJ = np.array([r[1] for r in rays], np.int64)
RS = np.array([r[2] for r in rays], np.int64)
D = (2 * RJ - RP) ** 2 - 8 * RP * (1 + RS)        # discriminant of 2 f


def f_val(P, j, s, m):
    return P * m * (m - 1) // 2 + j * m + 1 + s


def factor2(P, j, s):
    odd = sum(1 for m in range(8) if f_val(P, j, s, m) % 2 == 1)
    return 2.0 * odd / 8.0


def omega_vec(p, RPv, RJv, RSv, Dv, qr):
    """w(p) for all rays (vectorised), qr = Legendre table mod p."""
    w = 1 + qr[Dv % p].astype(np.int64)
    divP = (RPv % p == 0)
    if divP.any():
        jdiv = (RJv % p == 0)
        cdiv = ((1 + RSv) % p == 0)
        wl = np.where(~jdiv, 1, np.where(cdiv, p, 0))
        w = np.where(divP, wl, w)
    return w


def legendre_table(p):
    t = -np.ones(p, np.int8)
    t[0] = 0
    t[(np.arange(1, p, dtype=np.int64) ** 2) % p] = 1
    return t


def bh_constant(RPv, RJv, RSv, Dv, primes):
    C = np.array([factor2(int(a), int(b), int(c)) for a, b, c in zip(RPv, RJv, RSv)])
    for p in primes:
        p = int(p)
        w = omega_vec(p, RPv, RJv, RSv, Dv, legendre_table(p))
        C *= (1.0 - w / p) / (1.0 - 1.0 / p)
    return C


C = bh_constant(RP, RJ, RS, D, odd_primes)

count = np.zeros(len(rays), np.int64)
npts = np.zeros(len(rays), np.int64)
naive = np.zeros(len(rays))
var_cap = np.zeros(len(rays))
cnt3 = np.zeros(len(rays), np.int64); nai3 = np.zeros(len(rays))     # values f >= 10^3
cnt4 = np.zeros(len(rays), np.int64); nai4 = np.zeros(len(rays))     # values f >= 10^4
bern = np.zeros(len(rays))
red = np.array([d >= 0 and math.isqrt(int(d)) ** 2 == d for d in D])   # 2f has a square discriminant: f factors
for k, (P, j, s) in enumerate(rays):
    mmax = int(math.isqrt(2 * N // P)) + 3
    m = np.arange(1, mmax + 1, dtype=np.int64)
    f = f_val(P, j, s, m)
    f = f[(f >= 2) & (f <= N)]
    count[k] = int(isp[f].sum())
    npts[k] = len(f)
    inv = 1.0 / np.log(f.astype(float))
    naive[k] = inv.sum()
    q = np.minimum(C[k] * inv, 1.0)
    var_cap[k] = (q * (1 - q)).sum()
    b3, b4 = f >= 1000, f >= 10000
    cnt3[k] = int(isp[f[b3]].sum()); nai3[k] = inv[b3].sum()
    cnt4[k] = int(isp[f[b4]].sum()); nai4[k] = inv[b4].sum()
    bern[k] = math.sqrt(var_cap[k] / q.sum()) if q.sum() > 0 else 0.0

pred = C * naive
ok = C > 0
z = np.zeros(len(rays))
z[ok] = (count[ok] - pred[ok]) / np.sqrt(pred[ok])
z2 = np.zeros(len(rays))
ok2 = ok & (var_cap > 0)
z2[ok2] = (count[ok2] - pred[ok2]) / np.sqrt(var_cap[ok2])


def poly_str(P, j, s):
    c = 1 + s
    if P % 2 == 0:
        q = P // 2
        a, b = q, j - q
        t = f"{a}m^2"
        if b:
            t += f" {'+' if b > 0 else '-'} {abs(b) if abs(b) != 1 else ''}m"
        if c:
            t += f" {'+' if c > 0 else '-'} {abs(c)}"
        return t
    b2 = 2 * j - P
    t = f"({P}m^2"
    if b2:
        t += f" {'+' if b2 > 0 else '-'} {abs(b2) if abs(b2) != 1 else ''}m"
    if c:
        t += f" {'+' if 2 * c > 0 else '-'} {abs(2 * c)}"
    return t + ")/2"


def squarefree_kernel(n):
    if n == 0:
        return 0
    sgn = -1 if n < 0 else 1
    n = abs(n)
    k = 1
    p = 2
    while p * p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e % 2:
            k *= p
        p += 1
    return sgn * k * n


out_lines = []


def say(s=""):
    print(s)
    out_lines.append(s)


# Reference polynomials, same Euler product (brute-force roots for the reference ones)
def C_generic(fun, primes, period2=8):
    c = 2.0 * sum(1 for m in range(period2) if fun(m) % 2) / period2
    for p in primes:
        p = int(p)
        mm = np.arange(p, dtype=np.int64)
        w = int((fun(mm) % p == 0).sum())
        c *= (1 - w / p) / (1 - 1 / p)
    return c


say("Reference polynomials (p < 10^4, brute-force roots):")
say(f"  m^2 + m + 41 : C = {C_generic(lambda m: m*m + m + 41, odd_primes):.4f}")
say(f"  2m^2 - 199   : C = {C_generic(lambda m: 2*m*m - 199, odd_primes):.4f}")
say(f"  m^2 + 1      : C = {C_generic(lambda m: m*m + 1, odd_primes):.4f}")
say()

# Aggregate z statistics
say(f"rays: {len(rays)}; with C > 0: {int(ok.sum())}; C = 0: {int((~ok).sum())}")
say(f"total observed (C>0 rays) = {int(count[ok].sum())}, total predicted = {pred[ok].sum():.1f}, "
    f"ratio = {count[ok].sum() / pred[ok].sum():.4f}")
say(f"z = (obs-pred)/sqrt(pred): mean {z[ok].mean():+.4f}, sd {z[ok].std(ddof=1):.4f}, n {int(ok.sum())}")
say(f"z' = (obs-pred)/sqrt(sum q(1-q)), q = min(C/ln f, 1): mean {z2[ok2].mean():+.4f}, sd {z2[ok2].std(ddof=1):.4f}")
for lo, hi in [(0, 0.5), (0.5, 1), (1, 2), (2, 3), (3, 4), (4, 10)]:
    sel = ok & (C >= lo) & (C < hi)
    if sel.sum():
        say(f"  C in [{lo},{hi}): n {int(sel.sum()):5d}  obs/pred {count[sel].sum() / pred[sel].sum():.4f}  "
            f"z mean {z[sel].mean():+.3f} sd {z[sel].std(ddof=1):.3f}")
say(f"C = 0 rays: total primes on them {int(count[~ok].sum())} (small values only), max per ray {int(count[~ok].max())}")
say()
S_ = RS
say(f"reducible rays (disc of 2f a square): {int(red.sum())}, of which s = -1: {int((red & (S_ == -1)).sum())}; "
    f"primes on them {int(count[red].sum())}, max per ray {int(count[red].max())}; max C(<1e4) among them {C[red].max():.3f}")
oi = ok & ~red
say(f"IRREDUCIBLE rays with C > 0: {int(oi.sum())}")
for lab, cc, nn in [("all values", count, naive), ("values >= 10^3", cnt3, nai3), ("values >= 10^4", cnt4, nai4)]:
    pr = C * nn
    zz = (cc[oi] - pr[oi]) / np.sqrt(pr[oi])
    say(f"  {lab:15s}: obs {int(cc[oi].sum())}  pred {pr[oi].sum():.1f}  ratio {cc[oi].sum() / pr[oi].sum():.4f}  "
        f"z mean {zz.mean():+.4f}  sd {zz.std(ddof=1):.4f}")
say(f"  Bernoulli model (each value prime with prob min(C/ln f,1), independent): sd of z would be {bern[oi].mean():.3f} (mean over rays)")
for lo, hi in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 10)]:
    sel = oi & (C >= lo) & (C < hi)
    pr3 = C * nai3
    if sel.sum():
        say(f"  C in [{lo},{hi}): n {int(sel.sum()):5d}  obs/pred {count[sel].sum() / pred[sel].sum():.4f}  "
            f"(values >= 10^3: {cnt3[sel].sum() / pr3[sel].sum():.4f})  z mean {z[sel].mean():+.3f} sd {z[sel].std(ddof=1):.3f}")
say()
Crank = np.where(red, -1.0, C)     # C is meaningless for reducible rays: keep them out of the C ranking

# Per P table
sel_report = set()
say("Per P: brightest ray by C, and by observed count")
say(" P | by C: j  s  poly                       C     obs  pred   pts | by count: j  s  C     obs  pred  dens")
perP = {}
for P in range(3, 41):
    idx = np.nonzero(RP == P)[0]
    kc = idx[np.lexsort((-count[idx], -Crank[idx]))[0]]
    kn = idx[np.lexsort((-C[idx], -count[idx]))[0]]
    sel_report.update([int(kc), int(kn)])
    perP[P] = (int(kc), int(kn))
    say(f"{P:2d} | {RJ[kc]:2d} {RS[kc]:+d} {poly_str(P, RJ[kc], RS[kc]):26s} {C[kc]:5.3f} {count[kc]:4d} {pred[kc]:6.1f} {npts[kc]:4d} |"
        f" {RJ[kn]:2d} {RS[kn]:+d} {C[kn]:5.3f} {count[kn]:4d} {pred[kn]:6.1f} {count[kn] / npts[kn]:.3f}")
say()

top = np.lexsort((-C, -count))[:15]
sel_report.update(int(k) for k in top)
topC = np.argsort(-Crank, kind='stable')[:15]
sel_report.update(int(k) for k in topC)

# More accurate C for the reported rays: p < 10^6
idxs = np.array(sorted(sel_report))


def bh_constant_euler(P, j, s, Dk, primes):
    """Same product, Legendre symbol by Euler's criterion (scalar, for large p)."""
    c = factor2(P, j, s)
    for p in primes:
        p = int(p)
        if P % p == 0:
            w = 1 if j % p else (p if (1 + s) % p == 0 else 0)
        else:
            d = Dk % p
            w = 1 if d == 0 else (2 if pow(d, (p - 1) // 2, p) == 1 else 0)
        c *= (1.0 - w / p) / (1.0 - 1.0 / p)
    return c


Cbig_map = {int(k): bh_constant_euler(int(RP[k]), int(RJ[k]), int(RS[k]), int(D[k]), odd_primes_big) for k in idxs}

say("Top 15 rays by observed count (all P)")
say(" #  P  j  s  poly                         pts  obs   naive  C(<1e4) C(<1e6) pred   z     obs/naive")
for r, k in enumerate(top, 1):
    say(f"{r:2d} {RP[k]:2d} {RJ[k]:2d} {RS[k]:+d} {poly_str(RP[k], RJ[k], RS[k]):28s} {npts[k]:4d} {count[k]:4d} {naive[k]:7.2f} "
        f"{C[k]:7.4f} {Cbig_map[int(k)]:7.4f} {pred[k]:6.1f} {z[k]:+5.2f} {count[k] / naive[k]:.3f}")
say()
say("Top 15 rays by C (all P)")
say(" #  P  j  s  poly                         pts  obs   naive  C(<1e4) C(<1e6) pred   z     dens")
for r, k in enumerate(topC, 1):
    say(f"{r:2d} {RP[k]:2d} {RJ[k]:2d} {RS[k]:+d} {poly_str(RP[k], RJ[k], RS[k]):28s} {npts[k]:4d} {count[k]:4d} {naive[k]:7.2f} "
        f"{C[k]:7.4f} {Cbig_map[int(k)]:7.4f} {pred[k]:6.1f} {z[k]:+5.2f} {count[k] / npts[k]:.3f}")
say()

# Densest rays (primes per point) and per-P best density
dens = count / np.maximum(npts, 1)
topd = np.argsort(-dens)[:10]
say("Top 10 rays by density obs/pts")
for k in topd:
    say(f"  P {RP[k]:2d} j {RJ[k]:2d} s {RS[k]:+d} {poly_str(RP[k], RJ[k], RS[k]):28s} pts {npts[k]} obs {count[k]} dens {dens[k]:.3f} C {C[k]:.3f}")
say()

# Discriminant families
kern = np.array([squarefree_kernel(int(d)) for d in D])
say(f"rays with squarefree kernel of disc = -163 (Euler family): {int((kern == -163).sum())}")
say(f"rays with squarefree kernel of disc = 398 (2m^2-199 family): {int((kern == 398).sum())}")
say(f"disc range of 2f over all rays: {int(D.min())} .. {int(D.max())}; min |disc| of f*(2/gcd) etc. not needed")
say()

# per-P summary: max C, max count, sd of C (contrast), mean z
say(" P  maxC   maxObs  meanC  sdC   #C>3  mean z  pts(max)")
for P in range(3, 41):
    idx = np.nonzero(RP == P)[0]
    o = idx[C[idx] > 0]
    say(f"{P:2d} {C[idx].max():6.3f} {count[idx].max():6d} {C[idx].mean():6.3f} {C[idx].std():5.3f} {int((C[idx] > 3).sum()):4d} "
        f"{z[o].mean():+6.3f} {npts[idx].max():5d}")

json.dump({
    "N": N,
    "rays": [{"P": int(RP[k]), "j": int(RJ[k]), "s": int(RS[k]), "pts": int(npts[k]), "obs": int(count[k]),
              "naive": float(naive[k]), "C": float(C[k]), "pred": float(pred[k]), "z": float(z[k]),
              "Cbig": Cbig_map.get(k), "red": bool(red[k]), "obs3": int(cnt3[k]), "naive3": float(nai3[k]),
              "obs4": int(cnt4[k]), "naive4": float(nai4[k])} for k in range(len(rays))],
    "top_by_count": [int(k) for k in top],
    "top_by_C": [int(k) for k in topC],
    "perP": {str(P): v for P, v in perP.items()},
}, open("spiral_main.json", "w"))
open("spiral_main.txt", "w").write("\n".join(out_lines) + "\n")
