"""Independent recheck of spiral_main.py (layout 1, polyspiral).

Different methods on purpose:
  * the ray values are rebuilt from the geometry (T(m-1) + j*m + s, T(m) = 1 + P*m(m+1)/2)
    with plain Python integers, and tested with sympy.isprime;
  * w(p) is counted by brute force: the number of m in 0..p-1 with p | f(m) (no Legendre symbols),
    p = 2 by brute force over one full period m mod 4 (checked against m mod 8);
  * reducibility by sympy.factor_list, families by sympy.discriminant + squarefree part.
Then every per-ray number of spiral_main.json is compared, and the reported tables recomputed.
"""
import json
import math
import numpy as np
import sympy
from sympy import isprime, Rational, symbols, factor_list, discriminant, Poly
from sympy.ntheory.factor_ import core

N = 500_000
main = json.load(open("spiral_main.json"))
MR = main["rays"]


def T(P, m):
    return 1 + P * m * (m + 1) // 2


rays = [(P, j, s) for P in range(3, 41) for j in range(1, P + 1) for s in range(-6, 7)]
assert len(rays) == len(MR)
assert all((r["P"], r["j"], r["s"]) == x for r, x in zip(MR, rays))

# --- geometry check: ring m holds P*m numbers, corner j of ring m is T(m-1) + j*m ---
for P in (3, 4, 6, 40):
    for m in range(1, 50):
        assert T(P, m) - T(P, m - 1) == P * m
        for j in range(1, P + 1):
            c = T(P, m - 1) + j * m
            assert 2 * c == P * m * m + (2 * j - P) * m + 2
print("geometry: corner formula c_j(m) = (P/2)m^2 + (j-P/2)m + 1 verified (P=3,4,6,40, m<50)")

# --- counts with sympy.isprime ---
obs, pts, naive = [], [], []
o3l, n3l, o4l, n4l, logs = [], [], [], [], []
for P, j, s in rays:
    o = n = 0
    sm = 0.0
    o3 = o4 = 0
    s3 = s4 = 0.0
    lg = []
    m = 1
    while True:
        f = T(P, m - 1) + j * m + s
        if f > N:
            break
        if f >= 2:
            n += 1
            w = 1.0 / math.log(f)
            lg.append(w)
            sm += w
            pr = isprime(f)
            o += pr
            if f >= 1000:
                o3 += pr
                s3 += w
            if f >= 10000:
                o4 += pr
                s4 += w
        m += 1
    obs.append(o)
    pts.append(n)
    naive.append(sm)
    o3l.append(o3); n3l.append(s3); o4l.append(o4); n4l.append(s4); logs.append(lg)
obs = np.array(obs)
pts = np.array(pts)
naive = np.array(naive)

# --- C by brute-force root counting ---
plist = list(sympy.primerange(3, 10_000))
C = np.ones(len(rays))
index = {r: k for k, r in enumerate(rays)}
for P in range(3, 41):
    for p in plist:
        m = np.arange(p, dtype=np.int64)
        A = (P * (m * (m - 1) // 2)) % p
        for j in range(1, P + 1):
            cnt = np.bincount((A + j * m) % p, minlength=p)
            for s in range(-6, 7):
                w = int(cnt[(-1 - s) % p])
                C[index[(P, j, s)]] *= (1 - w / p) / (1 - 1 / p)
for k, (P, j, s) in enumerate(rays):
    odd4 = sum((T(P, m - 1) + j * m + s) % 2 for m in range(4))
    odd8 = sum((T(P, m - 1) + j * m + s) % 2 for m in range(8))
    assert 2 * odd4 == odd8
    C[k] *= 2 * odd4 / 4

pred = C * naive

# --- compare with main ---
mo = np.array([r["obs"] for r in MR])
mp = np.array([r["pts"] for r in MR])
mn = np.array([r["naive"] for r in MR])
mc = np.array([r["C"] for r in MR])
print("obs equal:", bool((mo == obs).all()), " pts equal:", bool((mp == pts).all()))
print("naive max rel diff: %.2e" % np.max(np.abs(mn - naive) / naive))
print("C max abs diff: %.2e" % np.max(np.abs(mc - C)))

# --- structure ---
x = symbols("x")
red = np.zeros(len(rays), bool)
kern = np.zeros(len(rays), np.int64)
for k, (P, j, s) in enumerate(rays):
    f = Rational(P, 2) * x ** 2 + (j - Rational(P, 2)) * x + 1 + s
    fl = factor_list(f)[1]
    red[k] = (len(fl) > 1) or any(e > 1 for _, e in fl) or Poly(fl[0][0], x).degree() < 2
    Dk = discriminant(2 * f, x)
    Dk = int(Dk)
    kern[k] = 0 if Dk == 0 else (1 if Dk > 0 else -1) * int(core(abs(Dk)))
S = np.array([r[2] for r in rays])
print("reducible rays:", int(red.sum()), " with s = -1:", int((red & (S == -1)).sum()), " primes on them:", int(obs[red].sum()),
      " max per ray:", int(obs[red].max()))
print("C = 0 rays:", int((C == 0).sum()), " primes on them:", int(obs[C == 0].sum()))
odd = np.array([r[0] % 2 == 1 for r in rays])
f2 = np.array([2 * sum((T(P, m - 1) + j * m + s) % 2 for m in range(4)) / 4 for P, j, s in rays])
print("odd P: 2-adic factor is 1 for every ray:", bool((f2[odd] == 1).all()), "; max C over odd P: %.4f" % C[odd].max(),
      "; max C over even P: %.4f" % C[~odd].max())
print("Euler family (core of disc = -163):", int((kern == -163).sum()), "; 2m^2-199 family (core 398):", int((kern == 398).sum()))

ok = (C > 0) & ~red
z = (obs[ok] - pred[ok]) / np.sqrt(pred[ok])
print("irreducible rays with C > 0: %d  obs %d  pred %.1f  ratio %.4f  z mean %+.4f sd %.4f"
      % (ok.sum(), obs[ok].sum(), pred[ok].sum(), obs[ok].sum() / pred[ok].sum(), z.mean(), z.std(ddof=1)))
for lab, oo, nn in [("values >= 10^3", o3l, n3l), ("values >= 10^4", o4l, n4l)]:
    oo = np.array(oo); pp = C * np.array(nn)
    zz = (oo[ok] - pp[ok]) / np.sqrt(pp[ok])
    print("  %s: obs %d pred %.1f ratio %.4f z mean %+.4f sd %.4f" % (lab, oo[ok].sum(), pp[ok].sum(), oo[ok].sum() / pp[ok].sum(),
                                                                    zz.mean(), zz.std(ddof=1)))
bsd = []
for k in np.nonzero(ok)[0]:
    q = np.minimum(C[k] * np.array(logs[k]), 1.0)
    bsd.append(math.sqrt((q * (1 - q)).sum() / q.sum()))
print("  Bernoulli-model sd of z, mean over rays: %.3f" % np.mean(bsd))
zfull = np.zeros(len(rays)); zfull[ok] = z
o3a = np.array(o3l); p3a = C * np.array(n3l)
for lo, hi in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 10)]:
    sel = ok & (C >= lo) & (C < hi)
    print("  C in [%d,%d): n %5d  obs/pred %.4f  (values >= 10^3: %.4f)  z mean %+.3f sd %.3f"
          % (lo, hi, sel.sum(), obs[sel].sum() / pred[sel].sum(), o3a[sel].sum() / p3a[sel].sum(), zfull[sel].mean(), zfull[sel].std(ddof=1)))
okall = C > 0
zall = (obs[okall] - pred[okall]) / np.sqrt(pred[okall])
print("all rays with C > 0 (incl. reducible): %d  ratio %.4f  z mean %+.4f sd %.4f"
      % (okall.sum(), obs[okall].sum() / pred[okall].sum(), zall.mean(), zall.std(ddof=1)))

# --- tables recomputed ---
order = sorted(range(len(rays)), key=lambda k: (-obs[k], -C[k]))[:15]
print("top 15 by count equal to main:", order == main["top_by_count"])
for k in order:
    print("  P %2d j %2d s %+d  obs %3d pts %3d C %.4f pred %.1f" % (*rays[k], obs[k], pts[k], C[k], pred[k]))
Crank = np.where(red, -1.0, C)
orderC = sorted(range(len(rays)), key=lambda k: -Crank[k])[:15]
print("top 15 by C equal to main (as sets):", set(orderC) == set(main["top_by_C"]))
for k in orderC:
    print("  P %2d j %2d s %+d  obs %3d pts %3d C %.4f pred %.1f dens %.3f" % (*rays[k], obs[k], pts[k], C[k], pred[k], obs[k] / pts[k]))
same = True
for P in range(3, 41):
    ks = [k for k in range(len(rays)) if rays[k][0] == P]
    kc = min(ks, key=lambda k: (-Crank[k], -obs[k]))
    kn = min(ks, key=lambda k: (-obs[k], -C[k]))
    same &= [kc, kn] == main["perP"][str(P)]
print("per-P brightest (by C, by count) equal to main:", same)

# centred polygonal numbers: ray j = P, s = 0
print("ray j=P, s=0 is the centred P-gonal numbers 1 + P m(m+1)/2:",
      all(T(P, m - 1) + P * m == 1 + P * m * (m + 1) // 2 for P in range(3, 41) for m in range(1, 30)))
k6 = index[(6, 6, 0)]
cub = [m for m in range(1, 1000) if (m + 1) ** 3 - m ** 3 <= N and isprime((m + 1) ** 3 - m ** 3)]
print("cuban primes (m+1)^3 - m^3 <= N:", len(cub), " ray (6,6,0) count:", int(obs[k6]))

# reference constants by brute force
def Cref(fun):
    c = 2 * sum(fun(m) % 2 for m in range(8)) / 8
    for p in plist:
        w = sum(1 for m in range(p) if fun(m) % p == 0)
        c *= (1 - w / p) / (1 - 1 / p)
    return c


print("C(m^2+m+41) = %.4f, C(2m^2-199) = %.4f (p < 10^4)" % (Cref(lambda m: m * m + m + 41), Cref(lambda m: 2 * m * m - 199)))
