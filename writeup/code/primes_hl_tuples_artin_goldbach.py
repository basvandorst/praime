import sys, time
import numpy as np
from math import log, sqrt, pi, atan2, gcd
sys.path.insert(0, "/tmp/claude-1000/-projects-novovibe/f816ae1c-f982-49ee-8973-7dc30d982cf1/scratchpad/twins")
from explib import sieve, primes_upto, hl_constant, lam_array, C2, GAMMA

t0 = time.time()
X = 10**8
isp = sieve(2 * X + 10)
pr = np.flatnonzero(isp[:X + 1])
print("E-D8 Hardy-Littlewood counts to 1e8 for six tuple types with constants from the root counts:")
t = np.linspace(11, X, 400001)
for name, offs in (("twins", (0, 2)), ("cousins", (0, 4)), ("sexy", (0, 6)), ("triplets 0,2,6", (0, 2, 6)), ("triplets 0,4,6", (0, 4, 6)), ("quadruplets", (0, 2, 6, 8))):
    m = np.ones(len(pr), dtype=bool)
    for o in offs[1:]:
        m &= isp[pr + o]
    cnt = int(m.sum())
    c = hl_constant(offs)
    hl = c * np.trapezoid(1 / np.log(t) ** len(offs), t)
    print(f"  {name:15s}: {cnt:7d}  HL {hl:9.0f}  constant {c:.5f}  ratio {cnt/hl:.4f}")
sg = int(isp[2 * pr + 1].sum())
csg = hl_constant((0, 2))
hl_sg = csg * np.trapezoid(1 / (np.log(t) * np.log(2 * t + 1)), t)
print(f"  Sophie Germain : {sg:7d}  HL {hl_sg:9.0f}  constant {csg:.5f}  ratio {sg/hl_sg:.4f}")

print("E-D9 twin Skewes point: first x where the twin count exceeds 2C2 int_2^x dt/ln^2 t (literature: 1369391, Wolf):")
tw = pr[isp[pr + 2]]
xs = np.arange(10**5, 3 * 10**6, 1)
li2 = np.concatenate([[0.0], np.cumsum(1 / np.log(np.arange(2.0, 3 * 10**6 + 1)) ** 2)])
cnt = np.searchsorted(tw, xs, side='right')
diff = cnt - 2 * C2 * li2[xs - 1]
first = xs[np.argmax(diff > 0)] if (diff > 0).any() else None
print(f"  first crossing above 1e5: x = {first}; count there {cnt[np.argmax(diff > 0)]}, integral {2*C2*li2[first-1]:.1f}; fraction of x in [1e5, 3e6] with count above the integral: {(diff > 0).mean():.3f}")

print("E-D10 the twin constant from the data: 2C2 fitted as count / int dt/ln^2 t at 1e6, 1e7, 1e8:")
for x in (10**6, 10**7, 10**8):
    tt = np.linspace(2, x, 200001)
    print(f"  x={x:.0e}: {np.searchsorted(tw, x, side='right') / np.trapezoid(1/np.log(tt)**2, tt):.5f}  (2C2 = {2*C2:.5f})")

print("E-G3 Gaussian primes: angles of p = a^2 + b^2 (a > b > 0) for primes p = 1 mod 4 below 1e6, equidistribution of the angle atan(b/a) on (0, pi/4):")
p1 = pr[(pr % 4 == 1) & (pr < 10**6)]
angles = []
for p in p1.tolist():
    n = 2
    while pow(n, (p - 1) // 2, p) != p - 1:
        n += 1
    x_ = pow(n, (p - 1) // 4, p)
    a, b = p, x_
    while b * b > p:
        a, b = b, a % b
    a2 = int(sqrt(p - b * b))
    angles.append(atan2(min(a2, b), max(a2, b)))
angles = np.array(angles)
hist = np.histogram(angles, bins=8, range=(0, pi / 4))[0]
print(f"  {len(angles)} primes; octile counts {hist.tolist()} (uniform: {len(angles)/8:.0f} each); max |deviation|/sqrt(expected) = {np.abs(hist - len(angles)/8).max()/sqrt(len(angles)/8):.2f} (Hecke: equidistributed)  [{time.time()-t0:.0f} s]")

print("E-G4 Artin: fraction of primes below 1e6 with 2 a primitive root, against Artin's constant 0.3739558:")
spf = np.arange(10**6 + 1)
for q in primes_upto(1000):
    m = spf[q * q::q]
    spf[q * q::q] = np.where(m == np.arange(q * q, 10**6 + 1, q), q, m)
cnt = 0; tot = 0
for p in pr[(pr < 10**6) & (pr > 2)].tolist():
    m = p - 1; fac = set()
    while m > 1:
        q = spf[m]
        if q == m: fac.add(m); break
        fac.add(q); m //= q
    tot += 1
    if all(pow(2, (p - 1) // q, p) != 1 for q in fac):
        cnt += 1
print(f"  {cnt}/{tot} = {cnt/tot:.5f} (Artin 0.37396; proved under GRH by Hooley)  [{time.time()-t0:.0f} s]")

print("E-G5 Wieferich primes below 1e8 (2^(p-1) = 1 mod p^2):")
w = [p for p in pr[pr < 10**7].tolist() if pow(2, p - 1, p * p) == 1]
print(f"  below 1e7: {w} (known: 1093, 3511, none other below 1e17)")

print("E-D11 Goldbach comet classes: mean r(N) for N = 0 mod 3 against N not 0 mod 3, N in [1e6, 2e6] (singular series ratio (3-1)/(3-2) = 2):")
M = 2 * 10**6
Lam = lam_array(M)
F = np.fft.rfft(Lam, 2 * M + 2)
r = np.fft.irfft(F * F, 2 * M + 2)[:M + 1]
N = np.arange(10**6, 2 * 10**6, 2)
r0 = r[N[N % 3 == 0]].mean(); r1 = r[N[N % 3 != 0]].mean()
print(f"  mean r(N): N = 0 mod 3: {r0:.1f}, N not 0 mod 3: {r1:.1f}, ratio {r0/r1:.4f} (predicted 2); for 5 | N against 5 not dividing: ratio {r[N[N % 5 == 0]].mean()/r[N[N % 5 != 0]].mean():.4f} (predicted 4/3)")

print("E-C17 Hardy-Littlewood's second conjecture pi(x+y) <= pi(x) + pi(y) for y <= 3000 and x <= 1e8 (false far out, Hensley-Richards):")
pic = np.cumsum(isp[:X + 1])
worst = (0, 0, 0)
for y in (100, 300, 1000, 3000):
    excess = (pic[y:X + 1] - pic[:X + 1 - y]).max() - pic[y]
    if excess > worst[0]: worst = (excess, y, 0)
    print(f"  y={y}: max_x pi(x+y) - pi(x) = {(pic[y:X+1] - pic[:X+1-y]).max()}, pi(y) = {pic[y]}, margin {pic[y] - (pic[y:X+1] - pic[:X+1-y]).max()}")
print(f"[{time.time()-t0:.0f} s]")
