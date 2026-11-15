"""Independent recheck of race_main.py (layout 2, chebrace).

Different methods on purpose:
  * primes from a pure-Python bytearray sieve, counted against pi(10^7) = 664579;
  * squares mod P by the local criterion (CRT): a is a square mod P iff it is a square mod every prime power
    p^k || P (odd p: Legendre symbol (a/p) = 1 by Euler's criterion; 2: a = 1 mod 2, 4, 8 for 2, 4, 2^k>=8);
  * h walked prime by prime in plain Python;
  * the average over ln x is an exact integral of the step function h against d(ln x)/Li(sqrt x),
    with Li from mpmath.li(offset=True) and a cumulative Simpson rule on a log grid (no sampling grid);
  * the prime-power identity behind the prediction checked exactly with Fractions.
Compares with race_main.json.
"""
import json
import math
from fractions import Fraction
import numpy as np
import mpmath as mp
from scipy.integrate import cumulative_simpson

XMAX = 10 ** 7
main = json.load(open("race_main.json"))

# --- primes ---
sv = bytearray([1]) * (XMAX + 1)
sv[0] = sv[1] = 0
for i in range(2, math.isqrt(XMAX) + 1):
    if sv[i]:
        sv[i * i::i] = bytes(len(range(i * i, XMAX + 1, i)))
primes = [i for i in range(XMAX + 1) if sv[i]]
assert len(primes) == 664579
print("pi(10^7) = 664579 ok; pi(5*10^5) =", sum(1 for p in primes if p <= 500000))


def factorize(n):
    f = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def is_square_local(a, P):
    for p, k in factorize(P).items():
        if p == 2:
            m = 2 if k == 1 else (4 if k == 2 else 8)
            if a % m != 1:
                return False
        elif pow(a % p, (p - 1) // 2, p) != 1:
            return False
    return True


# --- grid for G(x) = int_{10^3}^{x} d(ln t) / Li(sqrt t) ---
mp.mp.dps = 20
lo, hi = math.log(1e3), math.log(XMAX)
K = 40001
u = np.linspace(lo, hi, K)
inv = np.array([1.0 / float(mp.li(mp.e ** (uu / 2), offset=True)) for uu in u])
G = np.concatenate([[0.0], cumulative_simpson(inv, x=u)])


def Gx(x):
    return np.interp(np.log(x), u, G)


def avg_exact(hs, X, nN):
    """(1/(ln X - ln 10^3)) * int_{10^3}^{X} h(x) / ((1/2) nN Li(sqrt x)) d ln x, h piecewise constant."""
    ps = np.array(primes, dtype=float)
    hv = np.array(hs, dtype=float)
    a = np.clip(ps, 1e3, X)
    b = np.clip(np.append(ps[1:], np.inf), 1e3, X)
    w = Gx(b) - Gx(a)
    w[b <= a] = 0.0
    return float((hv * w).sum() / (0.5 * nN) / (math.log(X) - math.log(1e3)))


ok_all = True
rows = {}
for P in range(3, 41):
    units = [a for a in range(1, P) if math.gcd(a, P) == 1]
    R = [a for a in units if is_square_local(a, P)]
    Nn = [a for a in units if a not in R]
    r = sum(1 for b in units if b * b % P == 1)
    assert len(R) * r == len(units)
    inR = [False] * P
    inN = [False] * P
    for a in R:
        inR[a] = True
    for a in Nn:
        inN[a] = True
    nR, nN = len(R), len(Nn)
    h = 0
    hs = []
    stats = {}
    for X in (500_000, XMAX):
        stats[X] = {"neg": 0, "first": None, "last": None, "first1000": None}
    for p in primes:
        c = p % P
        if inN[c]:
            h += nR
        elif inR[c]:
            h -= nN
        hs.append(h)
        if h < 0:
            for X in (500_000, XMAX):
                if p <= X:
                    st = stats[X]
                    st["neg"] += 1
                    st["last"] = p
                    if st["first"] is None:
                        st["first"] = p
                    if st["first1000"] is None and p > 1000:
                        st["first1000"] = p
    m = main[str(P)]
    ok = (m["R"], m["N"], m["r"]) == (nR, nN, r)
    row = {"R": nR, "N": nN, "r": r}
    for X in (500_000, XMAX):
        mm = m[str(X)]
        nX = sum(1 for p in primes if p <= X)
        avg = avg_exact(hs, X, nN)
        st = stats[X]
        ok &= mm["n_primes_neg"] == st["neg"] and mm["first_neg"] == st["first"] and mm["last_neg"] == st["last"]
        ok &= mm["first_neg_after_1000"] == st["first1000"]
        ok &= abs(mm["share_primes_neg"] - st["neg"] / nX) < 1e-12
        ok &= mm["h_end"] == hs[nX - 1]
        ok &= abs(mm["avg_ratio"] - avg) < 5e-3
        row[X] = (avg, mm["avg_ratio"], st["neg"], st["first"], st["last"])
    rows[P] = row
    ok_all &= ok
    print(f"P {P:2d} |R| {nR:2d} |N| {nN:2d} r {r}  avg(10^7) exact {row[XMAX][0]:.4f} grid {row[XMAX][1]:.4f}  "
          f"avg(5e5) exact {row[500000][0]:.4f} grid {row[500000][1]:.4f}  #neg {row[XMAX][2]:6d} first {row[XMAX][3]} "
          f"last {row[XMAX][4]}  {'OK' if ok else 'MISMATCH'}")
print("all per-P numbers agree with race_main.json:", ok_all)
for X in (XMAX, 500_000):
    nv = [P for P in rows if rows[P][X][3] is None]
    print(f"P with h >= 0 at every prime <= {X}: {nv} ({len(nv)} moduli)")
print("max |exact - grid| average (10^7):", max(abs(rows[P][XMAX][0] - rows[P][XMAX][1]) for P in rows))
print("max |exact - grid| average (5e5):", max(abs(rows[P][500000][0] - rows[P][500000][1]) for P in rows))

# --- the prime-power identity: p^k is a square mod P iff k even or p is a square ---
bad = 0
for P in range(3, 41):
    for p in primes:
        if p * p > XMAX:
            break
        if P % p == 0:
            continue
        sq_p = is_square_local(p % P, P)
        k, q = 2, p * p
        while q <= XMAX:
            if is_square_local(q % P, P) != (k % 2 == 0 or sq_p):
                bad += 1
            k += 1
            q *= p
print("p^k (k >= 2, p^k <= 10^7, p !| P) square mod P exactly when k even or p square: violations =", bad)

# exact identity  h(x) = H_Pi(x) + |N| sum_{k even} pi*(x^(1/k))/k - sum_{k odd >= 3} h(x^(1/k))/k
def iroot(x, k):
    y = int(round(x ** (1.0 / k)))
    while y ** k > x:
        y -= 1
    while (y + 1) ** k <= x:
        y += 1
    return y


def h_at(P, x, nR, nN, inR, inN):
    t = 0
    for p in primes:
        if p > x:
            break
        c = p % P
        t += nR if inN[c] else (-nN if inR[c] else 0)
    return t


idbad = 0
for P in (3, 4, 8, 15, 19, 24, 40):
    units = [a for a in range(1, P) if math.gcd(a, P) == 1]
    R = [a for a in units if is_square_local(a, P)]
    inR = [a in R for a in range(P)]
    inN = [(a in units) and (a not in R) for a in range(P)]
    nR, nN = sum(inR), sum(inN)
    for x in (10 ** 3, 12345, 10 ** 5, 999_983):
        PiN = PiR = Fraction(0)
        for p in primes:
            if p > x:
                break
            if P % p == 0:
                continue
            k, q = 1, p
            while q <= x:
                if inR[q % P]:
                    PiR += Fraction(1, k)
                else:
                    PiN += Fraction(1, k)
                k += 1
                q *= p
        HPi = nR * PiN - nN * PiR
        rhs = HPi
        k = 2
        while 2 ** k <= x:
            y = iroot(x, k)
            if k % 2 == 0:
                rhs += Fraction(nN, k) * sum(1 for p in primes if p <= y and P % p)
            else:
                rhs -= Fraction(h_at(P, y, nR, nN, inR, inN), k)
            k += 1
        if rhs != h_at(P, x, nR, nN, inR, inN):
            idbad += 1
print("exact identity h = H_Pi + |N| sum_even pi*(x^(1/k))/k - sum_odd h(x^(1/k))/k: failures =", idbad,
      "(P in 3,4,8,15,19,24,40; x in 10^3, 12345, 10^5, 999983)")

# --- P = 4 and P = 3 classic values ---
c1 = c3 = 0
first4 = None
for p in primes:
    if p % 4 == 1:
        c1 += 1
    elif p % 4 == 3:
        c3 += 1
    if c1 > c3 and first4 is None:
        first4 = p
print("mod 4: first p with pi(p;4,1) > pi(p;4,3):", first4)
c1 = c2 = 0
first3 = None
for p in primes:
    if p % 3 == 1:
        c1 += 1
    elif p % 3 == 2:
        c2 += 1
    if c1 > c2 and first3 is None:
        first3 = p
print("mod 3: first p <= 10^7 with pi(p;3,1) > pi(p;3,2):", first3, " ; final pi(10^7;3,2) - pi(10^7;3,1) =", c2 - c1)

# --- exact log-scale integrals: refined ratio h/((1/2)|N| Pi*(sqrt x)) and the log density of {h < 0} ---
print()
parr = np.array(primes, dtype=np.int64)
grp = {True: [], False: []}
grpL = {True: [], False: []}
okx = True
for P in range(3, 41):
    units = [a for a in range(1, P) if math.gcd(a, P) == 1]
    R = [a for a in units if is_square_local(a, P)]
    inR = np.zeros(P, bool); inR[R] = True
    unit = np.zeros(P, bool); unit[units] = True
    nR, nN = len(R), len(units) - len(R)
    r = len(units) // nR
    hs = np.cumsum(np.where(inR[parr % P], -nN, np.where(unit[parr % P], nR, 0)))
    ppw = []
    for p in primes:
        if p * p > XMAX:
            break
        if P % p == 0:
            continue
        k, q = 1, p
        while q * q <= XMAX:
            ppw.append((q * q, 1.0 / k))      # Pi*(sqrt x) jumps by 1/k at x = (p^k)^2
            k += 1
            q *= p
    ppw.sort()
    jx = np.array([a for a, _ in ppw], dtype=np.int64)
    jw = np.cumsum([w for _, w in ppw])
    for X in (XMAX, 500_000):
        B = np.unique(np.concatenate([[1000], parr[(parr > 1000) & (parr < X)], jx[(jx > 1000) & (jx < X)], [X]]))
        a, b = B[:-1], B[1:]
        hseg = hs[np.searchsorted(parr, a, side="right") - 1]
        piseg = jw[np.searchsorted(jx, a, side="right") - 1]
        dl = np.log(b.astype(float)) - np.log(a.astype(float))
        L = math.log(X) - math.log(1000)
        refined = float((hseg / (0.5 * nN * piseg) * dl).sum() / L)
        logneg = float(dl[hseg < 0].sum() / L)
        mm = main[str(P)][str(X)]
        okx &= abs(refined - mm["avg_ratio_Pi"]) < 5e-3 and abs(logneg - mm["logdens_neg"]) < 2e-3
        if X == XMAX:
            has_even = (P - 1) in R or r >= 4         # a nonprincipal real character with chi(-1) = +1 exists
            grp[has_even].append(refined)
            grpL[has_even].append(main[str(P)][str(X)]["avg_ratio"])
            print(f"P {P:2d} refined avg exact {refined:.4f} (grid {mm['avg_ratio_Pi']:.4f})  log density h<0 exact {logneg:.5f} "
                  f"(grid {mm['logdens_neg']:.5f})  even real character: {has_even}")
print("refined averages and log densities agree with race_main (tol 5e-3, 2e-3):", okx)
for g in (False, True):
    print(f"P {'with' if g else 'without'} an even real character: n {len(grp[g])}, mean h/((1/2)|N|Li) {np.mean(grpL[g]):.4f}, "
          f"mean h/((1/2)|N|Pi*) {np.mean(grp[g]):.4f}")

# --- odd P and 2P run the same race except for the prime 2 ---
def hseq(P):
    units = [a for a in range(1, P) if math.gcd(a, P) == 1]
    R = [a for a in units if is_square_local(a, P)]
    inR = np.zeros(P, bool); inR[R] = True
    unit = np.zeros(P, bool); unit[units] = True
    return np.cumsum(np.where(inR[parr % P], -(len(units) - len(R)), np.where(unit[parr % P], len(R), 0))), 2 % P in R
for P in range(3, 21, 2):
    a, two_sq = hseq(P)
    b, _ = hseq(2 * P)
    d = a - b
    print(f"P {P:2d} vs {2*P:2d}: h_P - h_2P constant = {int(d[0])} for all primes <= 10^7: {bool((d == d[0]).all())}; 2 is a square mod P: {two_sq}")
