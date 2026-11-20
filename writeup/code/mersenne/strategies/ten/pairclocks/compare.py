"""Compare count_X.json with the frozen predictions of predict.py.

z = (observed share - predicted share) / sqrt(pred (1 - pred) / pairs)   (binomial scale)
Also: the single-prime classes (p mod 8, v2(p-1), p mod 5) against the single-prime model, and the
naive model that ignores the classes (both = (17/24)^2 or (2/3)^2).
Usage: python3 compare.py count_1000000000.json
"""
import json
import math
import sys
from fractions import Fraction as F
from predict import pair, single


def model_var(d, K=12):
    """Mean of P(1-P) over the classes of a pair: the variance per pair if, given the classes,
    each clock is an independent coin (the class mix itself is taken as exact)."""
    M = 1 << K
    cls5 = [c for c in range(1, 5) if (c + d) % 5 != 0]
    acc = {k: F(0) for k in ("H2", "H0", "L2", "L0")}
    n = 0
    for a in range(1, M, 2):
        for c in cls5:
            h1, l1 = single(a, c, K)
            h2, l2 = single((a + d) % M, (c + d) % 5, K)
            for k, p in (("H2", h1 * h2), ("H0", (1 - h1) * (1 - h2)),
                         ("L2", l1 * l2), ("L0", (1 - l1) * (1 - l2))):
                acc[k] += p * (1 - p)
            n += 1
    return {k: float(v / n) for k, v in acc.items()}

js = json.load(open(sys.argv[1]))
X = js["X"]
print("X = %d" % X)
zs = []
zero_ok = True
print(" d       pairs | quantity  predicted   observed    z(binomial)")
zm = []
for d in js["d"]:
    pr = pair(d, 12)
    mv = model_var(d)
    t = js["tally"][str(d)]
    n = t["N"][-1]
    for k in ("H2", "H0", "L2", "L0"):
        o = t[k][-1]
        p = float(pr[k])
        if pr[k] == 0:
            ok = o == 0
            zero_ok &= ok
            print("%2d %11d | %s  %-10s  %-10d  %s" % (d, n, k, "0 (proved)", o, "exact" if ok else "FAIL"))
            continue
        z = (o / n - p) / math.sqrt(p * (1 - p) / n)
        if mv[k] < 1e-15:
            print("%2d %11d | %s  %.6f    %.6f   %+6.2f      class-determined (pure equidistribution of pairs)"
                  % (d, n, k, p, o / n, z))
            continue
        zs.append((z, d, k))
        zm.append((d, k, p * (1 - p) - mv[k]))
        print("%2d %11d | %s  %.6f    %.6f   %+6.2f" % (d, n, k, p, o / n, z))
chi = sum(z * z for z, _, _ in zs)
print("\n%d nonzero coin predictions: sum z^2 = %.1f (mean %d if each z is standard; the z are correlated" % (
    len(zs), chi, len(zs)))
print("  (H2/H0 share coin pairs, a prime sits in pairs for several d), so no chi-square sd is quoted)")
print("max |z| = %.2f at d=%d %s" % max((abs(z), d, k) for z, d, k in zs))
print("all proved zeros exact:", zero_ok)
naive = {"H2": F(17, 24) ** 2, "L2": F(2, 3) ** 2}
print("\nnaive model (classes ignored): z for 'both'")
for d in js["d"]:
    t = js["tally"][str(d)]
    n = t["N"][-1]
    row = []
    for k in ("H2", "L2"):
        p = float(naive[k])
        row.append((t[k][-1] / n - p) / math.sqrt(p * (1 - p) / n))
    print("%2d  H %+9.1f   L %+9.1f" % (d, row[0], row[1]))

print("\nsingle primes by class (p mod 8, e = v2(p-1) capped at 11, p mod 5): model vs observed")
c = js["classes"]
zc = []
for key, n in enumerate(c["n"]):
    if n < 1000:
        continue
    m8, rest = divmod(key, 60)
    e, m5 = divmod(rest, 5)
    if e >= 11:
        continue
    q2 = m8 if e <= 2 else 1 + (1 << e)
    assert q2 % 8 == m8 and ((q2 - 1) & -(q2 - 1)).bit_length() - 1 == e, (m8, e)
    h, lu = single(q2, m5, 12)
    for name, p, o in (("H", h, c["h"][key]), ("L", lu, c["l"][key])):
        p = float(p)
        if p in (0.0, 1.0):
            if o != p * n:
                print("  class mod8=%d e=%d mod5=%d %s: predicted %d, observed %d of %d  FAIL" % (m8, e, m5, name, p, o, n))
            continue
        z = (o / n - p) / math.sqrt(p * (1 - p) / n)
        zc.append(z)
        if abs(z) > 3:
            print("  class mod8=%d e=%d mod5=%d %s: %.5f vs %.5f z=%+.1f" % (m8, e, m5, name, o / n, p, z))
print("  %d random classes: sum z^2 = %.1f, max |z| = %.2f; deterministic classes all exact unless FAIL above" % (
    len(zc), sum(z * z for z in zc), max(abs(z) for z in zc)))

if "EH2" in js["tally"][str(js["d"][0])]:
    print("\nsplit of each deviation: observed - predicted = (observed - E) + (E - predicted),")
    print("E = sum over the actual pairs of the model probability given their classes.")
    print("clock part z = (observed - E) / sqrt(sum P(1-P)): tests independence of the two clocks")
    print("class part = (E - predicted)/pairs: the mix of classes among the pairs (a prime race);")
    print("  its z uses the spread of the model chance over random classes, Var(P) = p(1-p) - mean P(1-P)")
    print("  (binomial variance = clock variance + class variance, exactly)")
    print(" d  q   clock z   class part (share)   class z")
    varP = {(d, k): v for d, k, v in zm}
    zk = []
    zc2 = []
    for d in js["d"]:
        pr = pair(d, 12)
        t = js["tally"][str(d)]
        n = t["N"][-1]
        for k in ("H2", "H0", "L2", "L0"):
            if pr[k] == 0:
                continue
            o, e, v = t[k][-1], t["E" + k][-1], t["V" + k][-1]
            cls = e / n - float(pr[k])
            vp = varP.get((d, k))
            if vp is None:
                vp = float(pr[k]) * (1 - float(pr[k]))
            zcl = cls / math.sqrt(vp / n)
            zk.append((zcl, d, k))
            if v < 1e-6:
                print("%2d %s   (no coin)  %+.2e          %+6.2f" % (d, k, cls, zcl))
                continue
            z = (o - e) / math.sqrt(v)
            zc2.append((z, d, k))
            print("%2d %s   %+6.2f    %+.2e          %+6.2f" % (d, k, z, cls, zcl))
    print("clock parts: %d values, sum z^2 = %.1f, max |z| = %.2f at d=%d %s" % (
        len(zc2), sum(z * z for z, _, _ in zc2), *max((abs(z), d, k) for z, d, k in zc2)))
    print("class parts: %d values, sum z^2 = %.1f, max |z| = %.2f at d=%d %s; largest |share| %.2e" % (
        len(zk), sum(z * z for z, _, _ in zk), *max((abs(z), d, k) for z, d, k in zk),
        max(abs(float(js["tally"][str(d)]["E" + k][-1]) / js["tally"][str(d)]["N"][-1] - float(pair(d, 12)[k]))
            for _, d, k in zk)))
