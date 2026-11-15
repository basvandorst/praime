"""Independent recheck of the 5 pre-registered features: plain-Python feature code (no numpy, no PARI),
sympy primes, exact discovery Poisson-binomial, stratified-sample expectations on validation and final
ranges, Poisson tails, and the final-range binomial test conditional on the 13 exponents there."""
import json, math, os, sys
import sympy
from recheck_common import w, pb_tail, poisson_tail, stratified_E, EXPS, B

HERE = os.path.dirname(os.path.abspath(__file__))
pre = json.load(open(os.path.join(HERE, "prereg.json")))


def cf_maxq_norm(k):
    a0 = math.isqrt(k)
    m, d, a, mx = 0, 1, a0, 0
    while True:
        m = d * a - m
        d = (k - m * m) // d
        a = (a0 + m) // d
        if a == 2 * a0:
            return mx / math.sqrt(k)
        mx = max(mx, a)


def digits(k, b):
    out = []
    while k:
        out.append(k % b); k //= b
    return out


def b4_distinct_frac(k):
    d = digits(k, 4)
    return len(set(d)) / min(len(d), 4)


def b35_digitsum_z(k):
    d = digits(k, 35)
    n = len(d)
    return (sum(d) - n * 17) / math.sqrt(n * (35 * 35 - 1) / 12)


def angle_x2_2y2(k):
    if sympy.jacobi_symbol(-8 % k, k) != 1 or k <= 3:
        return None
    for y in range(1, math.isqrt(k // 2) + 1):
        r = k - 2 * y * y
        x = math.isqrt(r)
        if x * x == r:
            return math.atan2(math.sqrt(2) * y, x) / (math.pi / 2)
    return None


def five_is_7th(k):
    if k <= 11 or (k - 1) % 7:
        return None
    return 1.0 if pow(5, (k - 1) // 7, k) == 1 else 0.0


FUN = {"cfsqrt_maxq_norm": cf_maxq_norm, "b4_distinct_frac": b4_distinct_frac, "b35_digitsum_z": b35_digitsum_z,
       "angle_x2_2y2": angle_x2_2y2, "5_is_7th_power": five_is_7th}

# sanity of the brute-force pieces
assert cf_maxq_norm(13) == 1 / math.sqrt(13)          # sqrt13 = [3; 1,1,1,1,6]
assert angle_x2_2y2(11) is not None                     # 11 = 9 + 2
disc = list(sympy.primerange(11, 10 ** 5))
print(f"discovery primes {len(disc)}, E_D total {sum(map(w, disc)):.4f}")
for t in pre["top5"]:
    nm, (typ, th) = t["feature"], t["class"]
    f = FUN[nm]

    def member(k, f=f, typ=typ, th=th):
        v = f(k)
        if v is None:
            return False
        return v <= th if typ == 0 else (v > th if typ == 1 else v == th)

    inD = [p for p in disc if member(p)]
    ED = sum(map(w, inD))
    OD = sum(1 for e in EXPS if 11 <= e < 10 ** 5 and member(e))
    pD = pb_tail([w(p) for p in inD], OD)
    OV = sum(1 for e in EXPS if 10 ** 5 <= e < 10 ** 7 and member(e))
    EV, sV = stratified_E(10 ** 5, 10 ** 7, member, per_bin=2500, seed=11)
    OF = sum(1 for e in EXPS if 10 ** 7 <= e < B and member(e))
    EF, sF = stratified_E(10 ** 7, B, member, per_bin=2500, seed=12)
    EFtot, _ = stratified_E(10 ** 7, B, lambda k: True, per_bin=2500, seed=12)
    q = EF / EFtot
    pcond = sum(math.comb(13, j) * q ** j * (1 - q) ** (13 - j) for j in range(OF, 14))
    print("%-18s %-14s D %d/%.3f (p %.3f; harness %d/%.3f)  V %d/%.2f+-%.2f (Poisson p %.3f; harness %d/%.2f)  "
          "F %d/%.2f+-%.2f (Poisson p %.3f; given 13 in F, share %.3f -> binomial p %.3f; harness %s)"
          % (nm, t["class_text"], OD, ED, pD, t["D"][0], t["D"][1], OV, EV, sV, poisson_tail(EV, OV), t["V"][0], t["V"][1],
             OF, EF, sF, poisson_tail(EF, OF), q, pcond, "see final_exam.out"))
    sys.stdout.flush()
