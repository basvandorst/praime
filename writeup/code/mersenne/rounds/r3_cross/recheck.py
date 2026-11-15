"""Independent recheck of the 5 pre-registered features on discovery and validation.
Own code throughout: sympy primes, mpmath / math feature values, plain-Python Collatz,
sympy Cornacchia, pure-Python Poisson-binomial DP, pure-Python deciles.
Shares nothing with harness.py/features.py except the frozen prereg.json (thresholds)."""
import json, math, sys, time
import mpmath, sympy
from sympy.solvers.diophantine.diophantine import cornacchia

T0 = time.time()
EXPS = {2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941,
        11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433, 1257787, 1398269,
        2976221, 3021377, 6972593}
G = 0.57721566490153286060651209
pre = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "prereg.json"))
disc = list(sympy.primerange(3, 10**5))
val = list(sympy.primerange(10**5, 10**7))
print(f"sympy primes: disc {len(disc)}, val {len(val)}")


def w(p):
    if p <= 7:
        return 1.0
    a = 6 if p % 4 == 1 else 2
    return min(1.0, math.exp(G) * math.log2(a * p) / p)


mpmath.mp.dps = 40
INV2PI = 1 / (2 * mpmath.pi)
FA = mpmath.mpf("2.502907875095892822283902873218")


def fracmp(p, a):
    x = p * a
    return float(x - mpmath.floor(x))


def glide(n):
    m, t = n, 0
    while True:
        m = 3 * m + 1 if m & 1 else m >> 1
        t += 1
        if m < n:
            return t


def ang2(p):
    if p % 8 not in (1, 3):
        return float("nan")
    sols = [s for s in cornacchia(1, 2, p) if s[0] > 0 and s[1] > 0]
    x, y = sols[0]
    assert x * x + 2 * y * y == p
    return math.atan2(math.sqrt(2) * y, x) / (math.pi / 2)


FEAT = {
    "frac(k*inv2pi)": lambda p: fracmp(p, INV2PI),
    "cos2pi(k*inv2pi)": lambda p: math.cos(p),             # cos(2 pi frac(k/2pi)) = cos k
    "frac(k*feigalpha)": lambda p: fracmp(p, FA),
    "Collatz glide": lambda p: float(glide(p)),
    "x^2+2y^2 angle": ang2,
}


def pb_tail(ws, o, L=60):
    d = [1.0] + [0.0] * L
    for q in ws:
        for n in range(L, 0, -1):
            d[n] = d[n] * (1 - q) + d[n - 1] * q
        d[0] *= (1 - q)
    return math.fsum(d[o:])


def quantile(xs, q):                      # numpy 'linear' rule, written out
    s = sorted(xs); h = (len(s) - 1) * q; lo = math.floor(h)
    return s[lo] + (h - lo) * (s[min(lo + 1, len(s) - 1)] - s[lo])


out = []
for f in pre:
    t1 = time.time()
    fn = FEAT[f["name"]]
    vd = [fn(p) for p in disc]
    # 1) thresholds: deciles of the defined discovery values
    xs = [v for v in vd if not math.isnan(v)]
    if f["kind"] == "cont":
        th = sorted(set(quantile(xs, k / 10) for k in range(1, 10)))
        th = [t for t in th if min(xs) <= t < max(xs)]
    else:
        th = [0.5]
    ok_thr = len(th) == len(f["thresholds"]) and all(abs(a - b) <= 1e-9 * max(1, abs(b)) for a, b in zip(th, f["thresholds"]))
    # 2) discovery: every candidate class, pure-Python Poisson-binomial; is the registered one the minimum?
    best = None
    for j, t in enumerate(th):
        for dr in ("<=", ">"):
            inc = [(not math.isnan(v)) and ((v <= t) if dr == "<=" else (v > t)) for v in vd]
            O = sum(1 for p, i in zip(disc, inc) if i and p in EXPS)
            pv = pb_tail([w(p) for p, i in zip(disc, inc) if i], O, 45)
            if best is None or pv < best[0] - 1e-15:
                best = (pv, j, dr, O, sum(w(p) for p, i in zip(disc, inc) if i))
    ok_choice = best[1] == f["chosen_j"] and best[2] == f["direction"]
    # 3) validation
    t = f["threshold"]; dr = f["direction"]
    Ov = 0; wv = []
    for p in val:
        v = fn(p)
        if math.isnan(v):
            continue
        if (v <= t) if dr == "<=" else (v > t):
            wv.append(w(p)); Ov += p in EXPS
    Ev = math.fsum(wv); pv = pb_tail(wv, Ov, 40)
    ln = (f"{f['name']}: thresholds match {ok_thr}; registered class is the discovery minimum {ok_choice} "
          f"(disc {best[3]}/{best[4]:.2f} p {best[0]:.3e} vs harness {f['disc_O']}/{f['disc_E']:.2f}); "
          f"validation {Ov}/{Ev:.3f} p {pv:.3e} vs harness {f['val_O']}/{f['val_E']:.3f} p {f['val_p']:.3e}  [{time.time()-t1:.0f}s]")
    print(ln, flush=True); out.append(ln)
# 4) final range [1e7, B): one streaming pass, Poisson tail with Le Cam bound (sum w^2) -- a different method
B = 83880077
FEX = {13466917, 20996011, 24036583, 25964951, 30402457, 32582657, 37156667, 42643801, 43112609, 57885161,
       74207281, 77232917, 82589933}
acc = {f["name"]: [0.0, 0.0, []] for f in pre}
Eall = 0.0; n = 0
fns = [(f["name"], FEAT[f["name"]], f["threshold"], f["direction"]) for f in pre]
for p in sympy.primerange(10**7, B):
    q = w(p); Eall += q; n += 1
    for nm, fn, t, dr in fns:
        v = fn(p)
        if v != v:
            continue
        if (v <= t) if dr == "<=" else (v > t):
            a = acc[nm]; a[0] += q; a[1] += q * q
            if p in FEX:
                a[2].append(p)
print(f"final range: {n} primes, E {Eall:.3f}")
out.append(f"final range [1e7,B): {n} primes, 13 exponents, E {Eall:.3f}")
for nm, (E, E2, ins) in acc.items():
    O = len(ins)
    pt = 1 - math.fsum(math.exp(-E) * E ** j / math.factorial(j) for j in range(O))
    pc = math.fsum(math.comb(13, j) * (E / Eall) ** j * (1 - E / Eall) ** (13 - j) for j in range(O, 14))
    ln = f"FINAL {nm}: {O}/{E:.3f} Poisson p {pt:.3e} (Le Cam error <= {E2:.1e}), binomial p given 13: {pc:.3f}; in class {ins}"
    print(ln, flush=True); out.append(ln)
open("recheck.txt", "w").write("\n".join(out) + f"\n({time.time()-T0:.0f}s total)\n")
