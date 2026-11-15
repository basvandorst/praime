"""Round 6 feature map: quantities of the exponent n and the base b only.

features(kind, b, n) -> (dense float64 matrix, list of names, list of groups, list of types)
kind 'm': R = (b^n - 1)/(b - 1); kind 'p': R = (b^n + 1)/(b + 1).  The Mersenne case is ('m', 2).
Groups:
  res    one-hot n = r (mod m), 3 <= m <= 60, gcd(r, m) = 1, m != 2 mod 4 (for odd n those repeat m/2)
  dig    the 33 round-5 digit statistics (r5_features.digit_features)
  near   the six round-3 near-miss classes (frozen thresholds) and their continuous values
  bits   bits 1..10 of n, the 2nd and 3rd leading bits, popcount/runs/digit-sum z-scores
  size   log2 n and log2 log2 R
  char   Jacobi symbols (a | n) for a = b, -b, b-1, -(b-1), b+1, -(b+1)
  rmod   for primes l = 3..47: Legendre (R | l), [R = 1], [R = -1], [b^n = b] (mod l); R mod 16 one-hot
Continuous columns are listed with type 'cont'; every continuous column also gets five threshold
indicators (type 'bin') at quantiles fixed by the caller (thr argument) for the tree model.
"""
import math
import numpy as np
import sympy
import r5_features as FT

SMALLP = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]


def glide_vec(k):
    """Collatz glide (vectorised copy of round 5's mersenne.py; equals FT.collatz_glide)."""
    k = np.asarray(k, np.int64)
    v = k.copy(); g = np.zeros(len(k), np.int64); act = np.arange(len(k)); t = 0
    while len(act):
        t += 1
        x = v[act]
        x = np.where(x & 1, 3 * x + 1, x >> 1)
        v[act] = x
        d = x < k[act]
        g[act[d]] = t
        act = act[~d]
    return g


def cf_vec(k):
    """Largest partial quotient of sqrt(k) in the period / sqrt(k) (vectorised; equals FT.cf_maxq_norm)."""
    k = np.asarray(k, np.int64)
    a0 = np.floor(np.sqrt(k.astype(np.float64))).astype(np.int64)
    a0 = np.where((a0 + 1) ** 2 <= k, a0 + 1, a0); a0 = np.where(a0 ** 2 > k, a0 - 1, a0)
    m = np.zeros(len(k), np.int64); d = np.ones(len(k), np.int64); a = a0.copy(); mx = np.zeros(len(k), np.int64)
    act = np.arange(len(k))
    while len(act):
        mm, dd, aa, kk, A0 = m[act], d[act], a[act], k[act], a0[act]
        mm = dd * aa - mm; dd = (kk - mm * mm) // dd; aa = (A0 + mm) // dd
        m[act], d[act], a[act] = mm, dd, aa
        fin = aa == 2 * A0
        mx[act[~fin]] = np.maximum(mx[act[~fin]], aa[~fin])
        act = act[~fin]
    return mx / np.sqrt(k.astype(np.float64))


def jacobi_vec(a, n):
    return np.array([sympy.jacobi_symbol(a % int(x), int(x)) if int(x) > 2 else 0 for x in n], float)


def r_mod(kind, b, n, M):
    """R mod M for each n (Python ints, exact)."""
    g = b - 1 if kind == "m" else b + 1
    s = -1 if kind == "m" else 1
    out = np.empty(len(n), np.int64)
    for i, x in enumerate(n.tolist()):
        v = (pow(b, x, M * g) + s) % (M * g)
        assert v % g == 0
        out[i] = (v // g) % M
    return out


def features(kind, b, n):
    n = np.asarray(n, np.int64)
    cols, names, groups, types = [], [], [], []

    def add(v, name, grp, typ="bin"):
        cols.append(np.asarray(v, float)); names.append(name); groups.append(grp); types.append(typ)

    # res
    for m in range(3, 61):
        if m % 4 == 2:
            continue
        r_ = n % m
        for r in range(1, m):
            if math.gcd(r, m) == 1:
                add(r_ == r, f"n={r} mod {m}", "res")
    # dig
    for k, v in FT.digit_features(n).items():
        add(v, "dig " + k, "dig")
    # near
    for name, fn, test in FT.FROZEN:
        if name.startswith("H5"):
            v = glide_vec(n).astype(float)
        elif name.startswith("H6"):
            v = cf_vec(n)
        else:
            v = fn(n)
        cls = np.where(np.isnan(v), False, test(np.nan_to_num(v, nan=-9)))
        add(cls, name, "near")
        if name.startswith("H1"):
            continue
        cv = np.nan_to_num(v, nan=0.0)
        if name.startswith("H5"):
            cv = np.log1p(cv)
        add(cv, name.split()[0] + " value", "near", "cont")
    # bits
    blen = np.array([int(x).bit_length() for x in n])
    for i in range(1, 11):
        add((n >> i) & 1, f"bit {i}", "bits")
    add([(int(x) >> (L - 2)) & 1 for x, L in zip(n, blen)], "2nd leading bit", "bits")
    add([(int(x) >> (L - 3)) & 1 if L >= 3 else 0 for x, L in zip(n, blen)], "3rd leading bit", "bits")
    bstr = [bin(int(x))[2:] for x in n]
    pop = np.array([t.count("1") for t in bstr], float)
    runs = np.array([1 + sum(t[i] != t[i + 1] for i in range(len(t) - 1)) for t in bstr], float)
    s10 = [str(int(x)) for x in n]
    dsum = np.array([sum(map(int, t)) for t in s10], float)
    dlen = np.array([len(t) for t in s10], float)
    add((pop - 1 - (blen - 1) / 2) / np.sqrt(np.maximum(blen - 1, 1) / 4), "popcount z", "bits", "cont")
    add((runs - 1 - (blen - 1) / 2) / np.sqrt(np.maximum(blen - 1, 1) / 4), "runs z", "bits", "cont")
    add((dsum - 4.5 * dlen) / np.sqrt(8.25 * dlen), "digit sum z", "bits", "cont")
    add(np.log2(n) - np.floor(np.log2(n)), "frac log2 n", "bits", "cont")
    # size
    add(np.log2(n), "log2 n", "size", "cont")
    g = b - 1 if kind == "m" else b + 1
    log2R = n * math.log2(b) - math.log2(g)
    add(np.log2(np.maximum(log2R, 1.0)), "log2 log2 R", "size", "cont")
    # char
    for a, nm in ((b, "b"), (-b, "-b"), (b - 1, "b-1"), (-(b - 1), "-(b-1)"), (b + 1, "b+1"), (-(b + 1), "-(b+1)")):
        add(jacobi_vec(a, n), f"({nm}|n)", "char", "cont")
    # rmod
    for l in SMALLP:
        R = r_mod(kind, b, n, l)
        bn = np.array([pow(b, int(x), l) for x in n])
        add([sympy.legendre_symbol(int(v), l) if v else 0 for v in R], f"(R|{l})", "rmod", "cont")
        add(R == 1, f"R=1 mod {l}", "rmod")
        add(R == l - 1, f"R=-1 mod {l}", "rmod")
        add(bn == b % l, f"b^n=b mod {l}", "rmod")
    R16 = r_mod(kind, b, n, 16)
    for r in range(1, 16, 2):
        add(R16 == r, f"R={r} mod 16", "rmod")
    X = np.column_stack(cols)
    return X, names, groups, types


def add_thresholds(X, names, groups, types, thr=None):
    """Append 5 threshold indicators x <= t for every continuous column; thr: dict name -> list."""
    newthr = {} if thr is None else thr
    extra, en, eg, et = [], [], [], []
    for j, (nm, ty) in enumerate(zip(names, types)):
        if ty != "cont":
            continue
        if thr is None:
            ts = sorted(set(np.quantile(X[:, j], [0.1, 0.25, 0.5, 0.75, 0.9]).tolist()))
            newthr[nm] = ts
        for t in newthr[nm]:
            extra.append(X[:, j] <= t); en.append(f"{nm}<={t:.4g}"); eg.append(groups[j]); et.append("bin")
    X2 = np.column_stack([X] + [np.asarray(e, float) for e in extra]) if extra else X
    return X2, names + en, groups + eg, types + et, newthr
