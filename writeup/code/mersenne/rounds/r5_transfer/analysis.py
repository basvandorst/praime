#!/usr/bin/env python3
"""Pooled tests of the round-3 near-misses (and a broad feature family) on generalized
repunit / Wagstaff families, with discovery and holdout halves of the families.

Universe per family: odd prime n <= N_b, log2 R >= 64, n not dividing b -/+ 1.
Nulls (each with one constant c_f per family, fitted so that sum w = observed count):
  G  growth only:        w = c_f ln(n)/n
  C  growth + classes:   w = c_f A_f(n mod 24) ln(2n)/n, A_f(class) = sum s(n)/sum ln(2n)/ln Q_n over
                         the class, s(n) = 1 if the sieve found no factor q = 2kn+1 <= Q_n = 2Kn+1
                         (derived from the factors, not from the hits)
  S  sieve-conditioned:  w = c_f s(n) ln(Q_n)/ln(R)   (theory: c_f = e^gamma)
Statistic per feature: O = hits in class, E = sum w in class, pooled over the families of a half;
p from the exact Poisson-binomial of the w in class (one-sided, upper and lower).
"""
import json, math, os, sys
import numpy as np
import features as FT

HERE = os.path.dirname(os.path.abspath(__file__))
EG = math.exp(0.5772156649015329)
SEED = 20261115


def load():
    fams = []
    for line in open(os.path.join(HERE, "data.jsonl")):
        d = json.loads(line)
        b, kind, K = d["b"], d["kind"], d["K"]
        g = b - 1 if kind == "m" else b + 1
        rows = [r for r in d["recs"] if r[1] != "ndiv"]
        n = np.array([r[0] for r in rows], np.int64)
        log2R = n * math.log2(b) - math.log2(g)
        keep = log2R >= 64
        st = [r[1] for r, k in zip(rows, keep) if k]
        n = n[keep]; log2R = log2R[keep]
        assert "mr" not in st
        hit = np.array([s == "PRP" for s in st], float)
        surv = np.array([s in ("PRP", "fermat") for s in st], float)
        Q = 2.0 * K * n + 1
        fams.append(dict(name=f"{kind}{b}", kind=kind, b=b, n=n, hit=hit, surv=surv, Q=Q,
                         lnR=log2R * math.log(2)))
    return fams


def nulls(f):
    n = f["n"].astype(float)
    out = {}
    wG = np.log(n) / n
    cls = f["n"] % 24
    base = np.log(2 * n) / np.log(f["Q"])
    A = np.zeros(len(n))
    for c in np.unique(cls):
        m = cls == c
        A[m] = f["surv"][m].sum() / base[m].sum()
    wC = A * np.log(2 * n) / n
    wS = f["surv"] * np.log(f["Q"]) / f["lnR"]
    O = f["hit"].sum()
    for nm, w in (("G", wG), ("C", wC), ("S", wS)):
        c = O / w.sum() if w.sum() > 0 else 0.0
        out[nm] = np.minimum(c * w, 0.95)
        if nm == "S":
            f["cS"] = c
    f["thS"] = (EG * wS).sum()
    return out


def pb_tails(w, O):
    """P(X >= O), P(X <= O) for X = sum of Bernoulli(w)."""
    w = w[w > 0]
    L = int(O) + 2
    pmf = np.zeros(L); pmf[0] = 1.0          # pmf[j] for j <= O, pmf[L-1] = mass above O
    for p in w:
        top = pmf[L - 1] + pmf[L - 2] * p
        pmf[1:L - 1] = pmf[1:L - 1] * (1 - p) + pmf[0:L - 2] * p
        pmf[0] *= (1 - p)
        pmf[L - 1] = top
    le = pmf[:L - 1].sum()          # P(X <= O)
    ge = 1 - le + pmf[L - 2]        # P(X >= O)
    return min(1.0, max(ge, 0.0)), min(1.0, max(le, 0.0))


def split(fams):
    rng = np.random.default_rng(SEED)
    disc, hold = [], []
    for kind in ("m", "p"):
        names = sorted([f["name"] for f in fams if f["kind"] == kind])
        rng.shuffle(names)
        h = len(names) // 2
        disc += names[:h]; hold += names[h:]
    return set(disc), set(hold)


def pooled(fams, half, featname, NULL):
    O = E = 0.0; ws = []
    for f in fams:
        if f["name"] not in half:
            continue
        m = f["feat"][featname]
        O += f["hit"][m].sum(); E += f["w"][NULL][m].sum(); ws.append(f["w"][NULL][m])
    w = np.concatenate(ws) if ws else np.zeros(0)
    pu, pl = pb_tails(w, O)
    return O, E, pu, pl


def main():
    fams = load()
    allfeat_names = None
    for f in fams:
        f["w"] = nulls(f)
        n = f["n"]
        feat = {}
        for name, fn, test in FT.FROZEN:
            v = fn(n)
            feat[name] = np.where(np.isnan(v), False, test(np.nan_to_num(v, nan=-9)))
        # positive control (theorem): q = 2n+1 prime and +-b a QR mod q  ->  q | R
        import sympy
        sgn = 1 if f["kind"] == "m" else -1
        g = f["b"] - sgn
        pc = np.array([sympy.isprime(2 * int(x) + 1) and g % (2 * int(x) + 1) != 0 and
                       sympy.jacobi_symbol(sgn * f["b"] % (2 * int(x) + 1), 2 * int(x) + 1) == 1
                       for x in n])
        feat["P0 control: 2n+1 prime, (+-b|2n+1)=1"] = pc
        for m in range(3, 61):
            for r in range(1, m):
                if math.gcd(r, m) == 1:
                    feat[f"n={r} mod {m}"] = n % m == r
        for k, v in FT.digit_features(n).items():
            feat["dig " + k] = v
        f["feat"] = feat
        allfeat_names = list(feat)
    disc, hold = split(fams)
    out = {"families": [], "disc": sorted(disc), "hold": sorted(hold), "tests": {}}
    tot = {"O": 0, "n": 0}
    for f in fams:
        out["families"].append(dict(name=f["name"], n=int(len(f["n"])), O=int(f["hit"].sum()),
                                    surv=int(f["surv"].sum()), cS=f["cS"], thS=f["thS"],
                                    half="D" if f["name"] in disc else "H"))
    for name in allfeat_names:
        res = {}
        for NULL in ("G", "C", "S"):
            for hn, half in (("D", disc), ("H", hold), ("A", disc | hold)):
                res[f"{NULL}{hn}"] = pooled(fams, half, name, NULL)
        out["tests"][name] = res
    with open(os.path.join(HERE, "results.json"), "w") as fh:
        json.dump(out, fh)
    print("saved", len(allfeat_names), "features")


if __name__ == "__main__":
    main()
