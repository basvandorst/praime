#!/usr/bin/env python3
"""Independent recheck of round 6.

A. Features: every column recomputed from its name with plain Python integers (own code, no feats.py,
   no round-5 features.py), for all 10,925 training rows and for the 52 Mersenne exponents plus 1,500
   population rows; standardisation recomputed; compared with train.npz / mersenne.npz.
B. The S null from data.jsonl in plain Python; the real out-of-family gain of every configuration
   recomputed from the stored out-of-fold scores with own gain code; the L2 configuration refitted fold by
   fold with Newton's method (exact Hessian, a different solver from L-BFGS) and its gain recomputed.
C. Mersenne: the final L2 model refitted on all families by Newton, the exponents scored on the
   independently computed features, the baseline, percentiles and the gain in bits recomputed;
   trial-factoring survival rechecked with Python pow on 300 population rows.
Output: recheck.txt."""
import json, math, os, sys, cmath
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "r5_transfer", "data.jsonl")
EG = math.exp(0.5772156649015329)
out = []
P = lambda s: (print(s, flush=True), out.append(s))


def isprime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):       # deterministic below 3.3e24
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def jacobi(a, n):
    a %= n; r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            r = -r
        a %= n
    return r if n == 1 else 0


def glide(n):
    v, t = n, 0
    while True:
        v = v // 2 if v % 2 == 0 else 3 * v + 1
        t += 1
        if v < n:
            return t


def cfmax(n):
    a0 = math.isqrt(n); m, d, a, mx = 0, 1, a0, 0
    while True:
        m = d * a - m; d = (n - m * m) // d; a = (a0 + m) // d
        if a == 2 * a0:
            return mx / math.sqrt(n)
        mx = max(mx, a)


def angle(n):
    if n % 8 not in (1, 3):
        return 0.0
    y = 1
    while 2 * y * y < n:
        x2 = n - 2 * y * y; x = math.isqrt(x2)
        if x * x == x2:
            return math.atan2(math.sqrt(2) * y, x) / (math.pi / 2)
        y += 1
    return 0.0


def Rmod(kind, b, n, M):
    g = b - 1 if kind == "m" else b + 1
    num = pow(b, n, M * g) + (-1 if kind == "m" else 1)
    return (num % (M * g)) // g % M


def raw_row(kind, b, n, names, thr):
    s10, s2 = str(n), bin(n)[2:]
    L2, L10 = len(s2), len(s10)
    pop, runs = s2.count("1"), 1 + sum(s2[i] != s2[i + 1] for i in range(L2 - 1))
    dsum = sum(map(int, s10))
    dz = (dsum - 4.5 * L10) / math.sqrt(8.25 * L10)
    pz = (pop - 1 - (L2 - 1) / 2) / math.sqrt(max(L2 - 1, 1) / 4)
    rz = (runs - 1 - (L2 - 1) / 2) / math.sqrt(max(L2 - 1, 1) / 4)
    H = {"H2": angle(n), "H3": (n / (2 * math.pi)) % 1.0, "H4": math.cos(n), "H5": math.log1p(glide(n)),
         "H6": cfmax(n)}
    g = b - 1 if kind == "m" else b + 1
    cont = {"H2 value": H["H2"], "H3 value": H["H3"], "H4 value": H["H4"], "H5 value": H["H5"],
            "H6 value": H["H6"], "popcount z": pz, "runs z": rz, "digit sum z": dz,
            "frac log2 n": math.log2(n) - math.floor(math.log2(n)), "log2 n": math.log2(n),
            "log2 log2 R": math.log2(max(n * math.log2(b) - math.log2(g), 1.0))}
    for a, nm in ((b, "b"), (-b, "-b"), (b - 1, "b-1"), (-(b - 1), "-(b-1)"), (b + 1, "b+1"), (-(b + 1), "-(b+1)")):
        cont[f"({nm}|n)"] = jacobi(a, n) if n > 2 else 0
    for l in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
        r = Rmod(kind, b, n, l)
        cont[f"(R|{l})"] = 0 if r == 0 else (1 if pow(r, (l - 1) // 2, l) == 1 else -1)
    dig = {"dsum_z<=-1": dz <= -1, "dsum_z>=1": dz >= 1, "popcount_z<=-1": pz <= -1, "popcount_z>=1": pz >= 1,
           "popcount odd": pop % 2 == 1, "binary runs_z<=-1": rz <= -1, "binary runs_z>=1": rz >= 1,
           "distinct dec digits = len": len(set(s10)) == L10, "has dec digit 0": "0" in s10,
           "all dec digits odd": all(c in "13579" for c in s10), "dec palindrome": s10 == s10[::-1],
           "bin palindrome": s2 == s2[::-1], "dec reversal prime": isprime(int(s10[::-1])),
           "bin reversal prime": isprime(int(s2[::-1], 2)), "binary length odd": L2 % 2 == 1,
           "dec length odd": L10 % 2 == 1, "leading two bits 10": s2[:2] == "10", "trailing bits 1111": n % 16 == 15,
           "n+2 prime": isprime(n + 2), "n-2 prime": isprime(n - 2), "2n+1 prime": isprime(2 * n + 1),
           "2n-1 prime": isprime(2 * n - 1)}
    for dd in range(1, 10):
        dig[f"leading dec digit {dd}"] = s10[0] == str(dd)
    near = {"H1 n=1 mod 7": n % 7 == 1, "H2 x^2+2y^2 angle >0.75 (n=1,3 mod 8)": H["H2"] > 0.75,
            "H3 frac(n/2pi)<=0.103": H["H3"] <= 0.103, "H4 cos n>0.803": H["H4"] > 0.803,
            "H5 Collatz glide>21": glide(n) > 21, "H6 cf sqrt maxq/sqrt n<=0.397": H["H6"] <= 0.397}
    row = []
    for nm in names:
        if nm in cont:
            row.append(cont[nm]); continue
        if nm in near:
            row.append(float(near[nm])); continue
        if nm.startswith("dig "):
            row.append(float(dig[nm[4:]])); continue
        if nm.startswith("n=") and " mod " in nm:
            r, m = nm[2:].split(" mod "); row.append(float(n % int(m) == int(r))); continue
        if nm.startswith("bit "):
            row.append(float((n >> int(nm[4:])) & 1)); continue
        if nm == "2nd leading bit":
            row.append(float(s2[1] == "1")); continue
        if nm == "3rd leading bit":
            row.append(float(len(s2) >= 3 and s2[2] == "1")); continue
        if nm.startswith("R=") and " mod " in nm:
            r, m = nm[2:].split(" mod "); m = int(m); r = int(r) % m
            row.append(float(Rmod(kind, b, n, m) == r)); continue
        if nm.startswith("b^n=b mod "):
            l = int(nm.split()[-1]); row.append(float(pow(b, n, l) == b % l)); continue
        base, t = nm.rsplit("<=", 1)
        tv = [x for x in thr[base] if f"{x:.4g}" == t]
        assert len(tv) == 1, nm
        row.append(float(cont[base] <= tv[0]))
    return row


def newton_l2(X1, y, off, lam, iters=30):
    """Penalised logistic regression, offset, intercept (column 0) unpenalised; Newton steps."""
    p = X1.shape[1]
    pen = np.full(p, lam); pen[0] = 0.0
    th = np.zeros(p)
    for it in range(iters):
        z = off + X1 @ th
        pr = 1 / (1 + np.exp(-z))
        g = X1.T @ (pr - y) + pen * th
        Hm = (X1 * (pr * (1 - pr))[:, None]).T @ X1 + np.diag(pen)
        step = np.linalg.solve(Hm, g)
        th -= step
        if np.max(np.abs(step)) < 1e-10:
            break
    return th


def gain(y, p, w):
    return sum((math.log2(pi / wi) if yi else math.log2((1 - pi) / (1 - wi))) for yi, pi, wi in zip(y, p, w))


def main():
    d = np.load(os.path.join(HERE, "train.npz"))
    names = [str(x) for x in d["names"]]; thr = json.loads(str(d["thr"])); types = d["types"]
    # ---- rows and S null from data.jsonl, plain Python
    rows = []
    for fi, line in enumerate(open(DATA)):
        f = json.loads(line)
        b, kind, K = f["b"], f["kind"], f["K"]
        g = b - 1 if kind == "m" else b + 1
        for n, st, _ in f["recs"]:
            if st == "ndiv":
                continue
            l2R = n * math.log2(b) - math.log2(g)
            if l2R >= 64 and st in ("PRP", "fermat"):
                rows.append([fi, kind, b, n, int(st == "PRP"), math.log(2 * K * n + 1) / (l2R * math.log(2))])
    O, Sw = {}, {}
    for r in rows:
        O[r[0]] = O.get(r[0], 0) + r[4]; Sw[r[0]] = Sw.get(r[0], 0.0) + r[5]
    w = [min(O[r[0]] / Sw[r[0]] * r[5], 0.95) for r in rows]
    P(f"rows {len(rows)}, hits {sum(r[4] for r in rows)}; equal to train.npz: n {np.array_equal([r[3] for r in rows], d['n'])}, "
      f"y {np.array_equal([r[4] for r in rows], d['y'])}")
    # ---- A. features
    Xraw = np.array([raw_row(r[1], r[2], r[3], names, thr) for r in rows])
    mu = np.where(types == "cont", Xraw.mean(0), 0.0); sd = np.where(types == "cont", Xraw.std(0), 1.0)
    Xr = (Xraw - mu) / sd; del Xraw
    dX = np.abs(Xr - d["X"].astype(np.float64))
    P(f"A. training features: {Xr.shape}, max |difference| {dX.max():.2e} (float32 storage), "
      f"columns with a difference > 1e-4: {int(np.sum(dX.max(0) > 1e-4))}")
    me = np.load(os.path.join(HERE, "mersenne.npz"))
    k = me["k"]; st = me["stratum"]; W = me["W"]; surv = me["surv"]
    rng = np.random.default_rng(1)
    chk = np.concatenate([np.flatnonzero(st < 0), rng.choice(np.flatnonzero(st >= 15), 1500, replace=False)])
    XmR = (np.array([raw_row("m", 2, int(k[i]), names, thr) for i in chk]) - mu) / sd
    dM = np.abs(XmR - me["X"][chk].astype(np.float64))
    P(f"A. Mersenne features ({len(chk)} rows incl. the 52 exponents): max |difference| {dM.max():.2e}")
    # ---- B. the real CV gains
    yv = np.array([r[4] for r in rows], float); wv = np.array(w); fam = np.array([r[0] for r in rows])
    act = wv > 0
    cv = json.load(open(os.path.join(HERE, "cv_results.json")))
    oof = np.load(os.path.join(HERE, "oof_real.npz"))
    worst = 0.0
    for key in oof.files:
        fs, c = key.split("|")
        z = np.log(wv[act] / (1 - wv[act])) + oof[key][act]
        pm = 1 / (1 + np.exp(-z))
        gg = gain(yv[act], pm, wv[act])
        worst = max(worst, abs(gg - cv["real"][fs][c][0]))
    P(f"B. gains of all {len(oof.files)} configurations recomputed from the stored out-of-fold scores: "
      f"largest difference {worst:.2e} bits")
    lam = float(json.loads(str(np.load(os.path.join(HERE, "final_models.npz"))["chosen"]))["L2"].split("=")[1])
    fold_of = {int(a): v for a, v in cv["folds"].items()}
    fo = np.array([fold_of[f] for f in fam])
    X1 = np.column_stack([np.ones(len(rows)), Xr])
    off = np.where(act, np.log(np.maximum(wv, 1e-300) / (1 - wv)), 0.0)
    zo = np.zeros(len(rows))
    for kf in range(5):
        tr = act & (fo != kf); te = act & (fo == kf)
        th = newton_l2(X1[tr], yv[tr], off[tr], lam)
        zo[te] = X1[te] @ th
    pm = 1 / (1 + np.exp(-(off[act] + zo[act])))
    gN = gain(yv[act], pm, wv[act])
    P(f"B. L2 lam={lam:g} refitted by Newton, fold by fold: out-of-family gain {gN:+.4f} bits "
      f"(run_cv: {cv['real']['all'][f'L2 lam={lam:g}'][0]:+.4f})")
    # ---- C. Mersenne
    th = newton_l2(X1[act], yv[act], off[act], lam)
    fm = np.load(os.path.join(HERE, "final_models.npz"))
    P(f"C. final L2 (all families): Newton vs L-BFGS coefficients, max |difference| {np.max(np.abs(th - fm['L2'])):.2e}")
    sM = th[0] + me["X"].astype(np.float64) @ th[1:]
    sE = th[0] + XmR[:52] @ th[1:]
    P(f"C. exponent scores from own features vs stored: max |difference| {np.max(np.abs(sE - fm['mscore_L2'][chk[:52]])):.2e}")
    Kd = 10000
    base = [min(EG * math.log2(2 * Kd * int(x)) / int(x), 0.95) if s_ else 0.0 for x, s_ in zip(k, surv)]
    pop = [i for i in range(len(k)) if st[i] != -2 and surv[i]]
    allw = sum(W[i] for i in pop)
    srt = sorted((sM[i], W[i]) for i in pop)
    import bisect
    keys = [s for s, _ in srt]; cum = [0.0]
    for _, ww in srt:
        cum.append(cum[-1] + ww)

    def pc(s):
        lo, hi = bisect.bisect_left(keys, s), bisect.bisect_right(keys, s)
        return 100 * (cum[lo] + 0.5 * (cum[hi] - cum[lo])) / allw
    ex = [i for i in range(len(k)) if st[i] == -1]
    mp = sum(pc(sM[i]) for i in ex) / len(ex)
    gM = 0.0
    for i in pop:
        wb = base[i]; zz = math.log(wb / (1 - wb)) + sM[i]; pmod = 1 / (1 + math.exp(-zz))
        if st[i] == -1:
            gM += math.log2(pmod / wb)
        else:
            gM += W[i] * math.log2((1 - pmod) / (1 - wb))
    newi = int(np.flatnonzero(st == -2)[0])
    P(f"C. L2 transfer recomputed: mean percentile of the 51 exponents {mp:.2f}, gain {gM:+.4f} bits, "
      f"52nd exponent percentile {pc(sM[newi]):.2f}, expected count {sum(W[i] * base[i] for i in pop):.2f}")
    # TF survival on 300 rows, plain pow
    sel = rng.choice(np.flatnonzero(st >= 15), 300, replace=False)
    bad = 0
    for i in sel:
        kk = int(k[i]); sv = True
        for j in range(1, Kd + 1):
            q = 2 * j * kk + 1
            if q % 8 in (1, 7) and (kk >= 64 or q < (1 << kk) - 1) and pow(2, kk, q) == 1:
                sv = False; break
        bad += sv != bool(surv[i])
    P(f"C. trial-factoring survival rechecked with pow on 300 population rows: {bad} disagreements")
    open(os.path.join(HERE, "recheck.txt"), "w").write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
