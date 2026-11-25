"""tz_test.py -- frozen test of the linearised explicit formula for prime patterns on the holdout (2e9, 2e10].

Run (after the holdout counts exist):
    python tz_counts.py hold 2e9 2e10 4e9 4e10        (repeat until it prints DONE; about 25 min in total)
    python tz_test.py                                  (real test: disc.npz + hold.npz, window 2e9 .. 2e10)
    python tz_test.py --dry                            (dry run on the discovery range only: window 2e8 .. 2e9)
The script refuses to run if predictions_twinzeros.json or any code/data file differs from the registered hashes.
"""
import sys, json, hashlib, os
import numpy as np

EXPECTED_JSON_SHA256 = "f525b850ae0bbffdda8b2d9c63cb0c9f61521c81fb5456abe9134c958ed51f2c"
HERE = os.path.dirname(os.path.abspath(__file__))
ZDIR = "/home/bas/praime/writeup/data/riemann/"


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def check_hashes(pred):
    jp = os.path.join(HERE, "predictions_twinzeros.json")
    if sha(jp) != EXPECTED_JSON_SHA256:
        sys.exit("REFUSED: predictions_twinzeros.json does not match the registered hash")
    for name, h in pred["hashes"].items():
        p = name if name.startswith("/") else os.path.join(HERE, name)
        if sha(p) != h:
            sys.exit("REFUSED: %s does not match the registered hash" % name)


def comb(gammas, uu, chunk=200):
    out = np.zeros(len(uu))
    for i in range(0, len(gammas), chunk):
        g = gammas[i:i + chunk]
        out += -2 * np.real(np.exp(1j * np.outer(uu, g)) @ (1 / (0.5 + 1j * g)))
    return out


def compute(D, j0, j1, tzlib, models, theory):
    out = {}
    band_full, band_low = (2.0, 6000.0), (2.0, 300.0)
    for name in ("h2", "h4", "h6", "h8", "h10", "h12", "h30", "SG", "T026", "T046"):
        Y, S = D.Ypat(name)
        k = models.kforms(name)
        for tag, band, null in (("full", band_full, models.NULL_SET), ("low", band_low, [])):
            if tag == "low" and name != "h2":
                continue
            S_, regs = models.build(D, name, null_set=null)
            b, cov, d = tzlib.fit(Y, [r[1] for r in regs], D.u, j0, j1, band[0], band[1], k)
            e = np.sqrt(np.diag(cov))
            res = {r[0]: [float(x), float(s), float(r[2])] for r, x, s in zip(regs, b, e)}
            idx = [i for i, r in enumerate(regs) if r[0].startswith("null")]
            if idx:
                bb = b[idx]
                res["null_chi2"] = [float(bb @ np.linalg.solve(cov[np.ix_(idx, idx)], bb)), len(idx)]
            res["chi2_dof"] = d["chi2_dof"]
            out["%s/%s/all" % (name, tag)] = res
        S_, regs = models.build(D, name, null_set=[], per_conductor=True)
        b, cov, d = tzlib.fit(Y, [r[1] for r in regs], D.u, j0, j1, band_full[0], band_full[1], k)
        e = np.sqrt(np.diag(cov))
        out["%s/full/cond" % name] = {r[0]: [float(x), float(s), float(r[2])] for r, x, s in zip(regs, b, e)}
    # raw chi_-3 coefficient of the triples (sign test): model zeta + chi3 alone + all other present chi
    for name in ("T026", "T046"):
        Y, S = D.Ypat(name)
        S_, rows = theory.prediction_table(name, D.chars)
        r3 = [r for r in rows if r["label"] == "3.1"][0]
        others = models.chi_combination(D, [r for r in rows if r["label"] != "3.1"], models.QS, 1)
        regs = [D.Yzeta(1), D.Ychi("3.1").real, others]
        b, cov, d = tzlib.fit(Y, regs, D.u, j0, j1, band_full[0], band_full[1], 3)
        out["%s/full/raw_chi3" % name] = [float(b[1]), float(np.sqrt(cov[1, 1])), float(r3["beta"][1].real)]
    # explicit zeros (approach B) for twins: combs of zeta, chi_5, chi_-3, chi_-4 zeros, band [2, 1900]
    zz = np.loadtxt(ZDIR + "zeros_40000_flint.txt"); zz = zz[zz <= 2500]
    z5 = np.loadtxt(ZDIR + "lzeros/D5.txt"); z5 = z5[(z5 > 0) & (z5 <= 2500)]
    z3 = np.loadtxt(ZDIR + "lzeros/D-3.txt"); z3 = z3[(z3 > 0) & (z3 <= 2500)]
    z4 = np.loadtxt(ZDIR + "lzeros/D-4.txt"); z4 = z4[(z4 > 0) & (z4 <= 2500)]
    uu = D.u.copy()
    sl = slice(j0, j1 + 1)
    combs = []
    for zl in (zz, z5, z3, z4):
        c = np.zeros(len(uu)); c[sl] = comb(zl, uu[sl]); combs.append(c)
    Y, S = D.Ypat("h2")
    b, cov, d = tzlib.fit(Y, combs, D.u, j0, j1, 2.0, 1900.0, 2)
    e = np.sqrt(np.diag(cov))
    preds = [2 * S, 2 * S / 3, 0.0, 0.0]
    out["h2/B/combs"] = {lab: [float(x), float(s), p] for lab, x, s, p in zip(("zeta comb", "chi5 comb", "chi-3 comb", "chi-4 comb"), b, e, preds)}
    return out


def z(v):
    return (v[0] - v[2]) / v[1]


def decide(st, S2):
    dec = {}
    a = st["h2/full/all"]
    zt = a["zeta(a=1)"]
    dec["P1 twins zeta coefficient = 2S(2)"] = abs(z(zt)) <= 3 and abs(zt[0] - S2) / zt[1] > 5
    lam = a["all present chi (loudness)"]
    dec["P2 twins joint L-function loudness = 1 and present"] = abs(z(lam)) <= 3 and lam[0] / lam[1] >= 5
    c3 = a["null 3.1"]
    dec["P3 twins: chi_-3 absent (and |beta| = 2S rejected)"] = abs(c3[0]) / c3[1] <= 3 and min(abs(c3[0] - 2 * S2), abs(c3[0] + 2 * S2)) / c3[1] > 5
    dec["P4 twins: null set chi2 below 99.9% quantile (34.53, 13 dof)"] = a["null_chi2"][0] <= 34.53
    zs = []
    for name in ("h4", "h6", "h8", "h10", "h12", "h30"):
        r = st["%s/full/all" % name]
        zs += [z(r["zeta(a=1)"]), z(r["all present chi (loudness)"])]
    dec["S1 other gaps: all 12 |z| <= 3.5"] = max(abs(x) for x in zs) <= 3.5
    sg = st["SG/full/all"]; sgc = st["SG/full/cond"]
    zz_ = [z(sg["zeta(a=1)"]), z(sg["zeta(a=2)"]), z(sg["all present chi (loudness)"]), z(sgc["chi q=3 (loudness)"]), z(sgc["chi q=4 (loudness)"])]
    dec["S2 Sophie Germain: zeta at both scales, chi_-3 and chi_-4 loud, all |z| <= 3"] = max(abs(x) for x in zz_) <= 3 and sgc["chi q=3 (loudness)"][0] / sgc["chi q=3 (loudness)"][1] >= 5
    ok = True
    for name in ("T026", "T046"):
        r = st["%s/full/cond" % name]["chi q=3 (loudness)"]
        raw = st["%s/full/raw_chi3" % name]
        ok = ok and abs(z(r)) <= 3 and np.sign(raw[0]) == np.sign(raw[2]) and abs(raw[0]) / raw[1] >= 2
    dec["S3 triples: chi_-3 with opposite signs for (0,2,6) and (0,4,6)"] = bool(ok)
    lo = st["h2/low/all"]
    dec["S4 twins low band [2,300]: zeta and loudness within 3 sigma"] = abs(z(lo["zeta(a=1)"])) <= 3 and abs(z(lo["all present chi (loudness)"])) <= 3
    B = st["h2/B/combs"]
    dec["S5 twins on explicit zeros (zeta, chi5, chi-3, chi-4 combs): all |z| <= 3"] = max(abs(z(v)) for v in B.values()) <= 3
    c = st["h2/full/cond"]
    dec["S6 twins per conductor 5, 7, 15, 21: loudness within 3 sigma"] = max(abs(z(c["chi q=%d (loudness)" % q])) for q in (5, 7, 15, 21)) <= 3
    return {k: bool(v) for k, v in dec.items()}


def main():
    dry = "--dry" in sys.argv
    pred = json.load(open(os.path.join(HERE, "predictions_twinzeros.json")))
    if not dry or EXPECTED_JSON_SHA256 != "__FILL__":
        check_hashes(pred)
    sys.path.insert(0, HERE)
    import tzlib, models, theory
    if dry:
        D = models.Data([os.path.join(HERE, "disc.npz")])
        j0, j1 = D.J(2 * 10 ** 8), D.J(2 * 10 ** 9, right=True) - 1
    else:
        if not os.path.exists(os.path.join(HERE, "hold.npz")):
            sys.exit("hold.npz missing: run tz_counts.py hold 2e9 2e10 4e9 4e10 until DONE")
        D = models.Data([os.path.join(HERE, "disc.npz"), os.path.join(HERE, "hold.npz")])
        j0, j1 = D.J(2 * 10 ** 9), D.J(2 * 10 ** 10, right=True) - 1
    print("window x = %d .. %d (u %.4f .. %.4f), %d grid points" % (D.X[j0], D.X[j1], D.u[j0], D.u[j1], j1 - j0 + 1))
    st = compute(D, j0, j1, tzlib, models, theory)
    S2 = theory.singular_const("h2")
    dec = decide(st, S2)
    for key in sorted(st):
        v = st[key]
        if isinstance(v, dict):
            for lab, x in v.items():
                if isinstance(x, list) and len(x) == 3:
                    print("  %-18s %-30s measured %9.4f +- %7.4f  predicted %8.4f  z %+6.2f" % (key, lab, x[0], x[1], x[2], z(x)))
                elif lab == "null_chi2":
                    print("  %-18s %-30s chi2 %.2f on %d dof" % (key, lab, x[0], x[1]))
        else:
            print("  %-18s raw chi_-3 coefficient %9.4f +- %7.4f  predicted %8.4f  z %+6.2f" % (key, v[0], v[1], v[2], z(v)))
    print("\nDECISIONS (%s):" % ("DRY RUN on 2e8..2e9, not the test" if dry else "HOLDOUT 2e9..2e10"))
    for k_, v in dec.items():
        print("  %-80s %s" % (k_, "PASS" if v else "FAIL"))
    primary = all(v for k_, v in dec.items() if k_.startswith("P"))
    print("PRIMARY VERDICT:", "CONFIRMED" if primary else "FALSIFIED")
    json.dump({"dry": dry, "stats": st, "decisions": dec, "primary_confirmed": primary},
              open(os.path.join(HERE, "tz_test_result%s.json" % ("_dry" if dry else "")), "w"), indent=1)


if __name__ == "__main__":
    main()
