"""Blind test of the Sym^2 family density.  Usage: test_holdout.py PREDICTIONS_JSON
Refuses to run unless the sha256 of the JSON is one of the two frozen hashes below.  Reads the data from data_dir (relative to the
JSON's folder), files data_k{k}_N{N}.npz written by gen_data.py; checks levels, prime lists and dimensions against the JSON;
computes the dim-weighted bin means of sqrt(p) (s/dim - 1) with empirical bin errors, and applies the frozen decisions."""
import sys, os, json, hashlib
import numpy as np

FROZEN = {"95c6e5a4ca0aa70c3d696e1070cdb5b6df86f22bd5ab8a9ca69ae3c6fb3f735d": "holdout",
          "75a1cceb5b88edc06c81fd7bf6f72b1153214476e72c0a793a68366b5fb8f57b": "dry"}
CODE = {"gen_data.py": "989c58a2e7cacabd9b9a9578f156ad9f05a3c3393c6ac307d48823841a06276d",
        "trace.py": "c28d958f0d94523d75dc8e45d6e617f18c81e7fb5ac30964f82fb8c84f78a3d0",
        "hurwitz.py": "fb81201f0691efe80f93362554f0cbfbc1750a74dfd03f769d7095f9d77cd51d"}
here = os.path.dirname(os.path.abspath(__file__))
for fn, hh in CODE.items():
    if hashlib.sha256(open(os.path.join(here, fn), "rb").read()).hexdigest() != hh:
        sys.exit("REFUSED: %s differs from the frozen data-generation code" % fn)
path = sys.argv[1]
raw = open(path).read()
h = hashlib.sha256(raw.encode()).hexdigest()
if h not in FROZEN:
    sys.exit("REFUSED: sha256 %s of %s is not a frozen prediction file" % (h, path))
doc = json.loads(raw)
assert doc["set"] == FROZEN[h]
ddir = os.path.join(os.path.dirname(os.path.abspath(path)), doc["data_dir"])
CHI2_995 = 91.95
res = dict(set=doc["set"], sha256=h, families={})
allpass = True
falsified = False
for fam in doc["families"]:
    k, edges = fam["k"], np.array(fam["bin_edges"]); nb = len(edges) - 1
    rows = []
    for N, dim0 in zip(fam["levels"], fam["dims"]):
        f = os.path.join(ddir, "data_k%d_N%d.npz" % (k, N))
        if not os.path.exists(f):
            sys.exit("MISSING data file %s" % f)
        d = np.load(f)
        assert int(d["N"]) == N and int(d["k"]) == k and int(d["dim"]) == dim0, (N, k)
        p = d["p"].astype(float)
        assert len(p) > 0 and int(p.max()) <= doc["ymax"] * N * N and N not in d["p"]
        y = p / N ** 2; v = np.sqrt(p) * (d["s"] / dim0 - 1.0)
        b = np.digitize(y, edges) - 1; ok = (b >= 0) & (b < nb)
        rows.append((b[ok], v[ok], np.full(ok.sum(), float(dim0))))
    B = np.concatenate([r[0] for r in rows]); V = np.concatenate([r[1] for r in rows]); W = np.concatenate([r[2] for r in rows])
    cnt = np.bincount(B, minlength=nb)
    assert list(cnt) == fam["n_pairs"], "prime lists differ from the frozen ones"
    sw = np.bincount(B, weights=W, minlength=nb); Db = np.bincount(B, weights=W * V, minlength=nb) / sw
    se = np.sqrt(np.bincount(B, weights=W ** 2 * (V - Db[B]) ** 2, minlength=nb)) / sw
    pred = {m: np.array(x) for m, x in fam["pred"].items()}
    base = pred["base_plancherel"]; F = pred["formula"] - base
    a = np.sum((Db - base) * F / se ** 2) / np.sum(F ** 2 / se ** 2); sa = 1 / np.sqrt(np.sum(F ** 2 / se ** 2))
    chi2 = {m: float(np.sum(((Db - pred[m]) / se) ** 2)) for m in pred if m == "formula" or m.startswith("C")}
    S = pred["smooth_part"] - base; C = pred["classnumber_part"] - base
    X = np.vstack([S, C]).T / se[:, None]; cov = np.linalg.inv(X.T @ X); ab = cov @ X.T @ ((Db - base) / se)
    r = dict(k=k, levels=len(fam["levels"]), pairs=int(cnt.sum()), amplitude=float(a), amplitude_sigma=float(sa),
             chi2=chi2, twocomp=[float(ab[0]), float(np.sqrt(cov[0, 0])), float(ab[1]), float(np.sqrt(cov[1, 1]))],
             se_ratio_median=float(np.median(se / np.array(fam["pred_se"]))),
             bins=dict(data=Db.tolist(), se=se.tolist()))
    dec = {}
    dec["P1/P2 amplitude"] = bool(abs(a - 1) <= 3 * sa and a / sa >= 8)
    dec["P3 chi2"] = bool(chi2["formula"] <= CHI2_995)
    if k == 4:
        dec["P4 two-component"] = bool(abs(ab[0] - 1) <= 3 * np.sqrt(cov[0, 0]) and abs(ab[1] - 1) <= 3 * np.sqrt(cov[1, 1])
                                       and ab[1] / np.sqrt(cov[1, 1]) >= 5)
    strong = [m for m, lam in fam["power"].items() if lam >= 25]
    dchi = {m: chi2[m] - chi2["formula"] for m in strong}
    dec["P5 comparators"] = bool(all(v >= 9 for v in dchi.values()))
    r["delta_chi2"] = dchi; r["decisions"] = dec
    r["S1 se ratio in [0.85,1.15]"] = bool(0.85 <= r["se_ratio_median"] <= 1.15)
    allpass &= all(dec.values())
    if abs(a - 1) > 3 * sa or chi2["formula"] > CHI2_995 or any(v < 0 for v in dchi.values()):
        falsified = True
    res["families"][fam["name"]] = r
    print("%s (k=%d, %d levels, %d pairs): amplitude %.3f +- %.3f (predicted sigma %.3f); chi2(formula) %.1f / %d bins; se ratio %.3f" % (
        fam["name"], k, len(fam["levels"]), cnt.sum(), a, sa, fam["amp_sigma"], chi2["formula"], nb, r["se_ratio_median"]))
    print("   two-component: a = %.3f +- %.3f, b = %.3f +- %.3f" % tuple(r["twocomp"]))
    for m in sorted(dchi):
        print("   delta chi2 vs %-34s %9.1f   (predicted %.1f)" % (m, dchi[m], fam["power"][m]))
    print("   decisions:", dec)
res["verdict"] = "CONFIRMED" if allpass else ("FALSIFIED" if falsified else "NOT CONFIRMED")
print("VERDICT (%s): %s" % (doc["set"], res["verdict"]))
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(path)), "test_result_%s.json" % doc["set"]), "w"), indent=1)
