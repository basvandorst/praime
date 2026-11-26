"""Own re-implementation (review): M_k(y) and nu(m) from the STUDY.md / derivation.md formulas, own binning, amplitude,
chi2, jackknife over levels, region split, class-number-part tests, halves of levels."""
import glob, os, json
import numpy as np
from scipy.special import eval_chebyu

DD = "/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/sym2murmur/out_holdout"
OUT = os.path.dirname(os.path.abspath(__file__))

# primes for Euler products
def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]
PR = primes_upto(2_000_000).astype(float)
u = 1.0 / ((PR - 1) * (PR ** 2 - 1))
C0 = float(np.exp(np.log1p(-u).sum()))

def nu_l(l, e):
    if e == 0: return 1 - 1.0 / ((l - 1) * (l * l - 1))
    v = 1 + 1.0 / (l * l - 1) + (l / (l + 1.0)) * sum(l ** (-j) for j in range(1, e // 2 + 1))
    if e % 2: v -= l ** (-(e - 1) / 2) / (l * l - 1.0)
    return v

def fac(m):
    f = {}; d = 2
    while d * d <= m:
        while m % d == 0: f[d] = f.get(d, 0) + 1; m //= d
        d += 1
    if m > 1: f[m] = f.get(m, 0) + 1
    return f

def nu(m):
    v = C0
    for l, e in fac(m).items(): v *= nu_l(l, e) / nu_l(l, 0)
    return v

def nu_sqf(m):  # only which primes divide m, exponent ignored (nu_l(e) := nu_l(1) = 1)
    v = C0
    for l in fac(m): v /= nu_l(l, 0)
    return v

def nu_zub(m):  # Zubrilina's local factor for l | r, transplanted to m (shape only; scale fitted)
    v = 1.0
    for l in fac(m): v *= 1 + l * l / (l ** 4 - 2 * l * l - l + 1.0)
    return v

MMAX = 12
NU = {"nu": [nu(m) for m in range(1, MMAX + 1)], "one": [1.0] * MMAX, "sqf": [nu_sqf(m) for m in range(1, MMAX + 1)],
      "zub": [nu_zub(m) for m in range(1, MMAX + 1)]}

def basis(y, k):
    """per-pair: smooth part S(y) and class-number basis B_m(y) (so CN part = sum_m nu(m) B_m)."""
    S = (24 * y ** 1.5 if k == 2 else 0 * y) - 6 * np.sqrt(y) / (k - 1)
    B = np.zeros((MMAX, len(y)))
    for m in range(1, MMAX + 1):
        msk = 4 * y > m
        B[m - 1, msk] = -12 / (np.pi * (k - 1) * np.sqrt(y[msk])) * np.sqrt(m * (4 * y[msk] - m)) * eval_chebyu(k - 2, 1 - m / (2 * y[msk]))
    return S, B

edges = np.arange(0, 3.0 + 1e-9, 0.05); NB = len(edges) - 1

def load(k):
    L = []
    for f in sorted(glob.glob(DD + "/data_k%d_N*.npz" % k)):
        d = np.load(f); N = int(d["N"]); dim = int(d["dim"]); p = d["p"].astype(float); y = p / N ** 2
        V = np.sqrt(p) * (d["s"] / dim - 1.0); rho = (k - 1) * N / (12.0 * dim)
        b = np.digitize(y, edges) - 1; ok = (b >= 0) & (b < NB)
        S, B = basis(y[ok], k)
        L.append(dict(N=N, dim=dim, b=b[ok], V=V[ok], w=np.full(ok.sum(), float(dim)), rho=rho, p=p[ok], y=y[ok],
                      base=rho / np.sqrt(p[ok]), S=rho * S, B=rho * B))
    return L

def binned(L):
    b = np.concatenate([x["b"] for x in L]); w = np.concatenate([x["w"] for x in L]); V = np.concatenate([x["V"] for x in L])
    sw = np.bincount(b, w, NB); D = np.bincount(b, w * V, NB) / sw
    se = np.sqrt(np.bincount(b, w ** 2 * (V - D[b]) ** 2, NB)) / sw
    base = np.bincount(b, w * np.concatenate([x["base"] for x in L]), NB) / sw
    S = np.bincount(b, w * np.concatenate([x["S"] for x in L]), NB) / sw
    Bm = np.array([np.bincount(b, w * np.concatenate([x["B"][m] for x in L]), NB) / sw for m in range(MMAX)])
    return dict(D=D, se=se, base=base, S=S, Bm=Bm, n=np.bincount(b, minlength=NB))

def amp(D, se, base, F):
    a = np.sum((D - base) * F / se ** 2) / np.sum(F ** 2 / se ** 2); return a, 1 / np.sqrt(np.sum(F ** 2 / se ** 2))

res = {}
for k in (2, 4):
    L = load(k); X = binned(L)
    D, se, base, S, Bm = X["D"], X["se"], X["base"], X["S"], X["Bm"]
    CN = NU["nu"] @ Bm; F = S + CN; pred = base + F
    a, sa = amp(D, se, base, F); chi2 = np.sum(((D - pred) / se) ** 2)
    # jackknife over levels
    aj = []
    for i in range(len(L)):
        Xi = binned(L[:i] + L[i + 1:]); aj.append(amp(Xi["D"], Xi["se"], Xi["base"], Xi["S"] + NU["nu"] @ Xi["Bm"])[0])
    aj = np.array(aj); n = len(aj); sj = np.sqrt((n - 1) / n * np.sum((aj - aj.mean()) ** 2))
    # per-level amplitude scatter (each level alone), normalized
    al = [amp(*[binned([x])[q] for q in ("D", "se", "base")], binned([x])["S"] + NU["nu"] @ binned([x])["Bm"]) for x in L]
    zl = np.array([(ai - a) / si for ai, si in al])
    r = dict(k=k, levels=[x["N"] for x in L], pairs=int(X["n"].sum()), amp=a, amp_se=sa, amp_jack=sj, chi2=chi2,
             level_amp_z_rms=float(np.sqrt(np.mean(zl ** 2))))
    # region split of delta chi2 vs null (model = base)
    lo = edges[1:] <= 0.25 + 1e-12
    dnull = ((D - base) ** 2 - (D - pred) ** 2) / se ** 2
    r["dchi2_null_total"] = float(dnull.sum()); r["dchi2_null_y<1/4"] = float(dnull[lo].sum()); r["dchi2_null_y>1/4"] = float(dnull[~lo].sum())
    # smooth-only on y>1/4 (formula and smooth-only coincide below 1/4)
    dsm = ((D - base - S) ** 2 - (D - pred) ** 2) / se ** 2
    r["dchi2_smoothonly_vs_formula"] = float(dsm.sum())
    # y<1/4 alone: amplitude of smooth part
    r["amp_y<1/4"] = amp(D[lo], se[lo], base[lo], S[lo])
    r["amp_y>1/4"] = amp(D[~lo], se[~lo], base[~lo], F[~lo])
    # class-number part tests on residual R = D - base - S (smooth fixed at prediction), bins y > 1/4
    R = (D - base - S)[~lo]; s2 = se[~lo]; Bh = Bm[:, ~lo]
    cn = {}
    for nm, wts in NU.items():
        G = np.array(wts) @ Bh
        c, sc = amp(R, s2, 0 * R, G)
        cn[nm] = dict(scale=c, scale_se=sc, chi2_fixed1=float(np.sum(((R - G) / s2) ** 2)), chi2_scaled=float(np.sum(((R - c * G) / s2) ** 2)))
    cn["none"] = dict(chi2=float(np.sum((R / s2) ** 2)))
    # free fit nu(1..11)
    A = (Bh[:11] / s2).T; cov = np.linalg.inv(A.T @ A); nh = cov @ A.T @ (R / s2)
    sig = np.sqrt(np.diag(cov))
    cn["free"] = dict(nu_hat=nh.tolist(), sigma=sig.tolist(), chi2=float(np.sum(((R - nh @ Bh[:11]) / s2) ** 2)))
    for nm in ("nu", "one", "sqf"):
        dlt = nh - np.array(NU[nm][:11]); cn["free"]["mahal_vs_" + nm] = float(dlt @ np.linalg.inv(cov) @ dlt)
    r["cn"] = cn
    # halves
    h = len(L) // 2
    for nm, part in (("low", L[:h]), ("high", L[h:])):
        Y = binned(part); FF = Y["S"] + NU["nu"] @ Y["Bm"]; aa, ss = amp(Y["D"], Y["se"], Y["base"], FF)
        rr = (Y["D"] - Y["base"] - FF) / Y["se"]
        r["half_" + nm] = dict(levels=[x["N"] for x in part], amp=aa, amp_se=ss, chi2=float(np.sum(rr ** 2)))
        r["half_" + nm]["resid_mean_by_third"] = [float(rr[i * 20:(i + 1) * 20].mean() * np.sqrt(20)) for i in range(3)]
    res[k] = r
    np.savez(os.path.join(OUT, "r159_bins_k%d.npz" % k), D=D, se=se, base=base, S=S, Bm=Bm)

def fmt(o):
    if isinstance(o, float): return round(o, 4)
    if isinstance(o, tuple): return [fmt(x) for x in o]
    if isinstance(o, list): return [fmt(x) for x in o]
    if isinstance(o, dict): return {kk: fmt(v) for kk, v in o.items()}
    return o
print("C0 = %.6f; nu(1..12) =" % C0, " ".join("%.4f" % x for x in NU["nu"]))
for k, r in res.items():
    print(json.dumps(fmt(r)))
json.dump(fmt({str(k): v for k, v in res.items()}), open(os.path.join(OUT, "r159_fit.json"), "w"), indent=1)
