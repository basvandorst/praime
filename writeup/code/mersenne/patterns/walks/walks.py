"""Walk layouts (mertens, race4, liouville): page check, measurements, random-walk models.
Output: walks_out.json and printed summary."""
import json, math, sys
import numpy as np
from common import series, load_pos, N, VIS

NB = 10 ** 6
S = series(NB)
out = {}

# ---------- 1. the page draws exactly these series (N = 500000) ----------
def page_series(layout, lo, hi):
    fx, fy = load_pos(layout)
    sp = max(1, hi - lo)
    Yp = hi - (fy.astype(np.float64) - 2) / 700 * sp
    n = np.arange(N + 1)
    xexp = (n - 1) / (N - 1) * 1400 + 2
    return np.rint(Yp).astype(np.int64), np.max(np.abs(fx[1:] - xexp[1:])), np.max(np.abs(Yp - np.rint(Yp))[1:])

for layout, key in [("mertens", "M"), ("race4", "R"), ("liouville", "L")]:
    W = S[key][: N + 1]
    lo, hi = int(W[1:].min()), int(W[1:].max())
    js = json.load(open(f"{VIS}/data/{layout}.json"))
    Yp, dxmax, frac = page_series(layout, lo, hi)
    out[f"{layout}_page_check"] = dict(lo=lo, hi=hi, page_lo=js["extra"]["walk"]["lo"],
                                       page_hi=js["extra"]["walk"]["hi"],
                                       mismatches=int(np.sum(Yp[1:] != W[1:])), max_dx=float(dxmax),
                                       max_round=float(frac))

# ---------- 2. measurements ----------
n_all = np.arange(NB + 1, dtype=np.float64)

def walk_stats(W, n_hi, norm, n0=1000):
    """W indexed by n; norm(n) the normalising scale."""
    n = n_all[n0: n_hi + 1]
    w = W[n0: n_hi + 1] / norm(n)
    lw = 1.0 / n
    mean = float(np.sum(w * lw) / np.sum(lw))
    rms = float(np.sqrt(np.sum(w * w * lw) / np.sum(lw)))
    sd = float(np.sqrt(max(0.0, rms ** 2 - mean ** 2)))
    Wf = W[1: n_hi + 1]
    zeros = int(np.sum(Wf == 0))
    s = np.sign(Wf)
    s = s[s != 0]
    changes = int(np.sum(s[1:] != s[:-1]))
    pos = float(np.mean(Wf > 0)); neg = float(np.mean(Wf < 0))
    i = int(np.argmax(np.abs(w)))
    return dict(logmean=mean, logrms=rms, logsd=sd, max_abs=float(abs(w[i])), argmax=int(n[i]),
                zeros=zeros, sign_changes=changes, frac_pos=pos, frac_neg=neg)

sq = lambda n: np.sqrt(n)
rnorm = lambda n: np.sqrt(n) / np.log(n)

for key, name, norm in [("M", "mertens", sq), ("L", "liouville", sq), ("R", "race4", rnorm)]:
    W = S[key]
    out[f"{name}_stats_N"] = walk_stats(W, N, norm)
    out[f"{name}_stats_1e6"] = walk_stats(W, NB, norm)
    # all n >= 2 for the max ratio
    w2 = W[2: N + 1] / norm(n_all[2: N + 1])
    i = int(np.argmax(np.abs(w2)))
    out[f"{name}_max_ratio_n>=2"] = dict(value=float(w2[i]), n=i + 2)
    out[f"{name}_at_powers10"] = {str(10 ** k): int(W[10 ** k]) for k in range(1, 7)}
    out[f"{name}_end"] = int(W[N])
    # scaling: RMS in dyadic blocks and running max
    js_, rms_, mx_ = [], [], []
    for j in range(8, 19):
        a, b = 2 ** j, min(2 ** (j + 1), N)
        js_.append(math.sqrt(a * b)); rms_.append(float(np.sqrt(np.mean(W[a:b].astype(float) ** 2))))
        mx_.append(float(np.max(np.abs(W[1:b]))))
    lx = np.log(js_)
    out[f"{name}_scaling"] = dict(slope_blockrms=float(np.polyfit(lx, np.log(rms_), 1)[0]),
                                  slope_runmax=float(np.polyfit(lx, np.log(mx_), 1)[0]),
                                  blocks=[[round(x), round(r, 2)] for x, r in zip(js_, rms_)])
    # structure function on [N/4, N]: RMS of W(n+h)-W(n)
    hs = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000]
    sf = []
    for h in hs:
        a = N // 4
        d = (W[a + h: N + 1] - W[a: N + 1 - h]).astype(float)
        sf.append(float(np.sqrt(np.mean(d * d))))
    out[f"{name}_structure"] = dict(h=hs, rms=sf,
                                    slope_small=float(np.polyfit(np.log(hs[:6]), np.log(sf[:6]), 1)[0]),
                                    slope_large=float(np.polyfit(np.log(hs[7:]), np.log(sf[7:]), 1)[0]))

# Mertens: zero counts vs A171910 (zeros of M in 0<x<10^k)
M = S["M"]
out["mertens_zero_counts"] = {str(10 ** k): int(np.sum(M[1: 10 ** k] == 0)) for k in range(1, 7)}
out["A171910"] = {"10": 1, "100": 6, "1000": 92, "10000": 406, "100000": 1549, "1000000": 5361}
out["A084237"] = {"10": -1, "100": 1, "1000": 2, "10000": -23, "100000": -48, "1000000": 212}
# Liouville: zeros, max for n>=2, A090410
L = S["L"]
out["liouville_zeros_le_1e6"] = [int(x) for x in np.nonzero(L[1:] == 0)[0] + 1]
out["liouville_max_n>=2_le_1e6"] = int(L[2:].max())
out["A090410"] = {"10": -2, "100": -14, "1000": -94, "10000": -288, "100000": -530, "1000000": -842}
out["one_over_zeta_half"] = -0.684765
# race: leader changes and ties vs A007350 / A007351
R = S["R"]
isp = S["isp"]
lead_changes = []
last = 0
nz = np.nonzero(R)[0]
sgn = np.sign(R[nz])
ch = np.nonzero(sgn[1:] != sgn[:-1])[0] + 1
lead_changes = [int(nz[0])] + [int(nz[i]) for i in ch]
out["race_leader_changes_le_1e6"] = lead_changes
out["race_leader_changes_count_le_1e6"] = len(lead_changes)
A007350 = [3, 26861, 26879, 616841, 617039, 617269, 617471, 617521, 617587, 617689, 617723, 622813,
           623387, 623401, 623851, 623933, 624031, 624097, 624191, 624241, 624259, 626929, 626963,
           627353, 627391, 627449, 627511, 627733, 627919, 628013, 628427, 628937, 629371]
out["race_A007350_prefix_match"] = lead_changes[:len(A007350)] == A007350
ties = [int(p) for p in S["primes"] if R[p] == 0]
A007351 = [2, 5, 17, 41, 461, 26833, 26849, 26863, 26881, 26893, 26921, 616769, 616793, 616829, 616843,
           616871, 617027, 617257, 617363, 617387, 617411, 617447, 617467, 617473, 617509, 617531, 617579,
           617681, 617707, 617719, 618437, 618521, 618593, 618637]
out["race_A007351_prefix_match"] = ties[:len(A007351)] == A007351
out["race_ties_le_N"] = [t for t in ties if t <= N]
out["race_n_with_R<0_le_N"] = int(np.sum(R[1: N + 1] < 0))
out["race_n_with_R<0_le_1e6"] = int(np.sum(R[1:] < 0))
out["race_n_with_R=0_le_N"] = int(np.sum(R[1: N + 1] == 0))
# race vs 1/2 li(sqrt n)
import mpmath
for x in [10 ** 4, 10 ** 5, N, 10 ** 6]:
    out[f"race_vs_half_li_sqrt_{x}"] = [int(R[x]), float(mpmath.li(math.sqrt(x)) / 2), float(math.sqrt(x) / math.log(x))]

# ---------- 3. random-walk models with the same step distribution ----------
rng = np.random.default_rng(20261115)
NS = 400
sqf = (S["mu"][: N + 1] != 0)
oddp = S["isp"][: N + 1].copy(); oddp[2] = False
models = {"mertens": (sqf, sq), "liouville": (np.ones(N + 1, bool), sq), "race4": (oddp, rnorm)}
models["liouville"][0][0] = False
for name, (mask, norm) in models.items():
    idx = np.nonzero(mask)[0]
    rec = {k: [] for k in ["logmean", "logrms", "max_abs", "zeros", "sign_changes", "frac_pos"]}
    for b in range(NS // 40):
        st = np.zeros((40, N + 1), dtype=np.int8)
        st[:, idx] = rng.choice(np.array([-1, 1], dtype=np.int8), size=(40, idx.size))
        if name == "mertens":
            st[:, 1] = 1  # mu(1)=1 fixed
        if name == "liouville":
            st[:, 1] = 1
        Wb = np.cumsum(st, axis=1, dtype=np.int32)
        del st
        for r in range(40):
            s = walk_stats(Wb[r].astype(np.int64), N, norm)
            for k in rec:
                rec[k].append(s[k])
        del Wb
    summ = {}
    real = out[f"{name}_stats_N"]
    for k, v in rec.items():
        v = np.array(v, float)
        summ[k] = dict(model_mean=float(v.mean()), model_sd=float(v.std()),
                       model_q05=float(np.quantile(v, 0.05)), model_q95=float(np.quantile(v, 0.95)),
                       real=real[k], frac_model_le_real=float(np.mean(v <= real[k])))
    out[f"{name}_rw_model"] = summ

json.dump(out, open("walks_out.json", "w"), indent=1)
for k, v in out.items():
    if isinstance(v, list) and len(v) > 40:
        v = v[:40]
    print(k, ":", v)
