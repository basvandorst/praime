"""Feature tests for Mersenne exponents against the Wagstaff/LPW null.

Null: each prime p < B is an exponent independently with prob P(p) = min(1, e^g log2(a p)/p).
For each binary feature F: O = #exponents in F, E = sum_{p in F} P(p), exact two-sided
Poisson-binomial p-value (exact DP for p < 1e5, Poisson for p >= 1e5 where P < 3e-4).
Continuous statistics: Monte Carlo of the null.
"""
import os, sys, math, json
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from common import EXPS, B, SPLIT, load_primes, null_prob, EULER_GAMMA

rng = np.random.default_rng(20261115)
NSIM = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
LMAX = 140

pr = load_primes()
nB = int(np.searchsorted(pr, B))           # primes < B
p_all = pr[:nB].astype(np.int64)
EXPS_A = np.array(EXPS, dtype=np.int64)
assert EXPS_A.max() < B and np.all(np.isin(EXPS_A, p_all))
P = null_prob(p_all, "mod4")
P2 = null_prob(p_all, "a2")
isexp = np.isin(p_all, EXPS_A)
small = p_all < SPLIT
log = []


def out(s=""):
    print(s)
    log.append(s)


def is_prime_arr(q):
    q = np.asarray(q, dtype=np.int64)
    idx = np.searchsorted(pr, q)
    idx = np.minimum(idx, len(pr) - 1)
    assert q.max() <= 2 * B + 1
    return pr[idx] == q


# ---------------- features ----------------
MODS = [3, 5, 7, 8, 12, 16, 24, 30, 60, 120]
NMC = set()
for k in range(0, 30):
    for v in (2**k - 1, 2**k + 1, 4**k - 3, 4**k + 3):
        if 2 <= v < B:
            NMC.add(v)
POW2pm1 = np.array(sorted({v for k in range(0, 30) for v in (2**k - 1, 2**k + 1) if 2 <= v < B}))
FOUR3 = np.array(sorted({v for k in range(1, 15) for v in (4**k - 3, 4**k + 3) if 2 <= v < B}))
NMC_A = np.array(sorted(NMC))

names = []
for m in MODS:
    for r in range(m):
        if math.gcd(r, m) == 1:
            names.append(f"p={r} mod {m}")
names += ["p=+-1 mod 8", "p=2^k+-1", "p=4^k+-3", "NMC form (2^k+-1 or 4^k+-3)",
          "log2 p within 0.1 of integer", "2p+1 prime", "2p+1 prime & p=3 mod 4 & p>3 [control]",
          "2p+1 prime & p=1 mod 4", "p-2 prime", "p+2 prime", "p-2 or p+2 prime"]
NF = len(names)


def popcount(x):
    return np.bitwise_count(np.asarray(x, dtype=np.uint64)).astype(np.int64)


def longest_run(x):
    x = np.asarray(x, dtype=np.uint32).copy()
    c = np.zeros(len(x), dtype=np.int8)
    while np.any(x):
        c += (x != 0)
        x &= (x << np.uint32(1))
    return c


def feature_matrix(p):
    p = np.asarray(p, dtype=np.int64)
    cols = []
    for m in MODS:
        rm = p % m
        for r in range(m):
            if math.gcd(r, m) == 1:
                cols.append(rm == r)
    m8 = p % 8
    cols.append((m8 == 1) | (m8 == 7))
    cols.append(np.isin(p, POW2pm1))
    cols.append(np.isin(p, FOUR3))
    cols.append(np.isin(p, NMC_A))
    fr = np.log2(p.astype(np.float64)) % 1.0
    cols.append(np.minimum(fr, 1 - fr) < 0.1)
    sg = is_prime_arr(2 * p + 1)
    cols.append(sg)
    cols.append(sg & (p % 4 == 3) & (p > 3))
    cols.append(sg & (p % 4 == 1))
    pm2 = is_prime_arr(np.maximum(p - 2, 0))
    pp2 = is_prime_arr(p + 2)
    cols += [pm2, pp2, pm2 | pp2]
    return np.stack(cols, axis=1)


# refined null: Euler's theorem (p=3 mod 4, p>3, 2p+1 prime => 2p+1 | 2^p-1) set to 0;
# other p=3 mod 4 rescaled per octave so the class total is unchanged (the LPW
# formula is an average over whether 2p+1 is a prime divisor).
killed = (p_all % 4 == 3) & (p_all > 3) & is_prime_arr(2 * p_all + 1)
_oct = np.floor(np.log2(p_all.astype(np.float64))).astype(np.int64)
_m3 = (p_all % 4 == 3)
_tot = np.bincount(_oct[_m3], weights=P[_m3], minlength=40)
_keep = np.bincount(_oct[_m3 & ~killed], weights=P[_m3 & ~killed], minlength=40)
fL = np.where(_keep > 0, _tot / np.maximum(_keep, 1e-300), 1.0)
Pref = np.minimum(1.0, np.where(killed, 0.0, np.where(_m3, P * fL[_oct], P)))
del _oct, _m3
assert not np.any(killed & isexp)

# expected counts per feature (full / discovery / holdout), chunked
E_full = np.zeros(NF); E_d = np.zeros(NF); E_h = np.zeros(NF)
E2_full = np.zeros(NF); Eref = np.zeros(NF); Eref_h = np.zeros(NF)
CH = 100000
for s in range(0, nB, CH):
    F = feature_matrix(p_all[s:s+CH]).astype(np.float32)
    Pc = P[s:s+CH]; sm = small[s:s+CH]
    E_full += (F.T @ Pc.astype(np.float64)) if False else F.T.astype(np.float64) @ Pc
    E2_full += F.T.astype(np.float64) @ P2[s:s+CH]
    E_d += F.T.astype(np.float64) @ (Pc * sm)
    E_h += F.T.astype(np.float64) @ (Pc * ~sm)
    Eref += F.T.astype(np.float64) @ Pref[s:s+CH]
    Eref_h += F.T.astype(np.float64) @ (Pref[s:s+CH] * ~sm)
    del F
FE = feature_matrix(EXPS_A)
O_full = FE.sum(0); O_d = FE[EXPS_A < SPLIT].sum(0); O_h = FE[EXPS_A >= SPLIT].sum(0)

# exact Poisson-binomial pmf for p < SPLIT, all features at once
ps_small = p_all[small]; Ps_small = P[small]
Fs = feature_matrix(ps_small).astype(np.float64)
dist = np.zeros((NF, LMAX)); dist[:, 0] = 1.0
for i in range(len(ps_small)):
    q = Ps_small[i] * Fs[i]                      # (NF,)
    sh = np.zeros_like(dist); sh[:, 1:] = dist[:, :-1]
    dist = dist * (1 - q)[:, None] + sh * q[:, None]
pmf_d = dist
Pr_small = Pref[small]
dist = np.zeros((NF, LMAX)); dist[:, 0] = 1.0
for i in range(len(ps_small)):
    q = Pr_small[i] * Fs[i]
    sh = np.zeros_like(dist); sh[:, 1:] = dist[:, :-1]
    dist = dist * (1 - q)[:, None] + sh * q[:, None]
pmf_ref_d = dist


def poisson_pmf(lam, L=LMAX):
    k = np.arange(L)
    from scipy.stats import poisson
    return poisson.pmf(k[None, :], np.asarray(lam)[:, None])


pmf_h = poisson_pmf(E_h)
pmf_full = np.array([np.convolve(pmf_d[f], pmf_h[f])[:LMAX] for f in range(NF)])
pmf_ref_h = poisson_pmf(Eref_h)
pmf_ref = np.array([np.convolve(pmf_ref_d[f], pmf_ref_h[f])[:LMAX] for f in range(NF)])


def two_sided(pmf, O):
    cdf = np.cumsum(pmf)
    lo = cdf[O]
    hi = 1.0 - (cdf[O - 1] if O > 0 else 0.0)
    return min(1.0, 2 * min(lo, hi)), lo, hi


# ---------------- Monte Carlo of the null ----------------
lam_cum = np.cumsum(P)                     # Lambda(x) = sum_{q <= x} P(q)
idx_large = np.flatnonzero(~small)
cumP_large = np.cumsum(P[idx_large]); LAM_L = cumP_large[-1]
idx_small = np.flatnonzero(small)


def simulate(nsim, Pvec, rng):
    """Return list of sorted int arrays of prime indices (exponent sets)."""
    sets = []
    Psm = Pvec[idx_small]
    cl = np.cumsum(Pvec[idx_large]); LL = cl[-1]
    CHS = 2000
    for s in range(0, nsim, CHS):
        n = min(CHS, nsim - s)
        U = rng.random((n, len(idx_small)), dtype=np.float32) < Psm[None, :].astype(np.float32)
        Nl = rng.poisson(LL, size=n)
        draws = np.searchsorted(cl, rng.random(Nl.sum()) * LL)
        draws = idx_large[np.minimum(draws, len(idx_large) - 1)]
        off = np.concatenate([[0], np.cumsum(Nl)])
        for j in range(n):
            a = idx_small[U[j]]
            b = np.unique(draws[off[j]:off[j+1]])
            sets.append(np.concatenate([a, b]))
    return sets


# continuous-statistic ingredients
fr_all = np.log2(p_all.astype(np.float64)) % 1.0
order = np.argsort(fr_all)
fr_sorted = fr_all[order]; Fw = np.cumsum(P[order]) / P.sum()
del order


def F0_frac(t):
    k = np.searchsorted(fr_sorted, t, side="right")
    return np.where(k > 0, Fw[np.maximum(k - 1, 0)], 0.0)


def ks_frac(ix):
    x = np.sort(fr_all[ix]); n = len(x)
    if n == 0:
        return 0.0
    F = F0_frac(x)
    # left limits of F0 at x
    kl = np.searchsorted(fr_sorted, x, side="left")
    Fl = np.where(kl > 0, Fw[np.maximum(kl - 1, 0)], 0.0)
    return max(np.max(np.arange(1, n + 1) / n - F), np.max(Fl - np.arange(n) / n))


bits = (np.floor(np.log2(p_all.astype(np.float64))) + 1).astype(np.int8)
pc = popcount(p_all).astype(np.int8); lr = longest_run(p_all)


def zscore(v):
    z = np.zeros(len(v), dtype=np.float32)
    for L in np.unique(bits):
        m = bits == L
        w = P[m]
        vm = v[m].astype(np.float64)
        mu = np.sum(w * vm) / w.sum()
        sd = math.sqrt(np.sum(w * (vm - mu) ** 2) / w.sum())
        z[m] = (vm - mu) / sd if sd > 1e-9 else 0.0
    return z


z_pc = zscore(pc); z_lr = zscore(lr)
octave = bits  # p in [2^(L-1), 2^L)


def gap_stats(ix):
    u = lam_cum[ix] - P[ix] / 2           # transformed position (mid-step)
    g = np.diff(np.sort(u))
    if len(g) < 3:
        return 0.0, 1.0
    gs = np.sort(g); n = len(gs)
    Fe = 1 - np.exp(-gs)
    ks = max(np.max(np.arange(1, n + 1) / n - Fe), np.max(Fe - np.arange(n) / n))
    cv = gs.std() / gs.mean()
    return ks, cv


E_oct = np.bincount(octave, weights=P, minlength=40)


def disp(ix):
    o = np.bincount(octave[ix], minlength=40)
    m = E_oct > 0
    return float(np.sum((o[m] - E_oct[m]) ** 2 / E_oct[m]))


def cont_stats(ix):
    ix = np.sort(ix)
    big = ix[p_all[ix] >= 1000]
    return np.array([ks_frac(ix), z_pc[ix].sum(), z_lr[ix].sum(),
                     *gap_stats(ix), gap_stats(big)[0], disp(ix)])


CNAMES = ["KS frac(log2 p) vs null", "popcount z-sum", "longest-1-run z-sum",
          "gap KS (Lambda-gaps vs Exp(1)), all", "gap CV (1 = Poisson)",
          "gap KS, exponents >= 1000", "octave-count dispersion chi2"]
# which are two-sided around null mean (z-sums, CV) vs upper-tail (KS, chi2)
TWO = [False, True, True, False, True, False, False]

ix_obs = np.flatnonzero(isexp)
ix_obs_d = ix_obs[p_all[ix_obs] < SPLIT]; ix_obs_h = ix_obs[p_all[ix_obs] >= SPLIT]

sims = simulate(NSIM, P, rng)
S_full = np.array([cont_stats(s) for s in sims])
S_d = np.array([cont_stats(s[p_all[s] < SPLIT]) for s in sims])
S_h = np.array([cont_stats(s[p_all[s] >= SPLIT]) for s in sims])
T_full = cont_stats(ix_obs); T_d = cont_stats(ix_obs_d); T_h = cont_stats(ix_obs_h)


def mc_p(Tobs, S, two):
    if two:
        c = S.mean()
        return (1 + np.sum(np.abs(S - c) >= abs(Tobs - c) - 1e-12)) / (1 + len(S))
    return (1 + np.sum(S >= Tobs - 1e-12)) / (1 + len(S))


K = NF + len(CNAMES)
out(f"B = {B} (GIMPS verified milestone 2026-11-15); primes < B: {nB}; exponents < B: {isexp.sum()}")
out(f"NSIM = {NSIM}; Bonferroni family size K = {K}; threshold 0.05/K = {0.05/K:.2e}")
out()

# ---------------- calibration ----------------
out("== Calibration: per decade of p ==")
out(f"{'range':>22} {'obs':>4} {'E(a=2/6)':>9} {'E(a=2)':>8} {'p(a=2/6)':>9}")
edges = [2, 10, 100, 10**3, 10**4, 10**5, 10**6, 10**7, B]
for lo, hi in zip(edges[:-1], edges[1:]):
    m = (p_all >= lo) & (p_all < hi)
    O = int(isexp[m].sum()); E = P[m].sum(); E2 = P2[m].sum()
    # exact pmf for the decade
    d = np.zeros(LMAX); d[0] = 1
    for q in P[m & small]:
        d = d * (1 - q) + np.concatenate([[0], d[:-1]]) * q
    d = np.convolve(d, poisson_pmf([P[m & ~small].sum()])[0])[:LMAX]
    pv = two_sided(d, O)[0]
    out(f"{f'[{lo}, {hi})':>22} {O:>4} {E:>9.2f} {E2:>8.2f} {pv:>9.3f}")
for lab, m in [("all p < B", np.ones(nB, bool)), ("p < 1e5", small), ("1e5 <= p < B", ~small)]:
    out(f"{lab:>22} {int(isexp[m].sum()):>4} {P[m].sum():>9.2f} {P2[m].sum():>8.2f}")
out()
out("== mod 4 split ==")
for lab, m in [("all", np.ones(nB, bool)), ("p<1e5", small), ("p>=1e5", ~small)]:
    r = []
    for res in (1, 3):
        mm = m & (p_all % 4 == res)
        r.append((int(isexp[mm].sum()), P[mm].sum(), P2[mm].sum()))
    out(f"{lab:>6}: p=1 mod 4: obs {r[0][0]}, E(a=6) {r[0][1]:.2f}, E(a=2) {r[0][2]:.2f} | "
        f"p=3 mod 4: obs {r[1][0]}, E {r[1][1]:.2f} (same both models)")


def class_p(mask, Pv):
    d = np.zeros(LMAX); d[0] = 1
    for q in Pv[mask & small]:
        d = d * (1 - q) + np.concatenate([[0], d[:-1]]) * q
    d = np.convolve(d, poisson_pmf([Pv[mask & ~small].sum()])[0])[:LMAX]
    O = int(isexp[mask].sum())
    return O, float(Pv[mask].sum()), two_sided(d, O)[0]


for res_ in (1, 3):
    mm = (p_all % 4 == res_)
    for lab, Pv in (("a=2/6", P), ("a=2", P2), ("refined", Pref)):
        O, E, pv = class_p(mm, Pv)
        out(f"  p={res_} mod 4 under {lab:>7}: obs {O}, E {E:.2f}, two-sided p {pv:.3f}")
for lab, Pv in (("a=2/6", P), ("a=2", P2), ("refined", Pref)):
    O, E, pv = class_p((p_all >= 10**7), Pv)
    out(f"  p in [1e7, B) under {lab:>7}: obs {O}, E {E:.2f}, two-sided p {pv:.4f}")
    O, E, pv = class_p((p_all >= SPLIT), Pv)
    out(f"  p in [1e5, B) under {lab:>7}: obs {O}, E {E:.2f}, two-sided p {pv:.4f}")
out(f"  refined null: {int(killed.sum())} Euler-killed primes, removed mass {P[killed].sum():.2f}; total E {Pref.sum():.2f}")


def loglik(Pv):
    Pv = np.clip(Pv, 1e-300, 1 - 1e-15)
    return float(np.sum(np.log(Pv[isexp])) + np.sum(np.log1p(-Pv[~isexp])))


dLL = loglik(P) - loglik(P2)
out(f"log-likelihood(a=2/6) - log-likelihood(a=2) = {dLL:+.3f}  (Bayes factor {math.exp(dLL):.2f})")
# power: distribution of dLL under each model
def _lterm(Pv):
    Pv = np.clip(Pv, 1e-300, 1 - 1e-15)
    return np.log(Pv) - np.log1p(-Pv)


base = float(np.sum(np.log1p(-np.minimum(P, 1 - 1e-15))) - np.sum(np.log1p(-np.minimum(P2, 1 - 1e-15))))


def dll_set(ix):
    return base + float(np.sum(_lterm(P[ix]) - _lterm(P2[ix])))


d6 = np.array([dll_set(s) for s in sims[:5000]])
sims2 = simulate(5000, P2, np.random.default_rng(7))
d2 = np.array([dll_set(s) for s in sims2])
out(f"  dLL under a=2/6 model: mean {d6.mean():+.2f} sd {d6.std():.2f}; under a=2 model: mean {d2.mean():+.2f} sd {d2.std():.2f}")
out(f"  P(dLL >= obs | a=2 model) = {np.mean(d2 >= dLL):.3f};  P(dLL <= obs | a=2/6 model) = {np.mean(d6 <= dLL):.3f}")
out()

# ---------------- binary features ----------------
res = []
for f in range(NF):
    pf = two_sided(pmf_full[f], int(O_full[f]))[0]
    pd_, lo_d, hi_d = two_sided(pmf_d[f], int(O_d[f]))
    ph, lo_h, hi_h = two_sided(pmf_h[f], int(O_h[f]))
    pr_ = two_sided(pmf_ref[f], int(O_full[f]))[0]
    res.append(dict(name=names[f], O=int(O_full[f]), E=E_full[f], E2=E2_full[f], p=pf, Eref=Eref[f], pref=pr_,
                    padj=min(1, K * pf), Od=int(O_d[f]), Ed=E_d[f], pd=pd_,
                    Oh=int(O_h[f]), Eh=E_h[f], ph=ph,
                    ph_dir=(hi_h if O_d[f] > E_d[f] else lo_h)))
out("== Binary features (two-sided exact Poisson-binomial; padj = Bonferroni K*p) ==")
out(f"{'feature':>42} {'O':>3} {'E':>6} {'p':>7} {'padj':>5} | {'Od':>3} {'Ed':>6} {'pd':>6} | {'Oh':>3} {'Eh':>6} {'ph':>6} | {'Eref':>6} {'pref':>6}")
for r in res:
    out(f"{r['name']:>42} {r['O']:>3} {r['E']:>6.2f} {r['p']:>7.4f} {r['padj']:>5.2f} | "
        f"{r['Od']:>3} {r['Ed']:>6.2f} {r['pd']:>6.3f} | {r['Oh']:>3} {r['Eh']:>6.2f} {r['ph']:>6.3f} | {r['Eref']:>6.2f} {r['pref']:>6.3f}")
out()

# omnibus per modulus (chi-square vs null, MC p-value)
out("== Omnibus per modulus (chi2 over coprime classes, MC p) ==")
col = 0
FEs = None
for m in MODS:
    ncl = sum(1 for r in range(m) if math.gcd(r, m) == 1)
    cols = list(range(col, col + ncl)); col += ncl
    Ec = E_full[cols]
    chi_obs = float(np.sum((O_full[cols] - Ec) ** 2 / Ec))
    cop = np.array([r for r in range(m) if math.gcd(r, m) == 1])
    chis = np.array([np.sum((np.bincount(p_all[s] % m, minlength=m)[cop] - Ec) ** 2 / Ec) for s in sims])
    out(f"  mod {m:>3}: chi2 = {chi_obs:6.2f} (df {ncl-1}), MC p = {(1+np.sum(chis>=chi_obs))/(1+len(chis)):.3f}")
out()

out("== Continuous statistics (MC p-values; full / discovery p<1e5 / holdout p>=1e5) ==")
cres = []
for k, nm in enumerate(CNAMES):
    pf = mc_p(T_full[k], S_full[:, k], TWO[k]); pd_ = mc_p(T_d[k], S_d[:, k], TWO[k]); ph = mc_p(T_h[k], S_h[:, k], TWO[k])
    cres.append(dict(name=nm, T=T_full[k], Tnull=S_full[:, k].mean(), p=pf, padj=min(1, K*pf), Td=T_d[k], pd=pd_, Th=T_h[k], ph=ph))
    plo = (1 + np.sum(S_full[:, k] <= T_full[k] + 1e-12)) / (1 + len(S_full))
    out(f"{nm:>40}: obs {T_full[k]:7.3f} null mean {S_full[:, k].mean():7.3f} p {pf:.4f} (lower tail {plo:.4f}) padj {min(1,K*pf):.2f} | "
        f"disc obs {T_d[k]:7.3f} p {pd_:.3f} | hold obs {T_h[k]:7.3f} (null {S_h[:, k].mean():.3f}) p {ph:.3f}")
# asymptotic check for gap KS on exponents >= 1000 and >= 1e5
from scipy.stats import kstest
for thr in (1000, SPLIT):
    ix = ix_obs[p_all[ix_obs] >= thr]
    u = lam_cum[ix] - P[ix] / 2
    g = np.diff(u)
    kt = kstest(g, "expon")
    out(f"  asymptotic KS of Lambda-gaps vs Exp(1), exponents >= {thr}: n={len(g)}, D={kt.statistic:.3f}, p={kt.pvalue:.3f}; "
        f"mean gap {g.mean():.3f}, CV {g.std()/g.mean():.3f}")
out()

# ---------------- out-of-sample protocol ----------------
out("== Out-of-sample: features with discovery p < 0.05 (p < 1e5), tested one-sided in the same direction on p >= 1e5 ==")
sel = [r for r in res if r["pd"] < 0.05] + [c for c in cres if c["pd"] < 0.05]
for r in sel:
    if "Od" in r:
        out(f"  {r['name']}: disc O={r['Od']} E={r['Ed']:.2f} pd={r['pd']:.4f} -> holdout O={r['Oh']} E={r['Eh']:.2f} one-sided p={r['ph_dir']:.3f}")
    else:
        out(f"  {r['name']}: disc obs={r['Td']:.3f} pd={r['pd']:.4f} -> holdout obs={r['Th']:.3f} two-sided p={r['ph']:.3f}")
out(f"  ({len(sel)} features carried forward; Bonferroni threshold on holdout {0.05/max(1,len(sel)):.3f})")

# family-wise min-p over binary features via MC (lookup in exact pmfs)
cdfs = np.cumsum(pmf_full, axis=1)
minp = []
NFW = min(NSIM, 5000)
flat = np.concatenate(sims[:NFW]); lens = np.array([len(s) for s in sims[:NFW]])
sid = np.repeat(np.arange(NFW), lens)
FM = feature_matrix(p_all[flat])
Fsim_counts = np.zeros((NFW, NF), dtype=np.int64)
for f in range(NF):
    Fsim_counts[:, f] = np.bincount(sid, weights=FM[:, f], minlength=NFW).astype(np.int64)
for j in range(len(Fsim_counts)):
    c = Fsim_counts[j]
    lo = cdfs[np.arange(NF), c]
    hi = 1 - np.where(c > 0, cdfs[np.arange(NF), np.maximum(c - 1, 0)], 0.0)
    minp.append(np.min(np.minimum(1, 2 * np.minimum(lo, hi))))
minp = np.array(minp)
obs_min = min(r["p"] for r in res)
out(f"\nFamily-wise (MC, 5000 null sets): smallest binary-feature p observed = {obs_min:.4f}; "
    f"fraction of null sets with a smaller min-p = {np.mean(minp <= obs_min):.3f}")

with open("results.txt", "w") as fh:
    fh.write("\n".join(log) + "\n")
json.dump(dict(binary=res, continuous=[{k: float(v) if isinstance(v, (np.floating, float)) else v for k, v in c.items()} for c in cres]),
          open("results.json", "w"), indent=1, default=float)
