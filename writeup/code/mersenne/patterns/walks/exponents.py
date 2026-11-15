"""Where the 31 Mersenne exponents <= 500000 sit on the walks, the turtle path and the rows,
against null sets of pseudo-exponents: each prime p <= N taken with probability
w(p) = min(1, e^gamma log2(a p) / p), a = 6 if p = 1 mod 4, a = 2 if p = 3 mod 4 (and p = 2),
conditioned on exactly 31 picks (exact conditional-Bernoulli sampling).
Statistics fixed before the run; two-sided rank p-values; family-wise p by the min-p method."""
import json, math
import numpy as np
from common import series, N, EXPONENTS

S = series(N)
primes = S["primes"]; M, L, R, X, Y, pi = S["M"], S["L"], S["R"], S["X"], S["Y"], S["pi"]
P = primes.astype(float)
a = np.where(primes % 4 == 1, 6.0, 2.0)
w = np.minimum(1.0, math.exp(0.5772156649015329) * np.log2(a * P) / P)
K = 31
NS = 3000
rng = np.random.default_rng(31)

# Q[i, r] = P(sum_{j>=i} B_j = r), r = 0..K
m = primes.size
Q = np.zeros((m + 1, K + 1)); Q[m, 0] = 1.0
for i in range(m - 1, -1, -1):
    Q[i] = (1 - w[i]) * Q[i + 1]
    Q[i, 1:] += w[i] * Q[i + 1, :-1]
need = np.full(NS, K)
picks = np.zeros((NS, K), dtype=np.int64)
cnt = np.zeros(NS, dtype=np.int64)
for i in range(m):
    act = need > 0
    if not act.any():
        break
    r = need[act]
    num = w[i] * Q[i + 1, r - 1]
    den = Q[i, r]
    pr = np.where(den > 0, num / np.where(den > 0, den, 1), 0.0)
    take = rng.random(r.size) < pr
    idx = np.nonzero(act)[0][take]
    picks[idx, cnt[idx]] = primes[i]
    cnt[idx] += 1
    need[idx] -= 1
assert np.all(cnt == K)

E = np.array(EXPONENTS)
sq = np.sqrt

names = ["T1 mean M(p)/sqrt p", "T2 mean |M(p)|/sqrt p", "T3 mean L(p)/sqrt p",
         "T4 mean |L(p)/sqrt p - 1/zeta(1/2)|", "T5 mean R(p) log p/sqrt p", "T6 #{R(p) <= 0}",
         "T7 mean |Z(p)|/sqrt(p log p)", "T8 chi2 of pi(p) mod 4", "T9 mean X(p)/sqrt(p log p)",
         "T10 mean Y(p)/sqrt(p log p)", "T11 mean omega/sieverows x = pi(p)/pi(N)", "T12 mean ruler x ~ p/N"]
one_sided_up = {7}
def T(P2):
    """P2: (S, K) array of primes -> (S, 12) statistics"""
    f = P2.astype(float); lg = np.log(f); sc = np.sqrt(f * lg)
    out = np.empty((P2.shape[0], 12))
    out[:, 0] = np.mean(M[P2] / sq(f), 1)
    out[:, 1] = np.mean(np.abs(M[P2]) / sq(f), 1)
    out[:, 2] = np.mean(L[P2] / sq(f), 1)
    out[:, 3] = np.mean(np.abs(L[P2] / sq(f) + 0.6847652800500), 1)
    out[:, 4] = np.mean(R[P2] * lg / sq(f), 1)
    out[:, 5] = np.sum(R[P2] <= 0, 1)
    out[:, 6] = np.mean(np.hypot(X[P2], Y[P2]) / sc, 1)
    c = np.stack([np.sum(pi[P2] % 4 == j, 1) for j in range(4)], 1)
    out[:, 7] = np.sum((c - K / 4) ** 2 / (K / 4), 1)
    out[:, 8] = np.mean(X[P2] / sc, 1)
    out[:, 9] = np.mean(Y[P2] / sc, 1)
    out[:, 10] = np.mean(pi[P2] / pi[N], 1)
    out[:, 11] = np.mean(P2 / N, 1)
    return out

t_obs = T(E[None, :])[0]
t_null = T(picks)
def pval(t, null, j):
    lo = (1 + np.sum(null <= t)) / (null.size + 1)
    hi = (1 + np.sum(null >= t)) / (null.size + 1)
    return hi if j in one_sided_up else min(1.0, 2 * min(lo, hi))
p_obs = np.array([pval(t_obs[j], t_null[:, j], j) for j in range(12)])
# p-values of each null set (against the others) for the min-p family-wise correction
ranks_p = np.empty_like(t_null)
for j in range(12):
    col = t_null[:, j]
    srt = np.sort(col)
    le = np.searchsorted(srt, col, side="right")
    ge = col.size - np.searchsorted(srt, col, side="left")
    lo = le / (col.size + 1); hi = ge / (col.size + 1)
    ranks_p[:, j] = hi if j in one_sided_up else np.minimum(1.0, 2 * np.minimum(lo, hi))
minp_null = ranks_p.min(1)
fw = (1 + np.sum(minp_null <= p_obs.min())) / (NS + 1)
res = {"n_null_sets": NS, "sum_w": float(w.sum()), "statistics": []}
for j in range(12):
    res["statistics"].append(dict(name=names[j], observed=float(t_obs[j]), null_mean=float(t_null[:, j].mean()),
                                  null_sd=float(t_null[:, j].std()), p=float(p_obs[j])))
res["min_p"] = float(p_obs.min()); res["family_wise_p"] = float(fw)
res["exponent_values"] = {int(p): dict(M=int(M[p]), L=int(L[p]), R=int(R[p]), X=int(X[p]), Y=int(Y[p]),
                                       pi_mod4=int(pi[p] % 4)) for p in E}
res["exponent_heading_classes"] = np.bincount(pi[E] % 4, minlength=4).tolist()
json.dump(res, open("exponents_out.json", "w"), indent=1)
print("sum w =", res["sum_w"])
for s in res["statistics"]:
    print(f"{s['name']:45s} obs {s['observed']:8.4f}  null {s['null_mean']:8.4f} +- {s['null_sd']:.4f}  p {s['p']:.4f}")
print("min p", res["min_p"], "family-wise p", fw)
print("heading classes of exponents", res["exponent_heading_classes"])
