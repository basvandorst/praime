"""Cross-class HL covariances: eps(c) = Cov(X_a, X_b)/Var X_a for c = b - a, sampled; noise inflation of the spread.
Var(SS) = 2 sigma^4 [ (K-1) + 2 tr(PE) + tr(PEPE) ],  tr(PE) = -(K-1) eps_bar, tr(PEPE) ~ K(K-1) var(eps)."""
import sys, math, json, random
import numpy as np
import theory as T

Q = int(sys.argv[1]); L0 = int(float(sys.argv[2])); W = int(float(sys.argv[3])); NS = int(sys.argv[4])
S = T.setup(Q, L0, W)
r = S['r']; mu = float(r.sum())
F0, _, _ = T.F(S, 0)
var = float(np.sum(r * (1 - r))) + F0
adm = [a for a in range(Q) if math.gcd(a * (a + 2), Q) == 1]
K = len(adm)
diffs = sorted({(b - a) % Q for a in adm for b in adm if a != b})
cnt = {}
for a in adm:
    for b in adm:
        if a != b:
            cnt[(b - a) % Q] = cnt.get((b - a) % Q, 0) + 1
random.seed(1)
samp = random.sample(diffs, min(NS, len(diffs)))
eps = []
for c in samp:
    Fc, _, _ = T.F(S, c)
    eps.append(Fc / var)
eps = np.array(eps); wts = np.array([cnt[c] for c in samp], float)
eb = float(np.sum(wts * eps) / wts.sum()); ve = float(np.sum(wts * (eps - eb) ** 2) / wts.sum())
infl = 1 - 2 * eb + K * ve
out = dict(Q=Q, L0=L0, W=W, K=K, v=var / mu, n_diffs=len(diffs), n_sampled=len(samp), eps_mean=eb, eps_sd=math.sqrt(ve),
           eps_min=float(eps.min()), eps_max=float(eps.max()), var_inflation=infl, sd_inflation=math.sqrt(infl),
           implied_vtot_over_v=1 + (K - 1) * eb)
print(json.dumps(out, indent=1))
json.dump(dict(out, samp=samp, eps=eps.tolist()), open(f'cross_{Q}_{L0}.json', 'w'))
