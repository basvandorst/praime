"""Analysis of the streamed LL / PRP statistics (runs_100_12000.pkl).
Pre-registered pipeline: see CRITERIA.txt.  Output: analysis.txt
"""
import sys, math, pickle, time
import numpy as np
from scipy.stats import rankdata, norm

EXP = {2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281,
       3217, 4253, 4423, 9689, 9941, 11213, 19937}
HOR = ['j32', 'j64', 'j128', 'j256', 'k/16', 'k/8', 'k/4', 'k/2']
SKIP = {'i0', 'a_n', 'chk_jac_s+2', 'chk_jac_s-2', 'v_dist_pm2m', 'a_hw_maxabs'}
NPERM = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
rng = np.random.default_rng(20261115)
out = open('analysis.txt', 'w')


def P(*a):
    s = ' '.join(str(x) for x in a)
    print(s); out.write(s + '\n'); out.flush()


R = pickle.load(open(sys.argv[1] if len(sys.argv) > 1 else 'runs_100_12000.pkl', 'rb'))
ks = np.array([r['k'] for r in R])
y = np.array([1 if k in EXP else 0 for k in ks])
tr = (ks >= 100) & (ks < 3000)
te = (ks > 3000) & (ks < 12000)
P('survivors: train', tr.sum(), 'pos', y[tr].sum(), '| test', te.sum(), 'pos', y[te].sum())


def snap(r, seq, h):
    """horizon j is capped at k//2: for k < 2j the fixed horizon is not reached
    before the end of the run, so the k//2 snapshot is used (k < 258 only)."""
    sn = r[seq]['snaps']
    return sn[h] if h in sn else sn['k/2']


def features(h):
    names, cols = [], []
    for seq in ('LL', 'PRP'):
        keys = [x for x in snap(R[0], seq, h).keys() if x not in SKIP]
        for key in keys:
            names.append(seq + ':' + key)
            cols.append([float(snap(r, seq, h)[key]) for r in R])
        # max |hw z| over n iterations -> its null CDF (k-free): (2Phi(m)-1)^n
        names.append(seq + ':a_hw_maxabs_u')
        col = []
        for r in R:
            f = snap(r, seq, h)
            n = max(f['a_n'], 1)
            col.append((2 * norm.cdf(f['a_hw_maxabs']) - 1) ** n if f['a_n'] else 0.5)
        cols.append(col)
        names.append(seq + ':v_absz')
        cols.append([abs(snap(r, seq, h)['v_hwz']) for r in R])
    return names, np.array(cols).T


def auc_from_ranks(rk, pos):
    """rk: ranks (n x F) within a split; pos: boolean (n,) -> AUC per feature."""
    n1 = pos.sum(); n0 = len(pos) - n1
    return (rk[pos].sum(0) - n1 * (n1 + 1) / 2) / (n1 * n0)


def auc(score, lab):
    rk = rankdata(score)
    n1 = lab.sum(); n0 = len(lab) - n1
    return (rk[lab == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


# ---------------- models ----------------
def single(Xtr, ytr, Xte, yte, rktr=None):
    if rktr is None:
        rktr = np.apply_along_axis(rankdata, 0, Xtr)
    a = auc_from_ranks(rktr, ytr == 1)
    f = int(np.argmax(np.abs(a - 0.5)))
    sg = 1.0 if a[f] >= 0.5 else -1.0
    s_tr = sg * Xtr[:, f]; s_te = sg * Xte[:, f]
    return f, sg, a[f] if sg > 0 else 1 - a[f], s_tr, s_te, a


def logreg(Xtr, ytr, Xte, a_tr, nf=3, lam=1.0):
    idx = np.argsort(-np.abs(a_tr - 0.5))[:nf]
    mu = Xtr[:, idx].mean(0); sd = Xtr[:, idx].std(0) + 1e-12
    A = np.c_[np.ones(len(Xtr)), (Xtr[:, idx] - mu) / sd]
    B = np.c_[np.ones(len(Xte)), (Xte[:, idx] - mu) / sd]
    w = np.zeros(A.shape[1])
    reg = lam * np.eye(A.shape[1]); reg[0, 0] = 0
    for _ in range(30):
        p = 1 / (1 + np.exp(-A @ w))
        g = A.T @ (p - ytr) + reg @ w
        Hs = (A * (p * (1 - p))[:, None]).T @ A + reg
        step = np.linalg.solve(Hs, g)
        w -= step
        if np.abs(step).max() < 1e-9:
            break
    return A @ w, B @ w, idx


def gini_split(X, yy, feats, minleaf=5):
    best = (1e9, None, None)
    n = len(yy)
    for f in feats:
        ts = np.unique(np.quantile(X[:, f], np.linspace(0.05, 0.95, 19)))
        L = X[:, f][:, None] <= ts[None, :]
        nl = L.sum(0); nr = n - nl
        pl = (L & (yy[:, None] == 1)).sum(0); pr = yy.sum() - pl
        ok = (nl >= minleaf) & (nr >= minleaf)
        with np.errstate(invalid='ignore', divide='ignore'):
            gl = 2 * pl / nl * (1 - pl / nl); gr = 2 * pr / nr * (1 - pr / nr)
            G = (nl * gl + nr * gr) / n
        G[~ok] = 1e9
        i = int(np.argmin(G))
        if G[i] < best[0]:
            best = (G[i], f, ts[i])
    return best


def tree2(Xtr, ytr, Xte, a_tr, nf=10):
    feats = np.argsort(-np.abs(a_tr - 0.5))[:nf]
    _, f0, t0 = gini_split(Xtr, ytr, feats)
    if f0 is None:
        sc = (ytr.sum() + 1) / (len(ytr) + 2)
        return np.full(len(Xtr), sc), np.full(len(Xte), sc)
    str_ = np.zeros(len(Xtr)); ste = np.zeros(len(Xte))
    for side in (0, 1):
        mtr = (Xtr[:, f0] <= t0) if side == 0 else (Xtr[:, f0] > t0)
        mte = (Xte[:, f0] <= t0) if side == 0 else (Xte[:, f0] > t0)
        Xs, ys = Xtr[mtr], ytr[mtr]
        sub = (None, None, None)
        if 0 < ys.sum() < len(ys):
            sub = gini_split(Xs, ys, feats)
        if sub[1] is None:
            v = (ys.sum() + 1) / (len(ys) + 2)
            str_[mtr] = v; ste[mte] = v
        else:
            f1, t1 = sub[1], sub[2]
            for s2 in (0, 1):
                m2 = (Xs[:, f1] <= t1) if s2 == 0 else (Xs[:, f1] > t1)
                v = (ys[m2].sum() + 1) / (m2.sum() + 2)
                idx_tr = np.nonzero(mtr)[0][m2]
                str_[idx_tr] = v
                Xe = Xte[mte]
                m2e = (Xe[:, f1] <= t1) if s2 == 0 else (Xe[:, f1] > t1)
                ste[np.nonzero(mte)[0][m2e]] = v
    return str_, ste


def evaluate(s_tr, ytr, s_te, yte):
    thr = s_tr[ytr == 1].min()
    rec = (s_te[yte == 1] >= thr).mean()
    rej = (s_te[yte == 0] < thr).mean()
    a = auc(s_te, yte)
    orc = (s_te[yte == 0] < s_te[yte == 1].min()).mean()   # oracle: rejection at full test recall
    return a, rec, rej, orc


def pipeline(X, ytr_, yte_, rktr=None):
    Xtr, Xte = X[tr], X[te]
    f, sg, atr, s_tr, s_te, a_all = single(Xtr, ytr_, Xte, yte_, rktr)
    res = {'single': (f, sg, atr) + evaluate(s_tr, ytr_, s_te, yte_)}
    l_tr, l_te, idx = logreg(Xtr, ytr_, Xte, a_all)
    res['logreg'] = (tuple(idx), 0, auc(l_tr, ytr_)) + evaluate(l_tr, ytr_, l_te, yte_)
    t_tr, t_te = tree2(Xtr, ytr_, Xte, a_all)
    res['tree2'] = (None, 0, auc(t_tr, ytr_)) + evaluate(t_tr, ytr_, t_te, yte_)
    return res


# ---------- null label generators ----------
idx_tr = np.nonzero(tr)[0]; idx_te = np.nonzero(te)[0]


def perm_global():
    a = np.zeros(len(idx_tr), int); a[rng.choice(len(idx_tr), y[tr].sum(), replace=False)] = 1
    b = np.zeros(len(idx_te), int); b[rng.choice(len(idx_te), y[te].sum(), replace=False)] = 1
    return a, b


def local_choice(idx, yy, m=30):
    kk = ks[idx]
    lab = np.zeros(len(idx), int)
    for p in np.nonzero(yy == 1)[0]:
        near = np.argsort(np.abs(kk - kk[p]), kind='stable')[:m]
        near = [q for q in near if lab[q] == 0]
        lab[rng.choice(near)] = 1
    return lab


def perm_local():
    return local_choice(idx_tr, y[tr]), local_choice(idx_te, y[te])


# ---------- run ----------
t0 = time.time()
Xs = {}
for h in HOR:
    names, X = features(h)
    Xs[h] = X
P('features per horizon:', len(names))
P()
P('controls (not in the selection pool): train/test AUC of k itself (exponents are denser at small k: Wagstaff)')
for nm, v in (('k', ks.astype(float)), ('-k', -ks.astype(float))):
    P('  %-4s train AUC %.3f  test AUC %.3f' % (nm, auc(v[tr], y[tr]), auc(v[te], y[te])))
P()

obs = {}
for h in HOR:
    obs[h] = pipeline(Xs[h], y[tr], y[te])
MODELS = ['single', 'logreg', 'tree2']
nulls = {nm: {(h, m): [] for h in HOR for m in MODELS} for nm in ('global', 'klocal')}
for nm, gen in (('global', perm_global), ('klocal', perm_local)):
    for it in range(NPERM):
        a, b = gen()
        for h in HOR:
            r = pipeline(Xs[h], a, b)
            for m in MODELS:
                nulls[nm][(h, m)].append(r[m][3:])
    P('null %s: %d permutations done, %.0f s' % (nm, NPERM, time.time() - t0))
for nm in nulls:
    for key in nulls[nm]:
        nulls[nm][key] = np.array(nulls[nm][key])   # cols: testAUC, recall, rej, oracle


def pval(arr, x):
    return (1 + (arr >= x - 1e-12).sum()) / (1 + len(arr))


P()
P('columns: horizon model | chosen feature (train AUC) | test AUC | test recall at train thr | test rejection at train thr | oracle test rejection at full recall | p(test AUC) global, k-local | p(rejection with full recall) global, k-local')
fam = {nm: [] for nm in nulls}
for h in HOR:
    for m in MODELS:
        f, sg, atr, a, rec, rej, orc = obs[h][m]
        if m == 'single':
            lab = ('+' if sg > 0 else '-') + names[f]
        elif m == 'logreg':
            lab = 'LR[' + ','.join(names[i] for i in f) + ']'
        else:
            lab = 'tree2'
        ps = []
        for nm in ('global', 'klocal'):
            N = nulls[nm][(h, m)]
            pa = pval(N[:, 0], a)
            stat_obs = rej if rec == 1 else 0.0
            stat_null = np.where(N[:, 1] == 1, N[:, 2], 0.0)
            pr = pval(stat_null, stat_obs)
            ps += [pa, pr]
        P('%-5s %-6s | %s (%.3f) | %.3f | %d/6 | %.3f | %.3f | %.3f %.3f | %.3f %.3f' % (
            h, m, lab, atr, a, round(rec * 6), rej, orc, ps[0], ps[2], ps[1], ps[3]))

# family-wise: min over the 24 (horizon, model) p-values of the test AUC
P()
for nm in nulls:
    keys = [(h, m) for h in HOR for m in MODELS]
    A = np.stack([nulls[nm][key][:, 0] for key in keys], 1)       # NPERM x 24
    # per-permutation p-values (rank within each column)
    rk = (A[None, :, :] >= A[:, None, :] - 1e-12).sum(1)  # rk[a,c] = #{b: A[b,c] >= A[a,c]}
    pmat = (1 + rk) / (1 + NPERM)
    minp_null = pmat.min(1)
    obs_p = min(pval(nulls[nm][key][:, 0], obs[key[0]][key[1]][3]) for key in keys)
    fw = (1 + (minp_null <= obs_p + 1e-12).sum()) / (1 + NPERM)
    P('family-wise (%s null): smallest observed p(test AUC) = %.3f, family-wise p = %.3f' % (nm, obs_p, fw))
    # same for the statistic 'test rejection if test recall = 6/6, else 0'
    A = np.stack([np.where(nulls[nm][key][:, 1] == 1, nulls[nm][key][:, 2], 0.0) for key in keys], 1)
    rk = (A[None, :, :] >= A[:, None, :] - 1e-12).sum(1)
    minp_null = ((1 + rk) / (1 + NPERM)).min(1)
    ob = [obs[key[0]][key[1]][5] if obs[key[0]][key[1]][4] == 1 else 0.0 for key in keys]
    obs_p = min(pval(A[:, c], ob[c]) for c in range(len(keys)))
    fw = (1 + (minp_null <= obs_p + 1e-12).sum()) / (1 + NPERM)
    P('family-wise (%s null): smallest observed p(rejection at full recall) = %.3f, family-wise p = %.3f' % (nm, obs_p, fw))

# post-hoc: best single feature on the TEST set itself (cherry-picked), against the same max under the null
P()
P('post-hoc (cherry-picked on the test set): max |test AUC - 0.5| over all features, observed vs global null')
for h in HOR:
    X = Xs[h][te]
    rk = np.apply_along_axis(rankdata, 0, X)
    aa = auc_from_ranks(rk, y[te] == 1)
    best = np.argmax(np.abs(aa - 0.5)); ob = abs(aa[best] - 0.5)
    nul = []
    for it in range(NPERM):
        b = np.zeros(len(idx_te), bool); b[rng.choice(len(idx_te), y[te].sum(), replace=False)] = True
        nul.append(np.abs(auc_from_ranks(rk, b) - 0.5).max())
    nul = np.array(nul)
    P('  %-5s %-28s test AUC %.3f  null 95%% of max|AUC-0.5| = %.3f  p = %.3f' % (
        h, names[best], aa[best], np.quantile(nul, 0.95), pval(nul, ob)))
P('time %.0f s' % (time.time() - t0))
pickle.dump({'obs': obs, 'names': names, 'nulls': nulls}, open('analysis.pkl', 'wb'))
