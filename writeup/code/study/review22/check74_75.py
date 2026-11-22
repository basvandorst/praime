"""Rounds 72, 74, 75 check, own code (the theory table of round 72 is read from the .npy only).
usage: check74_75.py p1e9.npy table.npy [p1e12.npy]"""
import sys, math, json
import numpy as np

P = np.load(sys.argv[1]); TH = np.load(sys.argv[2])
Q = 30030; QS = (3, 5, 7, 11, 13); KMAX = 150
CLS = np.array([a for a in range(Q) if math.gcd(a, Q) == 1]); IDX = np.full(Q, -1); IDX[CLS] = np.arange(len(CLS))
COP = np.zeros(Q, bool); COP[CLS] = True
out = {}
rs = TH.sum(1); out['theory_row_sum_minmax'] = [float(rs.min()), float(rs.max())]
TH = TH / rs[:, None]; NG = TH.shape[1]; gaps_th = np.arange(2, 2 * NG + 1, 2)
lam = math.log(1.05e9)

# ---------- empirical tables (backward: gap before p by class of p; forward: gap after p by class of p)
g = np.diff(P)
back_cls = IDX[P[1:] % Q]; fwd_cls = IDX[P[:-1] % Q]
GM = int(g.max()); NGE = GM // 2
Eb = np.zeros((len(CLS), NGE)); np.add.at(Eb, (back_cls, g // 2 - 1), 1)
Ef = np.zeros((len(CLS), NGE)); np.add.at(Ef, (fwd_cls, g // 2 - 1), 1)
out['min_per_class_back'] = int(Eb.sum(1).min())
# reflection in the data: backward table at c vs forward table at -c, chi^2 of homogeneity on mod-q histograms
refl = {}
negidx = IDX[(-CLS) % Q]
for q in QS:
    gm = (np.arange(2, 2 * NGE + 1, 2) % q)
    hb = np.stack([Eb[:, gm == r].sum(1) for r in range(q)], 1); hf = np.stack([Ef[negidx][:, gm == r].sum(1) for r in range(q)], 1)
    tot = hb + hf; nb = hb.sum(1, keepdims=True); nf = hf.sum(1, keepdims=True); n = nb + nf
    eb = tot * nb / n; ef = tot * nf / n
    m = tot > 0
    chi = float(np.sum(((hb - eb) ** 2 / np.where(eb > 0, eb, 1))[m] + ((hf - ef) ** 2 / np.where(ef > 0, ef, 1))[m]))
    dof = int(np.sum(m) - len(CLS))
    refl[q] = {'chi2': chi, 'dof': dof, 'z': (chi - dof) / math.sqrt(2 * dof)}
out['reflection_data_back_vs_fwd_neg'] = refl

# theory: class-averaged P(gap = 0 mod q) (uniform classes), against data
out['P_gap0_theory_uniform_classes'] = {q: float(TH[:, gaps_th % q == 0].sum(1).mean()) for q in QS}
out['P_gap0_data'] = {q: float(np.mean(g % q == 0)) for q in QS}

# ---------- chains
order_b = np.argsort(back_cls, kind='stable'); gs_b = g[order_b]; start_b = np.searchsorted(back_cls[order_b], np.arange(len(CLS)))
cnt_b = np.bincount(back_cls, minlength=len(CLS))
cumTH = np.cumsum(TH, 1); cumTH[:, -1] = 1.0; flat = (cumTH + np.arange(len(CLS))[:, None]).ravel()


def sample_theory(cur, rng):
    row = IDX[(-cur) % Q]                       # backward from cur = forward law of -cur
    u = rng.random(len(cur))
    pos = np.searchsorted(flat, row + u); j = pos - row * NG; j = np.clip(j, 0, NG - 1)
    return gaps_th[j]


def make_emp(gl, cl):
    o = np.argsort(cl, kind='stable'); gs = gl[o]; st = np.searchsorted(cl[o], np.arange(len(CLS))); c = np.bincount(cl, minlength=len(CLS))
    assert c.min() > 0

    def samp(cur, rng):
        r = IDX[cur]; return gs[st[r] + (rng.random(len(cur)) * c[r]).astype(np.int64)]
    return samp


def run_chain(sampler, nch, seed, start=None):
    rng = np.random.default_rng(seed); Ts = []; Tqs = {q: [] for q in QS}; rep = {q: 0 for q in QS}; ns = 0
    for b0 in range(0, nch, 20000):
        n = min(20000, nch - b0)
        cur = CLS[rng.integers(0, len(CLS), n)] if start is None else start[b0:b0 + n].copy()
        D = np.zeros(n, np.int64); seen = {q: np.zeros((n, q), bool) for q in QS}; f = {q: np.full(n, KMAX + 1) for q in QS}
        for q in QS:
            seen[q][:, 0] = True
        for k in range(1, KMAX + 1):
            gg = sampler(cur, rng); D += gg; cur = (cur - gg) % Q
            for q in QS:
                rep[q] += int(np.sum(gg % q == 0))
            ns += n
            for q in QS:
                seen[q][np.arange(n), D % q] = True
                f[q][(seen[q].sum(1) == q - 1) & (f[q] > KMAX)] = k
        assert COP[cur].all()
        Ts.append(np.max(np.stack([f[q] for q in QS]), 0))
        for q in QS:
            Tqs[q].append(f[q])
    T = np.concatenate(Ts)
    return {'mean_T': float(T.mean()), 'se': float(T.std() / math.sqrt(len(T))), 'unres': int(np.sum(T > KMAX)),
            **{'T_q%d' % q: float(np.concatenate(Tqs[q]).mean()) for q in QS}, **{'rep_q%d' % q: rep[q] / ns for q in QS},
            'cdf': [float(np.mean(T <= k)) for k in range(KMAX + 1)]}


NCH = int(float(sys.argv[4])) if len(sys.argv) > 4 else 100000
res = {}
res['theory'] = run_chain(sample_theory, NCH, 1)
res['empirical'] = run_chain(make_emp(g, back_cls), NCH, 2)
half = len(g) // 2
res['emp_first_half'] = run_chain(make_emp(g[:half], back_cls[:half]), NCH, 3)
res['emp_second_half'] = run_chain(make_emp(g[half:], back_cls[half:]), NCH, 4)
# chain started from the actual target classes of the study (seed 72 choice)
rng = np.random.default_rng(72); tg = rng.choice(np.arange(61, len(P) - 1), 20000, replace=False); tg = tg[tg > 150]
st = np.tile(P[tg] % Q, NCH // len(tg) + 1)[:NCH]
res['emp_target_starts'] = run_chain(make_emp(g, back_cls), NCH, 5, start=st)
if len(sys.argv) > 3 and sys.argv[3] != '-':
    P12 = np.load(sys.argv[3]); g12 = np.diff(P12); c12 = IDX[P12[1:] % Q]
    res['emp_1e12_table'] = run_chain(make_emp(g12, c12), NCH, 6)
    out['P_gap0_1e12'] = {q: float(np.mean(g12 % q == 0)) for q in QS}
out['chains'] = {k: {a: b for a, b in v.items() if a != 'cdf'} for k, v in res.items()}

# ---------- second-order memory beyond the class: P(D_2 = 0 mod q) in data vs the empirical first-order chain
mem = {}
for q in QS:
    gm = (np.arange(2, 2 * NGE + 1, 2) % q)
    F = np.stack([Eb[:, gm == r].sum(1) for r in range(q)], 1); F = F / F.sum(1, keepdims=True)   # P(back gap = r mod q | class)
    g1 = g[1:]; cprev = back_cls[:-1]          # p_n = P[i+2]: g1 = gap before p_n, class of p_{n-1} = P[i+1]
    pred = float(np.mean(F[cprev, (-g1) % q]))
    data = float(np.mean((g[1:] + g[:-1]) % q == 0))
    mem[q] = {'data': data, 'first_order_chain': pred, 'diff': data - pred, 'se': math.sqrt(data * (1 - data) / len(g1))}
out['D2_repeat_memory'] = mem

# ---------- round 72: gain vs I_model, with its standard error
mix = TH.mean(0); H = lambda p: float(-np.sum(p[p > 0] * np.log2(p[p > 0])))
out['I_model'] = H(mix) - float(np.mean([H(r) for r in TH]))
gn = P[tg + 1] - P[tg]; ok = gn <= 2 * NG; t2 = tg[ok]; gn = gn[ok]
lr = np.log2(TH[IDX[P[t2] % Q], gn // 2 - 1]) - np.log2(mix[gn // 2 - 1])
out['r72_gain'] = float(lr.mean()); out['r72_gain_se'] = float(lr.std() / math.sqrt(len(lr)))
out['r72_L0'] = float(np.mean(-np.log2(mix[gn // 2 - 1])))
# expected log-loss under the model vs data (calibration of the absolute level)
out['r72_Ltrue_data'] = float(np.mean(-np.log2(TH[IDX[P[t2] % Q], gn // 2 - 1])))
out['r72_Ltrue_model_expect'] = float(np.mean([H(r) for r in TH]))
out['r72_L0_model_expect'] = H(mix)
print(json.dumps(out, indent=1))
json.dump({'out': out, 'cdfs': {k: v['cdf'] for k, v in res.items()}}, open('check74_75.json', 'w'))
