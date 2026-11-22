"""Round 13, addendum 2: (a) Markov-chain null, (b) chi^2/dof against the window length."""
import os, json, math, pickle
import numpy as np
import theory2 as T2
import theory_r3 as R3
import common as C
import r7, r13

cache = pickle.load(open(r7.CACHE7, 'rb'))
Q = 210; U = r7.units(Q); idx = {a: i for i, a in enumerate(U)}
out = {}

# (a) Markov null at the train window
L, W = R3.WIN['train']
u = np.linspace(math.log(L), math.log(L + W), 401); w = np.exp(u) / u
hs = np.arange(2, r7.HM + 1, 2)
Ph = np.zeros((len(U), len(hs)))
for i, a in enumerate(U):
    S0 = float(T2.S_cond(a, Q)(np.array([[0]]))[0])
    for j, h in enumerate(hs):
        SB, la = cache[(Q, a, int(h))]
        if SB:
            Ph[i, j] = np.trapezoid(w * SB / S0 / u * np.exp(np.interp(u, R3.LAMG, la)), u) / np.trapezoid(w, u)
Ph /= Ph.sum(1, keepdims=True)
cum = np.cumsum(Ph, 1)
cls = np.minimum(hs // 2, 20) - 1
P1 = np.zeros((len(U), 20))
for j in range(len(hs)):
    P1[:, cls[j]] += Ph[:, j]
nxt = np.array([[idx[(a + int(h)) % Q] if math.gcd(a + int(h), Q) == 1 else -1 for h in hs] for a in U])
rng = np.random.default_rng(7)
n_steps = 5_096_875
vals = []
for run in range(10):
    s = 0
    r = rng.random(n_steps)
    N = np.zeros((len(U), 20))
    for k in range(n_steps):
        j = int(np.searchsorted(cum[s], r[k]))
        j = min(j, len(hs) - 1)
        N[s, cls[j]] += 1
        s = nxt[s, j]
    c, d = r13.chi2(N, P1); vals.append(c / d)
    print('markov run', run, round(c / d, 3), flush=True)
out['markov_null'] = {'mean': float(np.mean(vals)), 'sd': float(np.std(vals, ddof=1))}

# (b) chi^2/dof against the window length (train)
M = C.primes_mask(L, W); Pr = np.nonzero(M)[0].astype(np.int64) + L
for Wsub in (9 * 10 ** 6, 9 * 10 ** 5):
    vs = []
    for k in range(10):
        Ls = L + k * 9 * 10 ** 6
        sel = (Pr[:-1] >= Ls) & (Pr[:-1] < Ls + Wsub)
        st = np.array([idx[int(a)] for a in (Pr[:-1][sel] % Q)])
        cl = np.minimum(np.diff(Pr)[sel] // 2, 20) - 1
        N = np.zeros((len(U), 20)); np.add.at(N, (st, cl), 1)
        T = r7.trans(cache, Q, Ls, Wsub); Pa = T.sum(2).T; Pa = Pa / Pa.sum(1, keepdims=True)
        c, d = r13.chi2(N, Pa); vs.append(c / d)
    out['W=%d' % Wsub] = {'mean': float(np.mean(vs)), 'se': float(np.std(vs, ddof=1) / math.sqrt(10)), 'values': vs}
    print('W', Wsub, out['W=%d' % Wsub]['mean'], out['W=%d' % Wsub]['se'], flush=True)
json.dump(out, open(os.path.join(R3.OUT, 'r13c.json'), 'w'), indent=1)
