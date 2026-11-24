"""Intermediate quantities on the validation stretch [1.5e11, 1.5e11+1e8) only (no within part is assembled here).
q = 3: b_j = E[chi(p_n) chi(p_{n+j})] (chi = Legendre symbol mod 3), lag correlations rho_j of the stay indicators
s_n = [p_{n+1} = p_n mod 3] (k = 1) and t_n = [p_{n+2} = p_n mod 3] (k = 2), the connected 4-point parts.
q = 3, 5, 7: the triple law D3(a,b,c) and its Markov deviation; the lag covariances tau_j = E[y_n' Pi_W y_{n+j}] of the
within-projected transition increments (k = 1, 2).  Errors: batch means over batches of 20000 primes."""
import numpy as np, math, json
P = np.load('data/val_primes.npy'); X0 = 150_000_000_000; SPAN = 100_000_000
NW = int(np.searchsorted(P, X0 + SPAN)); LAM = math.log(X0 + SPAN / 2)
JMAX = 20; BATCH = 20000
def bm(v):
    """mean and batch-means error of a per-n series"""
    nb = len(v) // BATCH; m = v[:nb * BATCH].reshape(nb, BATCH).mean(1)
    return float(v.mean()), float(m.std(ddof=1) / math.sqrt(nb))
out = {'lam': LAM, 'nprimes': NW}
# ---------- q = 3 ----------
chi = np.where(P % 3 == 1, 1.0, -1.0)
b = {}; 
for j in range(1, JMAX + 3):
    b[j] = bm(chi[:NW] * chi[j:NW + j])
out['b'] = {j: b[j] for j in b}
s = (chi[:NW + JMAX] * chi[1:NW + JMAX + 1] > 0).astype(float)      # stays, k = 1
t = (chi[:NW + JMAX] * chi[2:NW + JMAX + 2] > 0).astype(float)      # k = 2
def lagcorr(x):
    m = x[:NW].mean(); v = x[:NW].var(); r = {}
    for j in range(1, JMAX + 1):
        r[j] = bm((x[:NW] - m) * (x[j:NW + j] - m) / v)
    return r
out['rho_k1'] = lagcorr(s); out['rho_k2'] = lagcorr(t)
# connected 4-point part of chi_n chi_{n+1} chi_{n+j} chi_{n+j+1} beyond the Wick (pairing) approximation, j >= 2
c4 = {}
for j in range(2, 11):
    e4 = bm(chi[:NW] * chi[1:NW + 1] * chi[j:NW + j] * chi[j + 1:NW + j + 1])
    wick = b[1][0] ** 2 + b[j][0] ** 2 + b[j + 1][0] * b[j - 1][0]
    c4[j] = {'E4': e4, 'wick': wick, 'conn': e4[0] - wick}
out['c4_k1'] = c4
# ---------- general q: triple law and tau_j ----------
def classes(q):
    r = [a for a in range(q) if math.gcd(a, q) == 1]; lut = np.full(q, -1); lut[r] = np.arange(len(r))
    return r, lut[P % q]
def proj_within(Pm, pi):
    C = len(pi); R = np.zeros((C * C, C * C)); V = np.zeros((C * C, C))
    for a in range(C):
        sp = np.sqrt(Pm[a]); R[a * C:(a + 1) * C, a * C:(a + 1) * C] = np.eye(C) - np.outer(sp, sp)
        for bb in range(C): V[a * C + bb, bb] = math.sqrt(pi[a] * Pm[a, bb])
    M = V.T @ R @ V
    return R - R @ V @ np.linalg.pinv(M) @ V.T @ R
res = {}
for q in (3, 5, 7):
    r, c = classes(q); C = len(r)
    for k in (1, 2):
        a0 = c[:NW + JMAX + 2]; a1 = c[k:NW + JMAX + 2 + k]
        X = np.zeros((C, C)); np.add.at(X, (a0[:NW], a1[:NW]), 1)
        pi = X.sum(1) / X.sum(); Pm = X / X.sum(1, keepdims=True)
        PiW = proj_within(Pm, pi)
        # y(a,b) in R^{C*C}: (1[c_n=a, c_{n+k}=b] - 1[c_n=a] P_ab)/sqrt(pi_a P_ab); only C*C distinct vectors
        Yv = np.zeros((C * C, C * C))
        for a in range(C):
            for bb in range(C):
                for b2 in range(C):
                    Yv[a * C + bb, a * C + b2] -= Pm[a, b2] / math.sqrt(pi[a] * Pm[a, b2])
                Yv[a * C + bb, a * C + bb] += 1 / math.sqrt(pi[a] * Pm[a, bb])
        G = (Yv @ PiW) @ Yv.T          # G[i, i'] = y(i)' Pi_W y(i')
        idx = a0 * C + a1
        tau = {0: bm(G[idx[:NW], idx[:NW]])}
        for j in range(1, 11):
            tau[j] = bm(G[idx[:NW], idx[j:NW + j]])
        d = {'C': C, 'P': Pm.tolist(), 'pi': pi.tolist(), 'tau': tau}
        if k == 1:
            a2 = c[2:NW + 2]
            D3 = np.zeros((C, C, C)); np.add.at(D3, (a0[:NW], a1[:NW], a2), 1); D3 /= NW
            D2 = D3.sum(2)
            Delta = D3 - D2[:, :, None] * Pm[None, :, :]       # D3(a,b,c) - D2(a,b) P(c|b)
            d['Delta'] = Delta.tolist(); d['D3'] = D3.tolist()
            # Delta normalised the way the first-order law predicts: phi^3 Delta(a,b,c) vs (1/(2 lam))(1 - phi delta_ac)
            d['Delta_scaled_same'] = float(np.mean([Delta[a, bb, a] for a in range(C) for bb in range(C)]) * C ** 3)
            d['Delta_scaled_diff'] = float(np.mean([Delta[a, bb, cc] for a in range(C) for bb in range(C) for cc in range(C) if cc != a]) * C ** 3)
        res['%d_%d' % (q, k)] = d
out['tables'] = res
json.dump(out, open('data/val_measure.json', 'w'), indent=1)
print('lam = %.4f, primes %d' % (LAM, NW))
print('b_j (q=3):', ' '.join('%d:%+.5f(%.5f)' % (j, b[j][0], b[j][1]) for j in range(1, 13)))
print('LOS: b1 = -ln(2 pi lam/3)/(2 lam) = %+.5f;  b_j = -1/(2(j-1)lam): ' % (-math.log(2 * math.pi * LAM / 3) / (2 * LAM)) + ' '.join('%+.5f' % (-1 / (2 * (j - 1) * LAM)) for j in range(2, 8)))
print('rho_j stays k=1:', ' '.join('%d:%+.5f(%.5f)' % (j, *out['rho_k1'][j]) for j in range(1, 11)))
print('rho_j k=2      :', ' '.join('%d:%+.5f(%.5f)' % (j, *out['rho_k2'][j]) for j in range(1, 11)))
print('identity check rho_1 vs (b2 - b1^2)/(1 - b1^2) = %+.5f' % ((b[2][0] - b[1][0] ** 2) / (1 - b[1][0] ** 2)))
print('connected 4-pt (k=1):', ' '.join('%d:%+.5f(%.5f)' % (j, c4[j]['E4'][0] - c4[j]['wick'], c4[j]['E4'][1]) for j in range(2, 8)))
for key, d in res.items():
    C = d['C']; tau = d['tau']
    print('q,k=%s  tau0/(C-1)^2 = %.4f   tau_j/(C-1)^2:' % (key, tau[0][0] / (C - 1) ** 2) + ' '.join(' %d:%+.5f(%.5f)' % (j, tau[j][0] / (C - 1) ** 2, tau[j][1] / (C - 1) ** 2) for j in range(1, 11)))
    if 'Delta_scaled_same' in d:
        print('   phi^3*Delta mean same(a=c) %+.5f  diff %+.5f   first-order law: %+.5f / %+.5f' % (d['Delta_scaled_same'], d['Delta_scaled_diff'], (1 - C) / (2 * LAM), 1 / (2 * LAM)))
