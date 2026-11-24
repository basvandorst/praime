"""Lag-1 term tau_1 of the k = 1 within part from the consecutive-prime triple law.
Exact algebra (proved, stationary sequence):  E[U_n(a,b) U_{n+1}(a',b')] = (delta_{a'b} - P_ab) Delta(a,a',b'),
Delta(a,a',b') = D3(a,a',b') - D2(a,a') P(b'|a'),  U_n(a,b) = 1[c_n=a, c_{n+1}=b] - 1[c_n=a] P_ab,
tau_1 = sum Pi_W[(ab),(a'b')] E[y_n(a,b) y_{n+1}(a',b')],  y = U/sqrt(pi_a P_ab).
HL input: D3 per fine class u mod Q in the order-1 conditioned expansion with the pair factors S_{>y} of every rough
site to every earlier point of the pattern (x, and the first neighbour x + t1)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, math
from hlcore import primorial, rough_classes, S_tail_table
from chainlaw import proj_within

def tau1_from_D3(D3):
    C = D3.shape[0]
    D3 = D3 / D3.sum()
    D2f = D3.sum(2); D2s = D3.sum(0)
    pi = D2f.sum(1); P1 = D2f / pi[:, None]
    P2 = D2s / D2s.sum(1, keepdims=True)
    Delta = D3 - D2f[:, :, None] * P2[None, :, :]
    P = P1
    E = np.zeros((C, C, C, C))       # E[U(a,b) U'(a',b')]
    for a in range(C):
        for b in range(C):
            for a2 in range(C):
                E[a, b, a2, :] = ((1.0 if a2 == b else 0.0) - P[a, b]) * Delta[a, a2, :]
    norm = np.sqrt(pi[:, None] * P)                      # (a,b)
    Ey = E / norm[:, :, None, None] / norm[None, None, :, :]
    PiW = proj_within(P, pi)
    tau1 = float(np.einsum('ij,ij->', PiW, Ey.reshape(C * C, C * C)))
    return tau1 / (C - 1) ** 2, Delta, P, pi

def D3_model(y, lam, qs, T1=420, T2=650, chunk=12):
    Q = primorial(y); U = rough_classes(Q); n = len(U); rho0 = (Q / n) / lam
    St = S_tail_table(T1 + T2 + 2, y)
    isr = np.zeros(Q, bool); isr[U] = True
    t1 = np.arange(1, T1 + 1); s = np.arange(1, T2 + 1)
    out = {q: np.zeros((q, q, q)) for q in qs}
    for c0 in range(0, n, chunk):
        uu = U[c0:c0 + chunk]; m = len(uu)
        pos1 = (uu[:, None] + t1[None, :]) % Q
        p1 = np.where(isr[pos1], rho0 * St[t1][None, :], 0.0)
        surv1 = np.cumprod(1 - p1, axis=1); surv1 = np.concatenate([np.ones((m, 1)), surv1[:, :-1]], 1)
        w1 = p1 * surv1                                                     # (m, T1)
        pos2 = (pos1[:, :, None] + s[None, None, :]) % Q                    # (m, T1, T2)
        p2 = np.where(isr[pos2], rho0 * St[s][None, None, :] * St[t1[:, None] + s[None, :]][None, :, :], 0.0)
        surv2 = np.cumprod(1 - p2, axis=2); surv2 = np.concatenate([np.ones((m, T1, 1)), surv2[:, :, :-1]], 2)
        w = w1[:, :, None] * p2 * surv2                                     # joint P(t1, s)
        for q in qs:
            a = (uu % q)[:, None, None]; b = (uu[:, None] + t1[None, :]) % q; c = (b[:, :, None] + s[None, None, :]) % q
            idx = (np.broadcast_to(a, w.shape) * q + np.broadcast_to(b[:, :, None], w.shape)) * q + c
            out[q] += np.bincount(idx.ravel(), weights=w.ravel(), minlength=q ** 3).reshape(q, q, q)
    res = {}
    for q in qs:
        red = [r for r in range(q) if math.gcd(r, q) == 1]
        res[q] = out[q][np.ix_(red, red, red)]
    return res, n

if __name__ == '__main__':
    import json, time
    lam = math.log(1.5e11 + 5e7)
    meas = json.load(open('data/val_measure.json'))
    for q in (3, 5, 7):     # algebra check on the measured triple law
        D3m = np.array(meas['tables']['%d_1' % q]['D3'])
        t1m = tau1_from_D3(D3m)[0]
        print('q=%d measured D3 -> tau1/(C-1)^2 = %+.5f   (direct lag covariance %+.5f)' % (q, t1m, meas['tables']['%d_1' % q]['tau']['1'][0] / (len(D3m) - 1) ** 2))
    for y in (11, 13):
        t0 = time.time(); D3s, n = D3_model(y, lam, (3, 5, 7))
        for q in (3, 5, 7):
            D3 = D3s[q]; C = D3.shape[0]
            t1, Delta, P, pi = tau1_from_D3(D3)
            same = np.mean([Delta[a, b, a] for a in range(C) for b in range(C)]) * C ** 3 / D3.sum()
            print('y=%d q=%d  model tau1/(C-1)^2 = %+.5f   phi^3 Delta(same) %+.5f   mass %.6f' % (y, q, t1, same, D3.sum() / n))
        print('  (%.1fs)' % (time.time() - t0), flush=True)
