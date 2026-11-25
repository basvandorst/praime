"""The spread law in closed form, checked against the exact Hardy-Littlewood sums on every seen window.
  v_X = [sum r2(1-r2) + 2 sum_J F(J) (w(J) - w(J+1))] / sum r2,    F(J) = -a ln^2 J + b ln J + c
  a = H_Q/4 (derived: H_Q = (phi(Q)/Q)^2 / prod_{p not | Q}(1 - 1/(p-1)^2), so c_Q H_Q = 1),  b, c: constants of the
  HL series itself (least squares on its exact partial sums for 300 <= J <= 4e6 with a fixed; asym_fit.json).
Constant-lambda form (window short against x): v_X = 1 - r2 - 2 r2 [a (L^2 - 2L + 2) - b (L - 1) - c],  L = ln(W/Q),
leading term 1 - L^2/(2 lambda^2) (universal).  Linear part (round 135 law for each member's class):
  P_lin = 2 kbar v_Y,  v_Y = 1 - [ln(W/Q) + S_Q - B]/lambda,  S_Q = sum_{p|Q} ln p/(p-1),  B = 1 - gamma - ln 2 pi."""
import sys; sys.path.insert(0, '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/twinclass')
import json, math
import numpy as np
import fourpoint as F
AF = json.load(open('asym_fit.json')); D = json.load(open('all_decomp.json'))
EG = 0.5772156649015329; B = 1 - EG - math.log(2 * math.pi)
WIN = {'1e7-1e8': (10**7, 9 * 10**7), '1e8-1e9': (10**8, 9 * 10**8), '1e9-1e10': (10**9, 9 * 10**9),
       'A': (10**10, 5 * 10**7), 'B': (10**12, 5 * 10**7), 'C': (10**15, 2 * 10**7)}
def v_closed(Q, L0, W):
    C = F.constants(Q, 2); a = AF[str(Q)]['H_Q/4 derived']; b, c = AF[str(Q)]['fixed-a fit b,c']
    N = W // Q; wf = F.profile_w(L0, Q, 2, C['c'], N)
    J = np.unique(np.concatenate([np.arange(1, 2000), np.logspace(np.log10(2000), np.log10(N - 1), 20000).astype(np.int64)]))
    # sum_J F(J)(w(J) - w(J+1)) ~ integral of F dw on a fine log grid (F smooth beyond 2000)
    FJ = -a * np.log(J) ** 2 + b * np.log(J) + c
    wJ = wf(J.astype(float)); wnext = np.append(wJ[1:], 0.0)
    off = 2 * float(np.sum(FJ * (wJ - wnext)))
    m = L0 + Q / 2 + Q * np.arange(0, N, max(1, N // 200000), dtype=float)
    r2 = C['c'] / (np.log(m) * np.log(m + 2)); mX = r2.mean() * N; diag = float((r2 * (1 - r2)).mean() * N)
    return (diag + off) / mX
def v_const(Q, L0, W):
    C = F.constants(Q, 2); a = AF[str(Q)]['H_Q/4 derived']; b, c = AF[str(Q)]['fixed-a fit b,c']
    lam = math.log(L0 + W / 2); r2 = C['c'] / lam ** 2; L = math.log(W / Q)
    return 1 - r2 - 2 * r2 * (a * (L * L - 2 * L + 2) - b * (L - 1) - c), 1 - L * L / (2 * lam * lam)
out = {}
print('window      Q     exact   closed(sum)  const-lambda  leading | P_lin exact  2 kbar v_Y(r135)')
for w, (L0, W) in WIN.items():
    for Q in (210, 2310, 30030):
        e = D[w]['h2_%d' % Q]; ex = e['four_point_only']
        cs = v_closed(Q, L0, W); cl, lead = v_const(Q, L0, W)
        lam = math.log(L0 + W / 2); SQ = sum(math.log(p) / (p - 1) for p in F.qfactors(Q))
        vY = 1 - (math.log(W / Q) + SQ - B) / lam; plin = 2 * e['HL']['kbar'] * vY
        out['%s|%d' % (w, Q)] = {'exact': ex, 'closed_sum': cs, 'const_lambda': cl, 'leading': lead, 'P_lin_exact': e['HL']['vZ'], 'P_lin_closed': plin}
        print('%-9s %6d  %.4f   %.4f      %.4f      %.3f  |  %.4f      %.4f' % (w, Q, ex, cs, cl, lead, e['HL']['vZ'], plin))
json.dump(out, open('closed_form.json', 'w'), indent=1)
# empirical scatter of the linear-part ratio over the 27 seen sub-windows
for Q in (210, 2310, 30030):
    r = np.array([D[w]['h2_%d' % Q]['measured_u2Z'] / D[w]['h2_%d' % Q]['pred_u2Z'] for w in D if '/sub' in w])
    K = D['A']['h2_%d' % Q]['K']
    print('u2Z ratio over 27 sub-windows, Q=%d: mean %.3f sd %.3f (nominal %.3f)' % (Q, r.mean(), r.std(ddof=1), math.sqrt(2 / (K - 1))))
