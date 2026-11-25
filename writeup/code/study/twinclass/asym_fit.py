"""Leading asymptotics of the partial sums F4(J) = sum_{k<=J} (rho(Qk) - 1) (partial_sums.json):
derived  F4(J) ~ -(H_Q/4) ln^2 J + ...,  H_Q = (phi(Q)/Q)^2 / prod_{p not | Q} (1 - 1/(p-1)^2)   (twins, h = 2);
linear part  Flin(J) ~ -(phi(Q)/Q) ln J + const.  Fit F4 = -a ln^2 J + b ln J + c on J >= 300 and compare a with H_Q/4."""
import sys; sys.path.insert(0, '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/twinclass')
import json, math
import numpy as np
import fourpoint as F
D = json.load(open('partial_sums.json'))
out = {}
for Q, rows in D.items():
    Q = int(Q)
    C = F.constants(Q, 2)
    HQ = (F.phi(Q) / Q) ** 2 / math.exp(C['logC'])
    J = np.array([r[0] for r in rows], float); F4 = np.array([r[1] for r in rows]); FL = np.array([r[2] for r in rows])
    m = J >= 300
    L = np.log(J[m])
    A = np.column_stack([-L ** 2, L, np.ones_like(L)])
    co, *_ = np.linalg.lstsq(A, F4[m], rcond=None)
    # with the ln^2 coefficient fixed at the derived value
    A2 = np.column_stack([L, np.ones_like(L)])
    co2, *_ = np.linalg.lstsq(A2, F4[m] + HQ / 4 * L ** 2, rcond=None)
    res2 = F4[m] + HQ / 4 * L ** 2 - A2 @ co2
    coL, *_ = np.linalg.lstsq(np.column_stack([-L, np.ones_like(L)]), FL[m], rcond=None)
    out[Q] = {'H_Q/4 derived': HQ / 4, 'a fitted': co[0], 'b': co[1], 'c': co[2], 'fixed-a fit b,c': list(co2), 'rms resid fixed a': float(res2.std()),
              'lin slope fitted': coL[0], 'phi/Q': F.phi(Q) / Q}
    print(Q, {k: (round(v, 4) if isinstance(v, float) else [round(x, 4) for x in v]) for k, v in out[Q].items()})
json.dump(out, open('asym_fit.json', 'w'), indent=1)
