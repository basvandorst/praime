"""Partial sums over k <= J of the four-point series (rho - 1), of its linear part (2(s0-1) + (s+-1) + (s--1)) and of
the remainder G = (rho - 1) - linear, for d = Q k; and of the three-point part.  Shows which part carries ln J."""
import sys; sys.path.insert(0, '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/twinclass')
import math, json
import numpy as np
import fourpoint as F
out = {}
for Q, Jmax in ((210, 4_000_000), (2310, 4_000_000), (30030, 1_000_000)):
    C = F.constants(Q, 2)
    Js = np.unique(np.logspace(1, math.log10(Jmax), 40).astype(int))
    cum4 = cumL = cumS = cumT = 0.0; res = []
    k0 = 1; chunk = 1 << 21
    tot4 = []; totL = []; totS = []; totT = []
    for a in range(1, Jmax + 1, chunk):
        b = min(Jmax + 1, a + chunk)
        rho, s0, sp, sm, t1, t2 = F.series(Q, 2, a, b, C)
        tot4.append(np.cumsum(rho - 1) + cum4); cum4 = tot4[-1][-1]
        lin = 2 * (s0 - 1) + (sp - 1) + (sm - 1)
        totL.append(np.cumsum(lin) + cumL); cumL = totL[-1][-1]
        totS.append(np.cumsum(2 * (s0 - 1)) + cumS); cumS = totS[-1][-1]
        totT.append(np.cumsum((t1 - 1) + (t2 - 1)) + cumT); cumT = totT[-1][-1]
    c4 = np.concatenate(tot4); cL = np.concatenate(totL); cS = np.concatenate(totS); cT = np.concatenate(totT)
    rows = [(int(J), float(c4[J - 1]), float(cL[J - 1]), float(cS[J - 1]), float(cT[J - 1]), float(c4[J - 1] - cL[J - 1])) for J in Js]
    out[Q] = rows
    print('Q =', Q, ' phi/Q =', F.phi(Q) / Q)
    print('   J        F4(J)      Flin(J)   Fsame(J)   F3(J)     G=F4-Flin   Fsame/lnJ')
    for r in rows[::3] + [rows[-1]]:
        print('%9d %10.4f %10.4f %10.4f %10.4f %10.4f %8.4f' % (r + (r[3] / math.log(r[0]),)))
    # fit F = -a ln J - b on the upper half (J >= 1e3)
    A = np.array([[math.log(r[0]), 1.0] for r in rows if r[0] >= 1000])
    for name, col in (('F4', 1), ('Flin', 2), ('Fsame', 3), ('F3', 4), ('G', 5)):
        y = np.array([r[col] for r in rows if r[0] >= 1000])
        co = np.linalg.lstsq(A, y, rcond=None)[0]
        print('   fit %-5s = %.4f ln J + %.4f' % (name, co[0], co[1]))
json.dump(out, open('partial_sums.json', 'w'))
