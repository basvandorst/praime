import sys; sys.path.insert(0, '/tmp/claude-1000/referee_martin/scripts')
from ref_core import *
# (a) eigenspace sums from PARI eigenforms + AL eigenvalues at small levels 3q
fails = 0; n = 0
for q in [5,7,11,13,17,19,23,29,31,37,41]:
    N = 3*q
    mf = pari.mfinit([N, 2], 0)
    B = pari.mfeigenbasis(mf)
    if len(B) == 0: continue
    AL3 = pari.mfatkineigenvalues(mf, 3); ALq = pari.mfatkineigenvalues(mf, q)
    for l in [p for p in [2,5,7,11,13,17,19,23,29,31,37,41,43,47] if N % p]:
        tq = {Q: tQ(N, Q, l) for Q in [1, 3, q, N]}
        for e3 in (1, -1):
            for eq in (1, -1):
                direct = Fraction(0)
                for i, F in enumerate(B):
                    a3 = AL3[i]; aq = ALq[i]
                    # each orbit: AL eigenvalue list (constant over orbit)
                    if int(a3[0]) != e3 or int(aq[0]) != eq: continue
                    c = pari.mfcoef(F, l)
                    direct += Fraction(str(pari.trace(c))) if str(pari.type(c)) == "t_POLMOD" else Fraction(str(c))
                pred = (tq[1] + e3*tq[3] + eq*tq[q] + e3*eq*tq[N]) / 4
                n += 1
                if pred != direct: fails += 1; print('FAIL', N, l, e3, eq, pred, direct)
print('Prop5.5 vs PARI eigenforms+AL:', n, 'cases, fails', fails)
# (b) big level 3*3989 against trace-form coefficients computed earlier with mftraceform
q = 3989; N = 3*q
ref = {(q,1201): 658, (N,1201): -134, (q,9001): -278, (N,9001): -6334, (1,1201): -90, (3,1201): 206, (1,9001): 58}
for (Q, l), a in ref.items():
    w = (-1)**len(primefactors(Q)) if Q > 1 else 1
    print('N=%d Q=%d l=%d  mine %s  traceform %s' % (N, Q, l, tQ(N, Q, l), w*a))
