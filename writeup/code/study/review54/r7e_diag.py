import math, numpy as np, cypari2, mpmath, time
pari = cypari2.Pari(); pari.allocatemem(1000 * 10 ** 6)
N = 13; PD = 110; NC = 3000
pari.set_real_precision(PD)
pari('mf = mfinit([13,2,Mod(4,13)],0); F = mfeigenbasis(mf)[1]; v = mfcoefs(F, %d); vt = variable(lift(v[3]));' % NC)
pari('z = exp(I*Pi/3); W = vector(%d, n, subst(lift(v[n+1]), vt, z)); Wc = conj(W);' % NC)
print('a(6)-a(2)a(3) =', pari('W[6]-W[2]*W[3]'), ' a(2)=', pari('W[2]'))
# root number from theta relation in GP
pari('th(y, A) = sum(n=1, #A, A[n]*exp(-2*Pi*n*y))'); pari('y0 = 0.9/sqrt(13)'); pari('EPS = th(1/(13*y0), W)/(13*y0^2*th(y0, Wc))'); pari('y1 = 1.4/sqrt(13)'); print('EPS check (y=1.4/sqrtN) diff', pari('abs(EPS - th(1/(13*y1), W)/(13*y1^2*th(y1, Wc)))'))
print('EPS =', pari('EPS'))
pari('Ld = lfuncreate([W, Wc, [0,1], 2, 13, EPS]);')
print('lfuncheckfeq =', pari('lfuncheckfeq(Ld)'))
for t in (7.3, 45.0, 51.7, 60.0, 100.0):
    t0 = time.time()
    print('t=%.1f L(1+it) = %s (%.1fs)' % (t, pari('lfun(Ld, 1 + %s*I)' % t), time.time() - t0))
