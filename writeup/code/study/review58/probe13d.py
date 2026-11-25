import time, cypari2
pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6)
L = pari('lfungenus2([(-1)*x^2+(-1)*x^3+(1)*x^4+(1)*x^5, x^0])')
t0 = time.time(); rr = pari.lfunrootres(L); eps = rr[2]; print('eps', eps, '%.0f s' % (time.time() - t0), flush=True)
pari('Lg = lfungenus2([(-1)*x^2+(-1)*x^3+(1)*x^4+(1)*x^5, x^0])'); print('L[6] before', pari('Lg[6]'))
pari('Lg[6] = %s' % eps); print('L[6] after', pari('Lg[6]'), flush=True)
t0 = time.time()
try:
    Li = pari('Li = lfuninit(Lg, [1, 153])'); print('lfuninit [1,153] with eps fixed ok %.0f s' % (time.time() - t0), flush=True)
    print('  feq check', pari('lfuncheckfeq(Li)'))
except cypari2.PariError as e:
    print('FAIL %.0f s: %s' % (time.time() - t0, str(e)[:70]))
