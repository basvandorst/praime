import sys, time, cypari2
pari = cypari2.Pari(); pari.allocatemem(600 * 10 ** 6)
L = pari('lfungenus2([(-1)*x^2+(-1)*x^3+(1)*x^4+(1)*x^5, x^0])')
for dom in ([0, 153], [0.5, 153]):
    t0 = time.time()
    try:
        Li = pari.lfuninit(L, dom); print('lfuninit', dom, 'ok %.0f s' % (time.time() - t0), flush=True)
        z = pari.lfunzeros(Li, [140, 151]); print('  zeros in [140,151]:', len(z), [round(float(x), 6) for x in z][:4], flush=True); del Li, z
    except cypari2.PariError as e:
        print('lfuninit', dom, 'FAIL %.0f s: %s' % (time.time() - t0, str(e)[:70]), flush=True)
