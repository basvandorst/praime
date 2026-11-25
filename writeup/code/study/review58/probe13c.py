import cypari2
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6)
print('realprecision', pari.default('realprecision'), 'bitprecision', pari('default(realbitprecision)'))
for N, f, h in ((1403, '(-1)*x^2+(-1)*x^3+(1)*x^4+(1)*x^5', 'x^0'), (1383, '(1)*x^4+(-1)*x^5', 'x^0'),
                (1549, '(-1)*x^0+(-1)*x^1+(1)*x^2+(1)*x^3+(-1)*x^4+(-1)*x^5', 'x^0+x^2+x^3'), (2981, None, None)):
    if f is None: continue
    L = pari('lfungenus2([%s,%s])' % (f, h))
    print(N, 'cost [1,120]', pari.lfuncost(L, [1, 120]), ' [1,153]', pari.lfuncost(L, [1, 153]), ' disc', pari('hyperelldisc([%s,%s])' % (f, h)), 'factor N', pari.factor(N))
    print('   ldata[2..6]', [str(pari('%s' % 0)) for _ in ()], pari('my(L=lfungenus2([%s,%s])); [L[3], L[4], L[5], L[6]]' % (f, h)))
