import cypari2
pari = cypari2.Pari(); pari.allocatemem(400*10**6); pari.set_real_precision(30)
for D,a,b in [(41,583.5,584.2),(-47,482.3,482.8)]:
    L = pari('lfuncreate(%d)'%D)
    print(D, pari.lfunzeros(L,[a,b]), [pari('lfunhardy(lfuncreate(%d), %s)'%(D,t)) for t in ['583.84','482.535']])
