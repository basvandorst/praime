import cypari2
pari = cypari2.Pari()
pari.allocatemem(400*10**6)
pari('bnf=bnfinit(y^2+1,1); gc=gcharinit(bnf,1); gc2=gcharinit(bnf,idealhnf(bnf,2+y)); gc5=gcharinit(bnf,5)')
for g in ['gc','gc2','gc5']:
    print(g, pari(g+'.cyc'))
    n = int(pari('#'+g+'.cyc'))
    for j in range(n):
        v = '[' + ','.join('1' if i==j else '0' for i in range(n)) + ']~'
        print('  basis', j, 'local at inf:', pari('gcharlocal(%s,%s,1)'%(g,v)))
