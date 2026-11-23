import cypari2, time, sys, math
pari=cypari2.Pari(); pari.allocatemem(800*10**6); pari.default("realprecision",28)
pols=open('/home/bas/praime/writeup/code/study/lowzero/f5list.txt').read().split('\n')
idx=[int(a) for a in sys.argv[1:]]
for i in idx:
    f=pols[i]; t=time.time()
    P=pari(f)
    print(f, 'galois', pari.polgalois(P), 'sig', pari.nfinit(P)[1], 'disc', pari.nfdisc(P))
    L=pari.lfundiv(pari.lfuncreate(P),pari.lfuncreate(1))
    pr=pari.lfunparams(L); print(' params',pr)
    z=pari.lfunzeros(L,4)
    print(' zeros',[float(x) for x in z][:3])
    # L'/L(1)
    try:
        v=pari.lfun(L,1); d=pari.lfun(L,1,1); print(' L1',float(v),'Lp/L',float(d)/float(v))
    except Exception as e: print('err',e)
    # check via symmetric difference of log zeta_K - log zeta
    h=1e-6
    print(' time',time.time()-t, flush=True)
