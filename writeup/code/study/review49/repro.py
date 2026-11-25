import cypari2, numpy as np, json
pari=cypari2.Pari(); pari.allocatemem(300*10**6)
Z=np.load(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/data/study/r144/fresh_zeros.npz')
for ainv,cid,(a,b) in (([0,-1,1,-3,1],'N2357_0_m1_1_m3_1',(78,90)),([1,0,1,-8,7],'N427_1_0_1_m8_7',(210,224))):
    L=pari.lfuncreate(pari.ellinit(ainv))
    for divz in (8,64):
        z=np.array([float(x) for x in pari.lfunzeros(L,[a,b],divz)])
        print(cid, 'divz',divz, 'zeros in',(a,b), len(z), np.round(z[(z>a+4)&(z<b-4)],5))
    g=Z[cid]; print(' stored', np.round(g[(g>a+4)&(g<b-4)],5))
    # Z-function sign check around the close pair
    pair=np.sort(g[(g>a+4)&(g<b-4)])
    d=np.diff(pair); j=np.argmin(d[d>0]) if (d>0).any() else None
