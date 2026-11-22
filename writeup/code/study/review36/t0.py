import cypari2, time
pari = cypari2.Pari(); pari.allocatemem(400*10**6)
for D in [-3,-4,5,-163,-4996+1, 4993, -4999]:
    if not pari.isfundamental(D): print(D,'notfund'); continue
    t=time.time()
    L = pari.lfuncreate(D)
    z = pari.lfunzeros(L, 30)
    v = pari.lfun(L,1); d = pari.lfun(L,1,1)
    print(D, len(z), [float(x) for x in z[:3]], float(v), float(d), round(time.time()-t,2))
