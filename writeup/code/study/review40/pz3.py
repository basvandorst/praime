import cypari2
pari = cypari2.Pari(); pari.allocatemem(400*10**6); pari.set_real_precision(40)
for D, ts in [(-47, [482.50,482.52,482.528,482.535,482.542,482.55,482.57]), (-59,[344.34,344.356,344.37,344.386,344.40]), (41,[583.80,583.82,583.84,583.86,583.88])]:
    L = pari('lfuncreate(%d)'%D)
    print(D, [(t, float(pari.lfunhardy(L, t))) for t in ts])
