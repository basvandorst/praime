import cypari2, math, time
pari=cypari2.Pari(); pari.allocatemem(800*10**6)
def g1(d):
    Lf=pari.lfuncreate(d); T=3*2*math.pi/math.log(abs(d)); z=[]
    while not z:
        z=[float(t) for t in pari.lfunzeros(Lf,T) if float(t)>1e-9]; T*=2
    return min(z)*math.log(abs(d))/(2*math.pi)
for d in (-163,7213,1398053,979373,-991027,5417453,73678233):
    t=time.time(); print(d, round(g1(d),4), round(time.time()-t,2), flush=True)
