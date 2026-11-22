import cypari2, math
from scipy.special import digamma
pari = cypari2.Pari(); pari.allocatemem(800*10**6)
for D,T in [(-4,600),(5,600),(-163,300)]:
    L=pari.lfuncreate(D); z=[float(x) for x in pari.lfunzeros(L,T)]
    q=abs(D); a=1 if D<0 else 0
    s=sum(2/(0.25+x*x) for x in z)
    tail=(math.log(q*T/(2*math.pi))+1)/(math.pi*T)
    S=2*float(pari.lfun(L,1,1))/float(pari.lfun(L,1))+math.log(q/math.pi)+digamma((1+a)/2)
    print(D,len(z),'zeros sum',s,'+tail',s+tail,'S',S, 'diff', s+tail-S)
