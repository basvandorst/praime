# Closed form of the diagonal via the Weil explicit formula (prime sum empty for delta<1):
# D(d) = ln(1/d) + B + 2 ln(1+d) + ((2+d)/d) ln(1+d/2) - 1
import mpmath as mp, json
mp.mp.dps=30
B=1-mp.euler-mp.log(2*mp.pi)
num=json.load(open('dnp.json')); rep=json.load(open('/home/bas/praime/writeup/data/study/r95.json'))['rows']
for row,rr in zip(num,rep):
    d=mp.mpf(row['delta'])
    Dc=mp.log(1/d)+B+2*mp.log1p(d)+(2+d)/d*mp.log1p(d/2)-1
    HL=mp.log(1/d)+B+(1+d)*mp.log1p(d)/d-1
    corr=Dc-HL; ser=mp.mpf(7)/4*d-mp.mpf(7)/8*d**2+mp.mpf(19)/32*d**3
    ref=mp.mpf(7)/4*mp.sqrt(1+d)*mp.log1p(d)**2/d; f0=7/d*(mp.sqrt(1+d)-1)**2
    print('d=%g closed D=%.10f  mine(10k zeros, S(T+) tail) %.10f diff %.2e | r95.json %.10f diff %.2e | (D-HL)/d=%.6f | corr-ser %.1e | corr-refined %.2e (d^3/12=%.2e) corr-(7/8)f0 %.2e'%(
      float(d),Dc,row['D_exact'],row['D_exact']-float(Dc),rr['D'],rr['D']-float(Dc),corr/d,corr-ser,corr-ref,d**3/12,corr-f0))
