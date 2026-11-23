import cypari2, math
pari=cypari2.Pari(); pari.allocatemem(500*10**6)
ds={}
for line in open('/home/bas/praime/writeup/code/study/polys/dist300.txt'):
    A,B,C,s,L=map(int,line.split())
    if L>=36:
        D=B*B-4*A*C; d=int(pari.coredisc(D)); ds.setdefault(d,[]).append(L)
for d,r in ds.items():
    chis=[int(pari.kronecker(d,p)) for p in (2,3,5,7,11,13,17,19,23)]
    print(d, r, d%8, len(pari.factor(abs(d))[0]), chis)
print(len(ds))
