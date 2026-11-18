"""A guess tested on every ring mistake below 1e10: does every prime factor p of the number n have a power p^j = n (mod r)?
(It would follow if a mistake could only arise from the Frobenius map.) Reads ring_psp2_to1e10.json."""
import json
from sympy import factorint, n_order
res=json.load(open('ring_psp2_to1e10.json'))
tot=0; viol=0
for r in range(5,41):
    L=res[str(r)]; v=[]
    for n in L:
        m=n%r; fs=list(factorint(n))
        if not all(any(pow(p,j,r)==m for j in range(1,n_order(p,r)+1)) for p in fs): v.append((n,m,[(p,p%r) for p in fs]))
    tot+=len(L); viol+=len(v)
    if v: print("ring of %d cells: %d mistakes, the guess fails for %d; first: %s"%(r,len(L),len(v),v[:2]))
    if r==5: print("ring of 5 cells: guess fails for %d of %d"%(len(v),len(L)))
print("rings of 5 to 40 cells: %d mistakes (number, ring), the guess fails for %d"%(tot,viol))
