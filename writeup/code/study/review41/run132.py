import json, random
from sympy import isprime
D='/home/bas/praime/writeup/code/study/lowzero/'
rows=json.load(open(D+'r132_part0.json'))+json.load(open(D+'r132_part1.json'))
random.seed(7); smp=random.sample(rows,12)+sorted(rows,key=lambda r:-r['run'])[:3]
bad=0
for r in smp:
    d=r['d']; best=0
    for A in range(1,301):
        for B in range(2*A):
            if (B*B-d)%(4*A): continue
            C=(B*B-d)//(4*A); cur=0
            for y in range(-80,81):
                v=abs(A*y*y+B*y+C)
                if v>1 and isprime(v): cur+=1; best=max(best,cur)
                else: cur=0
    bad+= best!=r['run']; print(d,r['run'],best)
print('mismatch',bad)
