import sys, itertools
from collections import defaultdict
X=int(float(sys.argv[2])); g=defaultdict(list)
for line in open(sys.argv[1]):
    q,j,p,a=map(int,line.split()); g[q].append((p,a,j))
sols=[]
for q,L in g.items():
    for r in range(1,len(L)+1):
        for S in itertools.combinations(L,r):
            if len({p for p,_,_ in S})<r: continue
            n=1
            for p,a,_ in S: n*=p**a
            k=sum(j for *_,j in S)
            if k>1 and n<=X: sols.append((n,q,k,S))
print(len(g),"q values;", "multi-component q:",sorted(q for q,L in g.items() if len(L)>1)[:20])
for s in sorted(sols): print(s)
