import sys
from sympy import primerange
from collections import Counter
k=int(sys.argv[1])
P=list(primerange(10**(k-1),10**k)); S=set(P); idx={p:i for i,p in enumerate(P)}
adj=[[] for _ in P]
for i,p in enumerate(P):
    s=str(p)
    for pos in range(k):
        for d in '0123456789':
            if d==s[pos] or (pos==0 and d=='0'): continue
            q=int(s[:pos]+d+s[pos+1:])
            if q in S and q>p: j=idx[q]; adj[i].append(j); adj[j].append(i)
deg=Counter(len(a) for a in adj)
print('k',k,'nodes',len(P),'edges',sum(len(a) for a in adj)//2,'degree histogram',sorted(deg.items())[:8])
leaves=[i for i,a in enumerate(adj) if len(a)==1]
from collections import defaultdict
par=defaultdict(list)
for l in leaves: par[adj[l][0]].append(l)
print('leaves',len(leaves),'leaf-parents with >=2 leaves',sum(1 for v in par if len(par[v])>=2))
with open(f'g{k}.txt','w') as f:
    f.write(f"{len(P)}\n")
    for i,a in enumerate(adj): f.write(str(P[i])+' '+' '.join(map(str,a))+'\n')
