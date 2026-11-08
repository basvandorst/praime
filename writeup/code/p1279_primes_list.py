import sys
from sympy import factorint, n_order, divisors
L=int(sys.argv[1]); B=2**32
out={}
for d in divisors(L):
    if 2<=d<=64:
        for p in factorint(2**d-1):
            if p<B and n_order(2,p)==d: out[p]=d
for line in open(sys.argv[2]):
    p,d=map(int,line.split()); out[p]=d
for p in sorted(out): print(p,out[p])
s=sum(d/p for p,d in out.items())
print('#',len(out),'primes; sum m/p =',s, file=sys.stderr)
