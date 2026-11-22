import gmpy2, json
L0,L1=10**9,int(1.01e9)
tw=[p for p in range(L0+1,L1,2) if gmpy2.is_prime(p) and gmpy2.is_prime(p+2)]
r=json.load(open('t_small.json'))
from collections import Counter
c=Counter(b-a for a,b in zip(tw,tw[1:]))
print(len(tw), r['twins'], {str(k):v for k,v in sorted(c.items())}==r['data'])
