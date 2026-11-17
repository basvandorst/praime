import json, numpy as np
from sympy import primerange, isprime
T=100000
PW=list(primerange(2,1584)); assert len(PW)==250
MXS={2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937,21701,23209,44497,86243}
cp=np.zeros(T+1,int); ce=np.zeros(T+1,int); cw=np.zeros(T+1,int)
for p in PW:
    cp[p::p]+=1
    if p in MXS: ce[p::p]+=1
    if isprime(p-2) or isprime(p+2): cw[p::p]+=1
cp[0]=ce[0]=cw[0]=0
r=json.load(open('best.json'))
fit=[x for x in r['rows'] if x[0].startswith('Wheels that fit')][0]
out=('Mersenne finder (11)',max(int(cp.max()),fit[1]),int(np.argmax(cp)) if cp.max()>=fit[1] else fit[2],max(int(ce.max()),fit[3]),int(np.argmax(ce)) if ce.max()>=fit[3] else fit[4],max(int(cw.max()),fit[5]),int(np.argmax(cw)) if cw.max()>=fit[5] else fit[6])
print(out,'| fit alone:',fit)
json.dump({'finder':out},open('best2.json','w'))
