# singular series S(h) for 0<=h<HMAX, own code
import numpy as np, sys
HMAX = int(float(sys.argv[1])) if len(sys.argv)>1 else 10**7
C2 = 0.66016181584686957392781211
def primes_upto(n):
    s = np.ones(n+1, bool); s[:2]=False
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0]
P = primes_upto(HMAX)
S = np.zeros(HMAX, np.float64)
S[0::2] = 2*C2
for p in P[1:]:
    S[p::p] *= (p-1)/(p-2)
S[0] = 0.0  # unused
np.save('S_%d.npy'%HMAX, S)
print('done', HMAX, S[2], S[6], S[30])
