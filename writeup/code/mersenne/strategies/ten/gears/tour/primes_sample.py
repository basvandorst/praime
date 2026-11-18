"""Writes primes_sample.txt: 3,000 primes of each bit length from 8 to 64 (next prime after a random start, seed 11)."""
import random, gmpy2
random.seed(11)
out=[]
for b in range(8,65):
    for _ in range(3000):
        lo=1<<(b-1); p=int(gmpy2.next_prime(random.randrange(lo,2*lo-1)))
        if p<2**64: out.append(p)
out+= [18446744073709551557,18446744073709551533,18446744073709551521]
open('primes_sample.txt','w').write('\n'.join(map(str,out))+'\n'); print(len(out))
