"""Did the range runs test every odd composite? Compare the 'tested' count of each log with
(odd numbers in the range) - (primes in the range), the primes counted by PARI's primepi."""
import re, os
import cypari2
pari=cypari2.Pari()
runs=[('single_doubling_1e8.txt',1,10**8),('single_fib_1e8.txt',1,10**8),('single_perrin_1e8.txt',1,10**8),('single_ring5_1e8.txt',1,10**8),
      ('single_frob_1e9.txt',1,10**9),('single_frob_1e10.txt',10**9,10**10),('single_frob_2e10a.txt',10**10,15*10**9),('single_frob_2e10b.txt',15*10**9,2*10**10)]
for fn,lo,hi in runs:
    if not os.path.exists(fn) or 'tested' not in open(fn).read(): continue
    t=int(re.search(r'tested (\d+) liars (\d+)',open(fn).read()).group(1))
    a=max(lo,3); odd=(hi+1)//2-a//2
    pr=int(pari.primepi(hi))-int(pari.primepi(a-1))
    print("%-26s range %d..%d: tested %d; odd numbers %d - primes %d = %d; equal: %s"%(fn,lo,hi,t,odd,pr,odd-pr,t==odd-pr))
