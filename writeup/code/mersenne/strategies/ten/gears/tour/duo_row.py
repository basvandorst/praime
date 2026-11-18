"""The three results of the two-cell wheel for 250 wheels (lamps stay on), in plain Python."""
import json
from sympy import isprime
from tour_check2_lib import home
MX={2,3,5,7,13,17,19,31,61,89,107,127}
lit=[2]; best={'p':(1,2),'m':(1,2),'w':(0,0)}
for t in range(3,251):
    if t%2 and home(t): lit.append(t)
    p=[k for k in lit if isprime(k)]; m=[k for k in p if k in MX]; w=[k for k in p if isprime(k-2) or isprime(k+2)]
    for key,a in (('p',p),('m',m),('w',w)):
        if len(a)>best[key][0]: best[key]=(len(a),t)
row=list(best['p']+best['m']+best['w'])
print("two-cell wheel, 250 wheels: primes %d (tick %d), Mersenne exponents %d (tick %d), twin primes %d (tick %d); lamps of numbers that are not prime: %d"%(tuple(row)+(sum(1 for k in lit if not isprime(k)),)))
json.dump({'row250':row},open('duo_row.json','w'))
