import numpy as np
exec(open('/home/bas/scratch/review-logic5/r129.py').read().split("# check integral")[0])
x = Nsm(T)
Lm = np.mean(np.log(T/(2*np.pi)))
P0 = primes.copy()
ls = np.arange(0.5, 14.01, 0.5)
D = []; B = []; Bno2 = []; G = []
for l in ls:
    D.append(sigma2(x, l, 200000)[0]); B.append(berry(l, Lm)); G.append(gue(l))
    primes = P0[1:]; Bno2.append(berry(l, Lm)); primes = P0
D, B, Bno2, G = map(np.array, (D, B, Bno2, G))
for row in zip(ls, D, B, D-B, Bno2, G):
    print('l=%.1f data=%.4f berry=%.4f diff=%+.4f berry_no2=%.4f gue=%.4f' % row)
m = ls <= 8
print('max |D-B| for l<=8:', np.abs(D-B)[m].max(), 'argmax data', ls[D.argmax()], 'argmax berry', ls[B.argmax()], 'argmax no2', ls[Bno2.argmax()])
