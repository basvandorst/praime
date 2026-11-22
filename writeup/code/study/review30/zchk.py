import mpmath as mp, sys
mp.mp.dps=25
Z=[l.strip() for l in open('/home/bas/praime/writeup/data/riemann/zeros_9999_flint.txt') if l.strip()]
print(len(Z))
for n in (1,1000,5000,9999,10000):
    z=mp.zetazero(n).imag
    print(n, Z[n-1], z, float(mp.mpf(Z[n-1])-z))
# N(T) check at T just after last
T=mp.mpf(Z[-1])
print('N(T+0.01) via backlund', mp.nzeros(T+mp.mpf('0.01')))
