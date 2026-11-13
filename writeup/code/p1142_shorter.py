"""Prime-representing constants f_{n+1} = floor(f_n) + R(floor(f_n)) * frac(f_n).

If every prime p is followed by a gap g(p) <= R(p) - 1, then
  f_1 = 2 + sum_{n>=1} g_n / (R(p_1) ... R(p_n))
gives floor(f_n) = p_n for every n (frac(f_n) = (g_n + frac(f_{n+1}))/R(p_n) < 1).
R(p) = p is Fridman, Garbulsky, Glecer, Grime, Tron Florentin, Amer. Math. Monthly 126 (2019) (Bertrand).
New: R(p) = floor(ln^2 p) + 4 for p < 4e18 (every gap there is below ln^2 p + 2: Firoozbakht's conjecture,
which implies g < ln^2 p - ln p for p > 11, is verified below 4e18, Kourbatov 2015; small p checked here), and
R(p) = floor(p / (25 ln^2 p)) + 1 for p >= 4e18 (Dusart 2010: a prime in (x, x(1 + 1/(25 ln^2 x))] for x >= 396738).
Measured: primes recovered from the first d digits, and digits needed for the first n primes.
"""
import math, sys
import mpmath as mp
import gmpy2

def R_fridman(p): return p
def R_new(p):
    if p < 4 * 10 ** 18:
        return int(math.floor(math.log(p) ** 2)) + 4
    return int(p / (25 * math.log(p) ** 2)) + 1

N = 2000
P = [2]
while len(P) < N + 1:
    P.append(int(gmpy2.next_prime(P[-1])))
gaps = [P[i + 1] - P[i] for i in range(N)]
for p, g in zip(P, gaps):          # the bound must hold for every prime used
    assert g + 1 <= R_new(p), (p, g, R_new(p))

def constant(Rf, dps):
    mp.mp.dps = dps
    f = mp.mpf(0); denom = mp.mpf(1)
    for p, g in zip(P, gaps):
        denom *= Rf(p); f += g / denom
        if denom > mp.mpf(10) ** (dps + 5): break
    return 2 + f

def decode(fd, Rf, dps):
    mp.mp.dps = dps
    out = []; f = fd
    for k in range(N):
        p = int(mp.floor(f))
        if p != P[k]: break
        out.append(p)
        f = p + Rf(p) * (f - p)
    return out

def R_half(p):          # radix for half-gaps of odd primes p (the gap after an odd prime is even)
    if p < 4 * 10 ** 18:
        return int(math.floor(math.log(p) ** 2 / 2)) + 3
    return int(p / (50 * math.log(p) ** 2)) + 1
for p, g in zip(P[1:], gaps[1:]):
    assert g // 2 + 1 <= R_half(p), (p, g)

def constant_half(dps):
    mp.mp.dps = dps
    f = mp.mpf(0); denom = mp.mpf(1)
    for p, g in zip(P[1:], gaps[1:]):
        denom *= R_half(p); f += (g // 2) / denom
        if denom > mp.mpf(10) ** (dps + 5): break
    return 1 + f                       # u_1 = (3 - 1)/2 + ...

def decode_half(ud, dps):
    mp.mp.dps = dps
    out = [2]; u = ud
    for k in range(1, N):
        m = int(mp.floor(u)); p = 2 * m + 1
        if p != P[k]: break
        out.append(p)
        u = m + R_half(p) * (u - m)
    return out

lines = []
for name, Rf in (('Fridman, R(p) = p', R_fridman), ('new, R(p) = floor(ln^2 p) + 4', R_new)):
    full = constant(Rf, 1200)
    lines.append(f'{name}: f_1 = {mp.nstr(full, 40)}')
    row = []
    for d in (20, 50, 100, 200, 1000):
        mp.mp.dps = 1200
        fd = mp.floor(full * mp.mpf(10) ** (d - 1)) / mp.mpf(10) ** (d - 1)     # first d significant digits
        got = decode(fd, Rf, 1200)
        row.append(f'{d} digits -> {len(got)} primes (up to {got[-1]})')
    lines.append('   ' + '; '.join(row))

full = constant_half(1200)
lines.append(f'new, odd primes: u_(n+1) = floor(u_n) + R(2 floor(u_n) + 1) frac(u_n), R(p) = floor(ln^2 p / 2) + 3, p_n = 2 floor(u_n) + 1 (2 given): u_1 = {mp.nstr(full, 40)}')
row = []
for d in (20, 50, 100, 200, 1000):
    mp.mp.dps = 1200
    ud = mp.floor(full * mp.mpf(10) ** (d - 1)) / mp.mpf(10) ** (d - 1)
    got = decode_half(ud, 1200)
    row.append(f'{d} digits -> {len(got)} primes (up to {got[-1]})')
lines.append('   ' + '; '.join(row))

# digits needed for the first n primes: log10 of the product of the radices
lines.append('digits needed for the first n primes (log10 of the product of radices), and the prime constant sum 2^-p (one bit per integer):')
p = 2; cnt = 0; sF = 0.0; sN = 0.0; sH = 0.0; marks = {10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6}
while cnt < 10 ** 6:
    q = int(gmpy2.next_prime(p)); g = q - p
    assert g + 1 <= R_new(p) and (p == 2 or g // 2 + 1 <= R_half(p))
    sF += math.log10(p); sN += math.log10(R_new(p)); cnt += 1
    if p > 2: sH += math.log10(R_half(p))
    if cnt in marks:
        lines.append(f'   first {cnt} primes (up to {p}): Fridman {sF:,.0f} digits, new {sN:,.0f} ({sF / sN:.2f}x fewer), odd variant {sH:,.0f} ({sF / sH:.2f}x fewer), prime constant {p * math.log10(2):,.0f} digits')
    p = q
print('\n'.join(lines))
open('p1142_shorter.txt', 'w').write('\n'.join(lines) + '\n')
