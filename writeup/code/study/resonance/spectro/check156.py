# Independent check of one reading of round 156: the resonance along twin primes for the Hecke character of Q(i).
# S = sum over twin primes (p, p+2), 3e7 < p <= 1e9, of a(p) ln p ln(p+2) p^(-1/4) W(sqrt p) e(-sqrt p),
# a(p) = 2 cos(4 theta_p) for p = a^2 + b^2 (0 for p = 3 mod 4), W the Hann window in v = sqrt p on [sqrt 3e7, sqrt 1e9].
# Law: S = C e^(-i pi/4) int W dv, C = 2 C2 c, c the constant of the single-prime law (round 153; from the registered file).
# Own code: segmented sieve (numpy), two squares by PARI. Reading = S / (C e^(-i pi/4) (V1 - V0)/2).
import numpy as np, math, cmath, sys, json, time
import cypari2
pari = cypari2.Pari()
X0, X1 = 30_000_000, 1_000_000_000
V0, V1 = math.sqrt(X0), math.sqrt(X1)
base = np.ones(31700, dtype=bool); base[:2] = False
for i in range(2, 179):
    if base[i]: base[i*i::i] = False
small = np.nonzero(base)[0]
def primes_seg(lo, hi):            # primes in [lo, hi)
    s = np.ones(hi - lo, dtype=bool)
    for q in small:
        if q * q >= hi: break
        st = max(q * q, (lo + q - 1) // q * q)
        s[st - lo::q] = False
    return np.nonzero(s)[0] + lo
t0 = time.time(); S = 0j; ntw = 0; n1 = 0; S_single = 0j
SEG = 10_000_000
twos = pari('(p) -> my(v = qfbsolve(Qfb(1,0,1), p)); [abs(v[1]), abs(v[2])]')
for lo in range(X0, X1, SEG):
    pr = primes_seg(lo + 1, lo + SEG + 3)          # primes in (lo, lo+SEG+2]
    isp = np.zeros(SEG + 4, dtype=bool); isp[pr - lo] = True
    p = pr[pr <= lo + SEG]
    tw = p[isp[p - lo + 2]]
    ntw += len(tw)
    t1 = tw[tw % 4 == 1]; n1 += len(t1)
    ab = np.array([[int(x) for x in twos(int(q))] for q in t1], dtype=np.float64)
    th = np.arctan2(ab[:, 1], ab[:, 0])
    a = 2 * np.cos(4 * th)
    v = np.sqrt(t1.astype(np.float64))
    W = np.sin(np.pi * (v - V0) / (V1 - V0)) ** 2
    w = a * np.log(t1) * np.log(t1 + 2.0) * t1 ** -0.25 * W
    S += np.sum(w * np.exp(-2j * np.pi * (v % 1.0)))
print('twin pairs with 3e7 < p <= 1e9:', ntw, ' with p = 1 mod 4:', n1, ' time %.0f s' % (time.time() - t0))
pred = json.load(open(sys.argv[1]))
C2 = 0.6601618158468696
c = 2.026446079853554           # single-prime constant of the label A4 as printed by the discovery run
C = 2 * C2 * c
MT = C * cmath.exp(-1j * math.pi / 4) * (V1 - V0) / 2
r = S / MT
print('S = %.4f %+.4fi   law = %.4f %+.4fi' % (S.real, S.imag, MT.real, MT.imag))
print('reading = %.4f %+.4fi   (registered run: +0.9370 +0.0295i)' % (r.real, r.imag))
ok1 = ntw == 3271614
ok2 = abs(r - complex(0.9370, 0.0295)) < 2e-3
print('%d checks, %d failed' % (2, (not ok1) + (not ok2)))
