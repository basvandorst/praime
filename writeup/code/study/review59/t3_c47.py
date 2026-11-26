"""Task 3: own C47 labels on primes in (L0, L1]: 0 inert ((p/47) = -1), else the reduced form of the prime form (p, b, c),
b^2 - 4pc = -47 (Tonelli-Shanks + Gauss reduction): 1 = (1,1,12), 2 = (2,+-1,6), 3 = (3,+-1,4).  Lag-1 statistics in the split
fibre: 3x3 table against Chebotarev (1/5, 2/5, 2/5), chi2 (4 dof), and r = mean f(p_n) g(p_n+1), f, g in {f0, f1},
f0 = 2cos(2 pi k/5)/sqrt2, f1 = 2cos(4 pi k/5)/sqrt2 with k = 0, 1, 2 for labels 1, 2, 3 (g = (2,1,6) taken as generator).
usage: t3_c47.py L0 L1"""
import sys, math
import numpy as np
L0, L1 = int(float(sys.argv[1])), int(float(sys.argv[2]))
D = -47
small = np.ones(31700, bool); small[:2] = False
for q in range(2, 179):
    if small[q]: small[q * q::q] = False
SP = np.nonzero(small)[0]

def primes_in(L, U):           # primes in (L, U]
    s = np.ones(U - L, bool)   # index i <-> L + 1 + i
    for q in SP:
        q = int(q)
        if q * q > U: break
        st = max(q * q, ((L + 1 + q - 1) // q) * q); s[st - L - 1::q] = False
    return np.nonzero(s)[0].astype(np.int64) + L + 1

def sqrt_mod(a, p):
    a %= p
    if p % 4 == 3: return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, tt = 0, t
        while tt != 1: tt = tt * tt % p; i += 1
        b = pow(c, 1 << (m - i - 1), p); m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r

def reduce_form(a, b, c):
    while True:
        if c < a or (c == a and b < 0): a, b, c = c, -b, a
        if not (-a < b <= a):
            k = (a - b) // (2 * a); b, c = b + 2 * k * a, a * k * k + b * k + c   # b -> b + 2ka into (-a, a]
            continue
        if c < a or (c == a and b < 0): continue
        return a, b, c

qr47 = np.zeros(47, bool); qr47[(np.arange(1, 47) ** 2) % 47] = True
lab_all = []
SEG = 10 ** 7
for L in range(L0, L1, SEG):
    pr = primes_in(L, min(L + SEG, L1)); lab = np.zeros(len(pr), np.int8)
    sp = qr47[pr % 47]
    for i in np.nonzero(sp)[0]:
        p = int(pr[i]); b = sqrt_mod(D, p)
        if b % 2 == 0: b = p - b
        c = (b * b - D) // (4 * p); a2, b2, c2 = reduce_form(p, b, c)
        assert b2 * b2 - 4 * a2 * c2 == D and a2 in (1, 2, 3), (p, a2, b2, c2)
        lab[i] = a2
    lab_all.append(lab)
lab = np.concatenate(lab_all); n = len(lab)
cnt = np.bincount(lab, minlength=4)
print('primes in (%.0e, %.0e]: %d; label counts %s; split shares %s (Chebotarev 0.2 0.4 0.4)'
      % (L0, L1, n, cnt.tolist(), np.round(cnt[1:] / cnt[1:].sum(), 5).tolist()))
x, y = lab[:-1], lab[1:]; m = (x > 0) & (y > 0); x, y = x[m], y[m]
T = np.zeros((4, 4), int); np.add.at(T, (x, y), 1); T = T[1:, 1:]
cheb = np.array([0.2, 0.4, 0.4]); E = T.sum() * np.outer(cheb, cheb)
chi2 = float(((T - E) ** 2 / E).sum())
# chi2 against the table's own margins (independence), 4 dof
Ei = np.outer(T.sum(1), T.sum(0)) / T.sum(); chi2i = float(((T - Ei) ** 2 / Ei).sum())
print('lag 1 split-split pairs %d; table %s' % (T.sum(), T.tolist()))
print('chi2 vs Chebotarev product %.2f (8 dof incl. margins); vs own margins %.2f (4 dof)' % (chi2, chi2i))
k = np.array([0, 0, 1, 2])
f0 = 2 * np.cos(2 * np.pi * k / 5) / math.sqrt(2); f1 = 2 * np.cos(4 * np.pi * k / 5) / math.sqrt(2)
for nm, F, G in (('f0f0', f0, f0), ('f0f1', f0, f1), ('f1f0', f1, f0), ('f1f1', f1, f1)):
    r = float(np.mean(F[x] * G[y])); print('r %s = %+.6f  z %+.2f  (sigma_r %.2e, 5 sigma = %.2e)' % (nm, r, r * math.sqrt(len(x)), 1 / math.sqrt(len(x)), 5 / math.sqrt(len(x))))
