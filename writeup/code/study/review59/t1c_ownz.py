"""Task 1(c): own matched-filter statistic for round 157 on 5 curves, by direct sums at the exact line positions.
z_naive = T / (sigma sqrt(sum w_r^2)) as registered (lines treated as independent);
z_cal   = T / sqrt(Var T) with the exact variance under random signs of lambda(p): Var = sum_p W_p^2 Re(ph G_p)^2,
          G_p = sum_r w_r e(-alpha_r sqrt p).  Also a permutation null (lambda permuted over primes, 20 draws).
Weights from the definitions (lines_lib.py docstring, spectrometer_blind.py header), written anew.
usage: t1c_ownz.py NAME[:Nalt] ..."""
import json, math, sys
import numpy as np
import cypari2
pari = cypari2.Pari(); pari.allocatemem(400 * 10 ** 6)
CUR = json.load(open('/home/bas/praime/writeup/code/study/resonance/spectro/testset157.json'))
OLD = json.load(open('/home/bas/praime/writeup/data/study/r157/out_blind157.json'))
X = 3 * 10 ** 6; P0 = 10 ** 4; AMAX = 3.0; AMAXNUM = 1500

def fac(n):
    out = {}; q = 2
    while q * q <= n:
        while n % q == 0: out[q] = out.get(q, 0) + 1; n //= q
        q += 1
    if n > 1: out[n] = out.get(n, 0) + 1
    return out

def weights(N, lam_small):
    bad = fac(N)
    def lpow(l, k):                    # good-prime Hecke lambda(l^k)
        t = [1.0, lam_small[l]]
        for j in range(1, k + 1): t.append(lam_small[l] * t[-1] - t[-2])
        return t[k]
    def F(l, ea, eb):
        lp = lam_small[l]
        if l in bad:
            if lp == 0.0: return 1.0 if ea == eb == 0 else 0.0           # additive
            if eb == 0: return lp ** ea
            if eb == 1: return (-lp / (l - 1.0)) / (1 + 1.0 / l)
            return 0.0
        L1 = 1 - lpow(l, 2) / l
        if eb == 0: Lv = lpow(l, ea) - lpow(l, ea + 2) / l
        elif eb == 1: Lv = -lp / l
        elif eb == 2: Lv = 1.0 / (l * (l - 1.0))
        else: Lv = 0.0
        return Lv / L1
    al, w = [], []
    for b in (1, 2, 3, 4):
        for a in range(1, AMAXNUM + 1):
            if math.gcd(a, b) != 1: continue
            r = a / b; alpha = 2 * math.sqrt(r / N)
            if alpha > AMAX: continue
            f = r ** -0.25
            for l in set(fac(a)) | set(fac(b)):
                ea = fac(a).get(l, 0); eb = fac(b).get(l, 0); f *= F(l, ea, eb)
            if f != 0.0: al.append(alpha); w.append(f)
    return np.array(al), np.array(w)

ph = complex(math.cos(3 * math.pi / 4), math.sin(3 * math.pi / 4))
for arg in sys.argv[1:]:
    name, *alt = arg.split(':'); ai = CUR[name]
    E = pari.ellinit(ai); Nt = int(pari.ellglobalred(E)[0]); et = int(pari.ellrootno(E))
    v = pari('my(E=ellinit(%s), Pr=primes([2,%d])); [Vecsmall(Pr), Vecsmall(apply(p->ellap(E,p), Pr))]' % (ai, X))
    p = np.array(list(v[0]), np.int64); a = np.array(list(v[1]), np.float64); lam = a / np.sqrt(p)
    small = {int(q): float(x) for q, x in zip(p[:2000], lam[:2000])}
    m = p >= P0; u = np.sqrt(p[m].astype(np.float64)); u0, u1 = math.sqrt(P0), math.sqrt(X)
    t = (u - u0) / (u1 - u0); hann = np.where((t > 0) & (t < 1), np.sin(np.pi * t) ** 2, 0.0)
    W = lam[m] * np.log(p[m].astype(np.float64)) * hann
    sig = math.sqrt(float(W @ W) / 2)
    rng = np.random.default_rng(1)
    for N in [Nt] + [int(x) for x in alt]:
        al, w = weights(N, small)
        G = np.zeros(len(u), complex)
        for s in range(0, len(u), 1500):
            G[s:s + 1500] = np.exp(-2j * math.pi * np.outer(u[s:s + 1500], al)) @ w
        RG = (ph * G).real
        for eps in ([et] if N == Nt else [1, -1]):
            T = eps * float(W @ RG)
            zn = T / (sig * math.sqrt(float(w @ w)))
            zc = T / math.sqrt(float((W * W) @ (RG * RG)))
            perm = [eps * float(rng.permutation(W) @ RG) for _ in range(20)]
            zp = T / float(np.std(perm))
            old = OLD[name]['z_true'] if N == Nt else None
            print('%-9s N=%4d eps=%+d lines=%4d  z_naive %.2f (stored %s)  z_cal %.2f  z_perm %.2f  ratio naive/cal %.2f'
                  % (name, N, eps, len(w), zn, ('%.2f' % old) if old else '-', zc, zp, zn / zc), flush=True)
