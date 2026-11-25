"""validate_sieved.py: the resonance on sieved INTEGERS (no primes), against the local-factor formula of spectro.md, section 1.

S(type, P, r) = sum_{X0 <= n < X1} lambda(n) F_P(n) n^{-1/4} W(sqrt n) e(-2 sqrt(r n / N)),  F_P = prod_{l <= P} f_l(n mod l)
types:  int (F = 1), prime (f_l = l/(l-1) [n != 0]), twin (f_l = (l/(l-1))^2 [n != 0, -2], f_2 = 4 [n odd]),
        shift2 (n - 2 'prime': f_l = l/(l-1) [n != 2], f_2 = 2 [n odd]), and residue classes n = c mod l (f = l [n = c]).
Prediction: S = c_pred * MT(r), MT(r) = (3/2) K1 r^{-1/4} int W (1 - i kappa sqrt N/(2 pi v sqrt r)) dv,
  c_pred = prod_{l <= P} L_l(r) * prod_{l > P} lambda(l^{v_l(a)}),  r = a/b (b | P# squarefree, else 0)
  L_l(r) = mean_l * lambda(l^{v_l(a)}) [l !| b] + (f_l(0) - mean_l) * chi(l) lambda(l^{v_l(a) + 2 - v_l(b)}) / l   (good l)
  bad l (l !| ab): L_l = mean_l (1 - rho_l) + rho_l f_l(0).
usage: validate_sieved.py LABEL X1      (window n in [1e5, X1), X1 <= 3e7)
"""
import sys, os, json, math, time
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import spectro_lib as S

label = sys.argv[1]; X1 = int(float(sys.argv[2])); X0 = 100000
assert X1 <= 30_000_000
I = S.label_info(label); N = I['N']
v0, v1 = math.sqrt(X0), math.sqrt(X1)
LS = [int(p) for p in S.primes_upto(100)]
PS = [2, 10, 30, 100]
RS = [Fr(1), Fr(1, 2), Fr(2), Fr(3), Fr(1, 3), Fr(2, 3), Fr(3, 2), Fr(5), Fr(1, 5), Fr(6), Fr(1, 6)]
CLASS_L = [2, 3, 5, 7, 11, 13] + ([23] if N == 23 else [])
TYPES = ['prime', 'twin', 'shift2']

def ftab(t, l):
    x = np.arange(l)
    if t == 'prime': return (l / (l - 1)) * (x != 0)
    if t == 'twin':
        if l == 2: return 4.0 * (x == 1)
        return (l / (l - 1)) ** 2 * ((x != 0) & (x != (l - 2) % l))
    if t == 'shift2':
        if l == 2: return 2.0 * (x == 1)
        return (l / (l - 1)) * (x != 2 % l)
    raise ValueError(t)

# ---------------------------------------------------------------- label on segments
if label.startswith('EC:'):
    ptab, atab = S.ap_table(label[3:], 30_000_000)
    lamseg = lambda L, U: S.lam_ec_seg(L, U, ptab, atab, N)
elif label.startswith('A'):
    lamseg = lambda L, U: S.lam_gauss_seg(L, U, int(label[1:]))
else:
    lamseg = lambda L, U: S.lam_c23_seg(L, U)

# ---------------------------------------------------------------- sums
t0 = time.time()
acc = {}; noise = {}
def add(key, z, nz):
    acc[key] = acc.get(key, 0j) + z; noise[key] = noise.get(key, 0.0) + nz
SEG = 1_000_000
for L in range(X0, X1, SEG):
    U = min(X1, L + SEG)
    n = np.arange(L, U, dtype=np.int64); nf = n.astype(np.float64); v = np.sqrt(nf)
    w0 = lamseg(L, U) * nf ** -0.25 * S.hann(v, v0, v1)
    res = {l: (n % l).astype(np.int64) for l in LS}
    Fs = {}
    for t in TYPES:
        F = np.ones(U - L)
        for l in LS:
            F *= ftab(t, l)[res[l]]
            if l in PS or (l == 97 and 100 in PS):
                pass
            for P in PS:
                if l == max(q for q in LS if q <= P): Fs[(t, P)] = F.copy()
    for r in RS:
        E = np.exp(-1j * S.TWO_PI * 2 * np.sqrt(float(r) * nf / N))
        wE = w0 * E
        add(('int', 0, str(r)), complex(wE.sum()), float(np.dot(w0, w0)))
        for (t, P), F in Fs.items():
            add((t, P, str(r)), complex(np.dot(F, wE)), float(np.dot(w0 * F, w0 * F)))
        if r == 1:
            for l in CLASS_L:
                re = np.bincount(res[l], weights=wE.real, minlength=l); im = np.bincount(res[l], weights=wE.imag, minlength=l)
                nz = np.bincount(res[l], weights=w0 * w0, minlength=l)
                for c in range(l): add(('class%d' % l, c, '1'), l * complex(re[c], im[c]), l * l * float(nz[c]))
    print('segment', L, '%.0fs' % (time.time() - t0), flush=True)

# ---------------------------------------------------------------- predictions
loc = S.local_lams(label, 200)
def lam_pk(l, k):
    lp, chi = loc[l]; return S.lam_pow_table(lp, chi if chi != 0 else 0.0, max(k, 2) + 2)[k]
def vval(x, l):
    k = 0
    while x % l == 0: x //= l; k += 1
    return k
def c_pred(t, P, r):
    a, b = r.numerator, r.denominator
    if t == 'int':
        if b != 1: return 0.0
        out = 1.0
        for l in S.primes_upto(max(a, 2)):
            l = int(l); e = vval(a, l)
            if e: out *= lam_pk(l, e)
        return out
    if t.startswith('class'):
        l = int(t[5:]); c = P
        f = np.zeros(l); f[c] = l; tabs = {l: f}; Ps = [l]
    else:
        tabs = {l: ftab(t, l) for l in LS if l <= P}; Ps = sorted(tabs)
    bb = b
    for l in Ps:
        if bb % l == 0: bb //= l
    if bb != 1: return 0.0
    out = 1.0
    for l in Ps:
        f = tabs[l]; mean = f.mean(); f0 = f[0]; lp, chi = loc[l]
        ea, eb = vval(a, l), vval(b, l)
        good = chi != 0.0
        if good:
            term1 = mean * lam_pk(l, ea) if eb == 0 else 0.0
            term2 = (f0 - mean) * chi * lam_pk(l, ea + 2 - eb) / l
            out *= term1 + term2
        else:
            if ea == 0 and eb == 0:
                out *= S.local_factor(label, l, lp, chi, f)
            else:
                return float('nan')
    for l in S.primes_upto(max(a, 2)):
        l = int(l)
        if l in Ps: continue
        e = vval(a, l)
        if e: out *= lam_pk(l, e)
    return out

rows = []
for key in sorted(acc, key=str):
    t, P, rs = key; r = Fr(rs)
    MT = S.main_term_unit(label, v0, v1, float(r))
    cp = c_pred(t, P, r)
    z = acc[key] / MT; err = math.sqrt(noise[key] / 2) / abs(MT)
    rows.append(dict(type=t, P=P, r=rs, read=[z.real, z.imag], pred=cp, err_random=err))
    print('%-8s P=%-3d r=%-4s read %+.4f %+.4fi  pred %+.4f  (random-model err %.4f)  dev/err %+.1f' % (t, P, rs, z.real, z.imag, cp, err, (z.real - cp) / err if cp == cp else float('nan')))
json.dump(dict(label=label, X0=X0, X1=X1, rows=rows), open(os.path.join(S.HERE, 'out_validate_%s_%d.json' % (label.replace(':', '_'), X1)), 'w'), indent=0)
print('done %.0fs' % (time.time() - t0))
