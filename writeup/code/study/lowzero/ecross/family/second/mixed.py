"""Mixed moments X_q = sum^h eps_f lambda_f(q^2) Z_f (q squarefree, (q,N)=1), C2(Q) = sum^h c_f^(Q)^2, Delta(1,1);
then V(Q) = M11 + 2 L Re X_Q + L^2 C2(Q),  X_Q = sum_{q<=Q} mu(q)/q X_q.
usage: mixed.py N k T Q MG model [lo hi] [M11 from moment npz files, comma separated]"""
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+"/")
import math, time, numpy as np
from smlib import *
from scipy.special import jv
N, k, T, Q, MG, model = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
lo_fac = float(sys.argv[7]) if len(sys.argv) > 7 else 0.3
hi_fac = float(sys.argv[8]) if len(sys.argv) > 8 else 3.5
npzs = sys.argv[9].split(',') if len(sys.argv) > 9 else []
t0 = time.time()
win = Window(N, k, T, h=0.05)
nlo = int(N*(lo_fac*T/TWO_PI)**2); nhi = int(N*(hi_fac*T/TWO_PI)**2)
n, b = prime_power_coeffs(win, nlo, nhi, model)
def mobius(x):
    r, m, p = 1, x, 2
    while p*p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0: return 0
            r = -r
        p += 1
    return -r if m > 1 else r
qs = [q for q in range(1, Q+1) if mobius(q) != 0 and q % N != 0]
inv_cache = {}
def G_vec(ell_arr, coef):
    """sum_i coef_i G(ell_i), G(l) = 2 pi N^{-1/2} sum_{m<=MG,(m,N)=1} S(l, Nbar; m)/m J_{k-1}(4 pi sqrt(l/N)/m)."""
    ell_arr = np.asarray(ell_arr, np.int64); tot = 0j
    y0 = 4*math.pi*np.sqrt(ell_arr/N)
    for m in range(1, MG+1):
        if m % N == 0:
            if N >= 1000: continue
            # general modulus c = N m: G-term -2 pi N^{-1/2} S(N l, 1; c)/m J (direct, Lemma 2 before Lemma 3(i))
            c = N*m; w = coef*jv(k-1, y0/m)
            if c not in inv_cache: inv_cache[c] = inv_table(c)
            inv = inv_cache[c]; x = np.nonzero(inv)[0]
            A = np.bincount((N*ell_arr) % c, weights=w.real, minlength=c) + 1j*np.bincount((N*ell_arr) % c, weights=w.imag, minlength=c)
            Ah = np.fft.ifft(A)*c
            tot += -np.sum(Ah[x]*np.exp(2j*math.pi*inv[x]/c))/m
            continue
        w = coef*jv(k-1, y0/m)
        if m == 1:
            tot += np.sum(w); continue
        if m not in inv_cache: inv_cache[m] = inv_table(m)
        inv = inv_cache[m]; x = np.nonzero(inv)[0]
        A = np.bincount(ell_arr % m, weights=w.real, minlength=m) + 1j*np.bincount(ell_arr % m, weights=w.imag, minlength=m)
        Ah = np.fft.ifft(A)*m                                   # sum_b A(b) e(bx/m)
        Nb = pow(N, -1, m)
        tot += np.sum(Ah[x]*np.exp(2j*math.pi*((Nb*inv[x]) % m)/m))/m
    return 2*math.pi/math.sqrt(N)*tot
def Hecke_terms(q2, nn):
    """lambda(q^2) lambda(nn) = sum_{d | (q^2, nn)} lambda(q^2 nn/d^2)  (all (d,N)=1 here)."""
    g = math.gcd(q2, int(nn)); return [q2*int(nn)//(d*d) for d in range(1, g+1) if g % d == 0]
def Hh(xi): return np.sum(win.H*np.exp(-1j*xi*win.tt))*win.h
def ram(m, l): return sum(d*mobius(m//d) for d in range(1, m+1) if m % d == 0 and l % d == 0)
def phi(m): return sum(1 for a in range(1, m+1) if math.gcd(a, m) == 1)
print('N=%d k=%d T=%g model=%s: %d indices; L = %.5f; Q = %d, MG = %d' % (N, k, T, model, len(n), win.L, Q, MG), flush=True)
X = {}
jobs = {}
for q in qs:
    q2 = q*q; ells = []; cfs = []
    for v, bv in zip(n, b):
        if math.gcd(q2, int(v)) == 1:
            ells.append(q2*int(v)); cfs.append(bv)
        else:
            for e in Hecke_terms(q2, v): ells.append(e); cfs.append(bv)
    jobs[q] = (np.array(ells, np.int64), np.array(cfs, complex))
    X[q] = 0j
t1 = time.time()
for m in range(1, MG+1):
    if m % N == 0 and N >= 1000: continue
    gen = (m % N == 0)
    c = N*m if gen else m
    inv = inv_table(c) if c > 1 else None
    x = np.nonzero(inv)[0] if c > 1 else None
    if not gen and m > 1:
        ph = np.exp(2j*math.pi*((pow(N, -1, m)*inv[x]) % m)/m)
    elif gen:
        ph = np.exp(2j*math.pi*inv[x]/c)
    for q in qs:
        ell_arr, coef = jobs[q]
        w = coef*jv(k-1, 4*math.pi*np.sqrt(ell_arr/N)/m)
        if m == 1:
            X[q] += 2*math.pi/math.sqrt(N)*np.sum(w); continue
        key = (N*ell_arr) % c if gen else ell_arr % c
        A = np.bincount(key, weights=w.real, minlength=c) + 1j*np.bincount(key, weights=w.imag, minlength=c)
        Ah = np.fft.ifft(A)*c
        X[q] += (-1 if gen else 1)*2*math.pi/math.sqrt(N)*np.sum(Ah[x]*ph)/m
for q in qs:
    val = X[q]
    lines = -sum(mobius(mm)*ram(mm, q*q)/(phi(mm)*q)*Hh(math.log(mm*mm/(q*q))) for mm in range(1, 41) if mobius(mm) != 0)/TWO_PI
    print('q=%3d: X_q = sum^h eps lam(q^2) Z = %+.5f%+.5fi ; Theorem D main -L mu(q)/q = %+.5f ; uniform line sum %+.5f%+.5fi ; diff %+.5f%+.5fi'
          % (q, val.real, val.imag, -win.L*mobius(q)/q, lines.real, lines.imag, (val-lines).real, (val-lines).imag), flush=True)
print('G sums: %.0fs' % (time.time()-t1))
# C2(Q) via Petersson Delta(q^2, q'^2); Delta(1,1)
def Delta(a, bb, mmax=60):
    s = 0.0
    for m in range(1, mmax+1):
        if m % N == 0: continue
        c = N*m; s += kloost(a, bb, c)/c*jv(k-1, 4*math.pi*math.sqrt(a*bb)/c)
    return (1.0 if a == bb else 0.0) + 2*math.pi*(-1)**(k//2)*s
D11 = Delta(1, 1)
print('Delta(1,1) = sum^h 1 = %.10f' % D11)
res = []
for Qc in sorted(set([1, 2, 3, 5, 6, 10, 15, 20, 30, Q])):
    if Qc > Q: continue
    qq = [q for q in qs if q <= Qc]
    XQ = sum(mobius(q)/q*X[q] for q in qq)
    C2 = 0.0
    for i, q in enumerate(qq):
        for q_ in qq:
            if q_ < q: continue
            d = Delta(q*q, q_*q_, 6) if q_ != q or q > 1 else D11
            C2 += (1 if q == q_ else 2)*mobius(q)*mobius(q_)/(q*q_)*d
    res.append((Qc, XQ, C2))
    print('Q=%3d: X_Q = %+.5f%+.5fi (-L C2 = %+.5f) ; C2(Q) = %.5f (sum mu^2/q^2 = %.5f)' % (Qc, XQ.real, XQ.imag, -win.L*C2, C2, sum(1/q**2 for q in qq)), flush=True)
if npzs:
    M11 = 0j
    for f in npzs:
        z = np.load(f, allow_pickle=True)
        if int(z['rows'][0, 0]) == 1:
            M11 += complex(z['D11']) + abs(complex(z['b1']))**2 + 2*(np.conj(complex(z['b1']))*complex(z['G1'])).real
        M11 += np.sum(z['rows'][:, 3])
    print('M11 = %.5f (from %s)' % (M11.real, npzs))
    for (Qc, XQ, C2) in res:
        V = M11.real + 2*win.L*XQ.real + win.L**2*C2
        print('Q=%3d: V = M11 + 2L Re X_Q + L^2 C2 = %.5f ; V/T^2 = %.5f ; V/D11 = %.4f' % (Qc, V, V/T**2, V/float(np.load(npzs[0])['D11'])))
np.savez('mixed_%d_%d_%g_%s.npz' % (N, k, T, model), qs=np.array(qs), Xq=np.array([X[q] for q in qs]), res_Q=np.array([r[0] for r in res]),
         res_XQ=np.array([r[1] for r in res]), res_C2=np.array([r[2] for r in res]), L=win.L, D11=D11)
print('total %.0fs' % (time.time()-t0))
