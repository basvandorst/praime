"""Form side at a small prime level: every newform of S_k(Gamma_0(N)) from PARI (mfeigenbasis, all embeddings), harmonic
weight omega_f from the Petersson norm (mfpetersson), root number eps_f = -i^k lambda_f(N) sqrt N, and
Z_f(H) = sum_n beta_n lambda_f(n) (Proposition 1 with the same beta as the Petersson side, smlib.prime_power_coeffs).
Outputs sum_f omega_f |Z_f|^2, sum omega Z_f^2, X_q = sum omega eps lambda(q^2) Z, C2(Q), V(Q), and sum omega c_f^2 (c_f = B_f/L(1,Sym^2 f)
via omega_f: L(1,Sym^2 f) = 2 pi^2/((k-1) N omega_f)).   usage: forms_side.py N k T [lo hi]"""
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+"/")
import math, time, numpy as np, cypari2
from smlib import *
N, k, T = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
lo_fac = float(sys.argv[4]) if len(sys.argv) > 4 else 0.3
hi_fac = float(sys.argv[5]) if len(sys.argv) > 5 else 3.5
pari = cypari2.Pari(); pari.allocatemem(900*10**6); pari.set_real_precision(38)
t0 = time.time()
win = Window(N, k, T, h=0.05)
nlo = int(N*(lo_fac*T/TWO_PI)**2); nhi = int(N*(hi_fac*T/TWO_PI)**2)
n, b = prime_power_coeffs(win, nlo, nhi, 'primes')
Q = 10
P = max(int(n.max()), Q*Q, 5000)
pari('mf=mfinit([%d,%d],0); B=mfeigenbasis(mf);' % (N, k))
forms = []
for j in range(1, int(pari('#B'))+1):
    co = pari('mfembed(B[%d], mfcoefs(B[%d], %d))' % (j, j, P))
    pet = pari('mfpetersson(mfsymbol(mf,B[%d]))' % j) if N <= 61 else None
    isv = str(pari('type(mfembed(B[%d], mfcoefs(B[%d],2))[1])' % (j, j))) == 't_VEC'
    if pet is None: embs = [(co[i], 1.0) for i in range(len(co))] if isv else [(co, 1.0)]
    else: embs = [(co[i], float(pet[i, i].real())) for i in range(len(co))] if isv else [(co, float(pet.real()))]
    for (c, pv) in embs:
        a = np.array([float(complex(x).real) for x in c])
        lam = np.zeros(len(a)); lam[1:] = a[1:]/np.arange(1, len(a))**((k-1)/2)
        om = math.gamma(k-1)/((4*math.pi)**(k-1)*(N+1)*pv)
        eps = -(1 if k % 4 == 0 else -1)*lam[N]*math.sqrt(N)
        forms.append(dict(lam=lam, om=om, eps=eps))
if N > 61:   # omega_f by least squares on sum_f omega_f lambda_f(n) = Delta(1,n), n <= 60 (validated against mfpetersson at N = 61 by the referee)
    from scipy.special import jv
    ns = np.arange(1, 61); rhs = np.zeros(60)
    for c in range(N, N*300 + 1, N):
        rhs += np.array([kloost(1, int(v), c) for v in ns])/c*jv(k-1, 4*math.pi*np.sqrt(ns)/c)
    rhs = (ns == 1).astype(float) + 2*math.pi*(-1)**(k//2)*rhs
    A = np.array([[f['lam'][v] for f in forms] for v in ns])
    oms = np.linalg.lstsq(A, rhs, rcond=None)[0]
    for f, o in zip(forms, oms): f['om'] = float(o)
    print('omega by least squares: max residual %.1e' % np.max(np.abs(A @ oms - rhs)))
print('N=%d k=%d T=%g: %d forms (coefficients to %d), %d indices; sum omega = %.10f  (%.0fs)' % (N, k, T, len(forms), P, len(n), sum(f['om'] for f in forms), time.time()-t0), flush=True)
L = win.L
def mobius(x):
    r, m, p = 1, x, 2
    while p*p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0: return 0
            r = -r
        p += 1
    return -r if m > 1 else r
qs = [q for q in range(1, Q+1) if mobius(q) != 0]
pr = primes_upto(P); pr = pr[pr != N]
M11 = M2 = 0j; Xq = {q: 0j for q in qs}; C2 = {}; Sc2 = 0.0; Sc = 0.0; eps_round = []
cQ_all = []
for f in forms:
    lam, om, eps = f['lam'], f['om'], f['eps']
    eps_round.append(eps)
    Z = np.sum(b*lam[n])
    M11 += om*abs(Z)**2; M2 += om*Z*Z
    for q in qs: Xq[q] += om*eps*lam[q*q]*Z
    cQ = [sum(mobius(q)/q*lam[q*q] for q in qs if q <= Qc) for Qc in range(1, Q+1)]
    cQ_all.append(cQ)
    # c_f = B_f / L(1,Sym^2 f), L(1,Sym^2) from omega
    l = lam[pr]; e1 = l*l - 1; p = pr.astype(float)
    num = 1 - e1/p; den = 1 - e1/p + e1/p**2 - 1/p**3
    Bf = np.prod(np.sign(num))*math.exp(np.sum(np.log(np.abs(num)) - np.log(den)) - math.log1p(-1/N))
    cf = Bf/(2*math.pi**2/((k-1)*N*om))
    f['cf'] = cf; f['Z'] = Z
    Sc += om*cf; Sc2 += om*cf*cf
om = np.array([f['om'] for f in forms]); cQ_all = np.array(cQ_all)
print('eps_f: max |eps - round(eps)| = %.1e, sum^h eps = %+.5f' % (max(abs(e - round(e)) for e in eps_round), sum(f['om']*f['eps'] for f in forms)))
print('FORMS: M11 = sum omega |Z|^2 = %.6f ; M2 = sum omega Z^2 = %.6f%+.6fi ; L = %.5f' % (M11.real, M2.real, M2.imag, L))
for q in qs: print('FORMS: X_%d = %+.6f%+.6fi' % (q, Xq[q].real, Xq[q].imag))
for Qc in [1, 2, 3, 5, 6, 10]:
    XQ = sum(mobius(q)/q*Xq[q] for q in qs if q <= Qc)
    C2 = float(np.sum(om*cQ_all[:, Qc-1]**2))
    V = sum(f['om']*abs(f['Z'] + f['eps']*cQ_all[i, Qc-1]*L)**2 for i, f in enumerate(forms))
    print('FORMS: Q=%2d  X_Q = %+.6f%+.6fi  C2 = %.6f  V = %.6f  (M11 + 2L Re X_Q + L^2 C2 = %.6f)' % (Qc, XQ.real, XQ.imag, C2, V, M11.real + 2*L*XQ.real + L*L*C2))
Vc = sum(f['om']*abs(f['Z'] + f['eps']*f['cf']*L)**2 for f in forms)
print('FORMS: with the full c_f = B_f/L(1,Sym^2 f): sum^h c_f = %.4f, sum^h c_f^2 = %.4f, V(c_f) = %.6f, D-like ref (1/2pi) int H^2 = %.5f' % (Sc, Sc2, Vc, win.intH2/TWO_PI))
np.savez('forms_%d_%d_%g.npz' % (N, k, T), M11=M11, M2=M2, qs=np.array(qs), Xq=np.array([Xq[q] for q in qs]), Vc=Vc, Sc2=Sc2, L=L)
print('total %.0fs' % (time.time()-t0))
