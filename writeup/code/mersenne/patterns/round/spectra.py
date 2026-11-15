"""Angular power spectra of primes vs all integers in annuli (blocks of consecutive n; every
round layout has radius monotone in n).  For layouts whose angle is linear in n
(theta = 2*pi*alpha*n), harmonic k of the primes is the exponential sum sum_p e(k*alpha*p); by
the prime number theorem in progressions (Siegel-Walfisz) its normalised amplitude is
|mu(q)|/phi(q) * |F(beta)| where k*alpha = a/q + beta and F is the block's Dirichlet-kernel
factor.  We compare observed amplitude A_obs = sqrt(sum(|S|^2-Np)/sum(Np^2)) with that."""
import numpy as np, common as C, sys, json
from sympy import totient, mobius
isp = C.sieve()
M = 1 << 16

def harm(theta, K):
    """sum_n e^{i k theta_n}, k = 1..K, via linear-interpolated histogram + FFT (deconvolved)."""
    u = np.mod(theta, 2*np.pi)/(2*np.pi)*M
    i0 = np.floor(u).astype(np.int64); f = u - i0
    h = np.bincount(i0 % M, 1-f, M) + np.bincount((i0+1) % M, f, M)
    F = np.conj(np.fft.fft(h))[1:K+1]  # sum e^{+i k theta}
    k = np.arange(1, K+1); s = np.sinc(k/M)**2
    return F/s

def harmw(theta, wts, K):
    u = np.mod(theta, 2*np.pi)/(2*np.pi)*M
    i0 = np.floor(u).astype(np.int64); f = u - i0
    h = np.bincount(i0 % M, wts*(1-f), M) + np.bincount((i0+1) % M, wts*f, M)
    F = np.conj(np.fft.fft(h))[1:K+1]
    k = np.arange(1, K+1); s = np.sinc(k/M)**2
    return F/s

def blockspec(theta, blocks, K, mask):
    """Primes: A^2 = sum_b(|S_b|^2 - noise_b)/sum_b Np_b^2.  Model (PNT, no residue structure):
    each block's integers weighted by 1/ln n, normalised, combined with the same weights Np_b^2."""
    num = np.zeros(K); den = 0.0; cnt = 0; modnum = np.zeros(K)
    for a, b in blocks:
        sel = np.arange(a, b)
        selp = sel[mask[a:b]]
        if len(selp) < 3: continue
        S = harm(theta[selp], K)
        P2 = np.abs(S)**2
        noise = np.median(P2)/np.log(2) if len(P2) >= 20 else len(selp)
        num += P2 - noise; den += len(selp)**2; cnt += len(selp)
        wl = 1/np.log(np.maximum(sel, 2).astype(float))
        Sm = harmw(theta[sel], wl, K)/wl.sum()
        modnum += len(selp)**2*np.abs(Sm)**2
    return np.sqrt(np.maximum(num, 0)/den), cnt, np.sqrt(modnum/den)

def best_rational(x, Delta, Q=400):
    best = (0, 1, 0.0, 0.0)
    for q in range(1, Q+1):
        mu = int(mobius(q))
        if mu == 0: continue
        a = round(x*q); beta = x - a/q
        z = np.pi*beta*Delta
        Fb = 1.0 if abs(z) < 1e-12 else abs(np.sin(z)/(Delta*np.sin(np.pi*beta)))
        amp = Fb/int(totient(q))
        if amp > best[3]: best = (a % q, q, beta, amp)
    return best

ints = np.ones(C.N+1, bool); ints[0] = False
res = {}
EXTRA = {'polar': [4, 22, 44, 333], 'sunflower': [17, 21, 34, 55, 72, 89, 144, 305, 610], 'wheel': [], 'sacks': [1, 2], 'theodorus': [1, 2], 'doubling': [1, 2]}
def run(name, alpha_fn, groups, K, top=8, Delta_fn=None, exact_q=None):
    L = C.layout(name); th = L['theta']
    out = []
    for gname, blocks in groups:
        Ap, npr, Ai = blockspec(th, blocks, K, isp)
        Delta = int(np.median([b-a for a, b in blocks]))
        idx = list(np.argsort(-Ap)[:top]) + [k-1 for k in EXTRA.get(name, []) if k <= K]
        idx = sorted(set(int(i) for i in idx))
        print(f"\n{name} {gname}: {len(blocks)} blocks of ~{Delta} numbers, {npr} primes")
        print("   k   A_primes  A_model(ints/ln n)   k*alpha~a/q   mu/phi*F (pred)")
        rows = []
        for i in idx:
            k = i+1
            if exact_q is not None:
                from math import gcd
                q = exact_q//gcd(k, exact_q); amp = abs(int(mobius(q)))/int(totient(q))
                pred = f"exact q={q}  {amp:.3f}"
            elif alpha_fn is not None:
                al = alpha_fn(blocks)
                a, q, beta, amp = best_rational((k*al) % 1, Delta)
                pred = f"{a}/{q} (beta={beta:+.1e})  {amp:.3f}"
            else:
                pred = "-"; q = None; amp = None
            print(f"{k:5d}  {Ap[i]:.3f}   {Ai[i]:.3f}   {pred}")
            rows.append(dict(k=int(k), A_primes=float(Ap[i]), A_ints=float(Ai[i]), q=None if q is None else int(q), pred=None if amp is None else float(amp)))
        out.append(dict(group=gname, nblocks=len(blocks), Delta=Delta, peaks=rows,
                        noise=float(np.median(Ap))))
        print(f"   median A_primes over k = {np.median(Ap):.4f} (noise level)")
    res[name] = out

def nblocks(n0, n1, D):
    return [(a, min(a+D, n1)) for a in range(n0, n1, D) if min(a+D, n1)-a >= D//2]

which = sys.argv[1:] or ['polar','sunflower','wheel','sacks','theodorus','doubling','binring','decrings','ternrings']
for name in which:
    if name == 'polar':
        al = lambda b: 1/(2*np.pi)
        run(name, al, [('inner n<50000, D=300', nblocks(1000, 50000, 300)),
                       ('outer n>=50000, D=300', nblocks(50000, C.N+1, 300)),
                       ('outer n>=50000, D=25000', nblocks(50000, C.N+1, 25000))], 1500)
    if name == 'sunflower':
        al = lambda b: (3-np.sqrt(5))/2
        run(name, al, [('n<20000, D=400', nblocks(1000, 20000, 400)),
                       ('n>=100000, D=400', nblocks(100000, C.N+1, 400)),
                       ('n>=100000, D=3000', nblocks(100000, C.N+1, 3000))], 1500)
    if name == 'wheel':
        al = lambda b: 1/708
        run(name, al, [('all rings, 10 rings per block', nblocks(1, C.N+1, 7080))], 708, top=12, exact_q=708)
    if name == 'sacks':
        bl = [(m*m, (m+1)**2) for m in range(30, 707)]
        run(name, None, [('one turn per block, m=30..706', bl)], 60)
    if name == 'theodorus':
        run(name, None, [('D=2000', nblocks(1000, C.N+1, 2000))], 200)
    if name == 'doubling':
        bl = [(2**b, min(2**(b+1), C.N+1)) for b in range(6, 19)]
        run(name, None, [('one doubling per block, b=6..18', bl)], 40)
    if name in ('binring', 'decrings', 'ternrings'):
        B = {'binring': 2, 'decrings': 10, 'ternrings': 3}[name]
        L = C.layout(name)
        d = 1
        rings = []
        while B**(d-1) <= C.N:
            a = B**(d-1); b = min(B**d, C.N+1)
            if B == 2: cnt = a
            else: cnt = (B-1)*a
            rings.append((a, b, cnt)); d += 1
        for a, b, cnt in rings:
            if cnt < 200 or cnt > 400000: continue
            run(name, (lambda bb, cnt=cnt: 1/cnt), [(f'ring {a}..{b-1}', [(a, b)])], min(2048, cnt//2+1), top=5)
json.dump(res, open('out/spectra_%s.json' % '_'.join(which), 'w'), indent=1)
