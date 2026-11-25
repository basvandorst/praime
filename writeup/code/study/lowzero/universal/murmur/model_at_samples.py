"""Transform of the MODEL densities sampled at exactly the data points (same levels, same primes p <= YMAX*N, p not dividing N,
same weights dim_N, ln p/p) as trace_density_data.py: separates the prime-sampling effect from the data-vs-density comparison.
Usage: model_at_samples.py FAMILY N0 N1 YMAX [LASTN]   (LASTN: stop after this level, to match a chunk)"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import sympy as sp
import cypari2
from densities import M_prime, M_zub, G_kernel, H_int

pari = cypari2.Pari()
fam, N0, N1, YMAX = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
LAST = int(sys.argv[5]) if len(sys.argv) > 5 else N1
WINDOWS = ((10.0, 2.0), (14.0, 2.0))
ygrid = np.linspace(1e-4, YMAX + 0.01, 100001)
Gg = {w: G_kernel(ygrid, *w, nt=1201) for w in WINDOWS}
PR = np.array(list(sp.primerange(2, int(YMAX * N1) + 2)), dtype=float)
acc = {(w, m): 0j for w in WINDOWS for m in ("prime", "zub")}
wsum = 0.0
nlev = 0
for N in range(N0, min(N1, LAST) + 1):
    if fam == "prime" and not sp.isprime(N):
        continue
    if fam == "squarefree" and any(e > 1 for e in sp.factorint(N).values()):
        continue
    dim = int(pari('mfdim([%d, 2], 0)' % N))
    if dim == 0:
        continue
    P = PR[(PR <= np.floor(YMAX * N)) & (np.mod(N, PR) != 0)]
    y = P / N
    for w in WINDOWS:
        G = np.interp(y, ygrid, Gg[w].real) + 1j * np.interp(y, ygrid, Gg[w].imag)
        for m, M in (("prime", M_prime), ("zub", M_zub)):
            acc[(w, m)] += dim * (-(1 / (2 * np.pi)) * np.sum(np.log(P) / P * M(y) * G))
    wsum += dim
    nlev += 1
print(f"{fam}: {nlev} levels, total dim {wsum:.0f}")
for w in WINDOWS:
    norm = -H_int(*w) / (2 * np.pi)
    a = acc[(w, "prime")] / wsum / norm
    b = acc[(w, "zub")] / wsum / norm
    print(f"   t0={w[0]}: model sampled at the data points / (-int H/2pi):  M_prime {a.real:+.4f}{a.imag:+.4f}i   M_zub {b.real:+.4f}{b.imag:+.4f}i")
