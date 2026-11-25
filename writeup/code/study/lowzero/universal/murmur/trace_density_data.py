"""Natural-weight murmuration data from the (validated, trace_validate.py) class-number formula, k = 2:
    T(p,N) = sum_{f in H_2^new(N)} eps_f a_f(p) = H(4pN)/2 + sum_{r>=1, r^2 N < 4p} H(4pN - r^2 N^2) - (p+1)
Usage: trace_density_data.py FAMILY N0 N1 YMAX OUT.npz     FAMILY in {prime, squarefree}; env BUDGET = wall seconds (default 450),
  results saved after every level (sums are additive over chunks)
  prime      : prime levels N in [N0, N1]
  squarefree : all squarefree N in [N0, N1] (Zubrilina's family)
Primes p <= YMAX*N, p not dividing N. Saves per-bin sums (bins of width 0.25 in y = p/N) of T, dim, dim*M_zub, dim*M_prime,
and the transform sums  sum_N dim_N * [-(1/2pi) sum_p (ln p/p) (T/dim_N) G(p/N)]  for two Gaussian windows.
"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import cypari2
from densities import M_prime, M_zub, G_kernel

pari = cypari2.Pari()
pari.allocatemem(500 * 10 ** 6)
pari('T(p,N) = my(s = qfbhclassno(4*p*N)/2, r = 1); while(r^2*N < 4*p, s += qfbhclassno(4*p*N - r^2*N^2); r++); s - (p+1)')
pari('V(N, Y) = my(P = select(q -> N % q, primes([2, floor(Y*N)]))); [P, apply(q -> T(q, N), P)]')

fam, N0, N1, YMAX, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
WINDOWS = ((10.0, 2.0), (14.0, 2.0))
if fam == "prime":
    levels = [N for N in range(N0, N1 + 1) if pari('isprime(%d)' % N)]
else:
    levels = [N for N in range(N0, N1 + 1) if pari('issquarefree(%d)' % N)]

ygrid = np.linspace(1e-4, YMAX + 0.01, 100001)
Ggrid = {w: G_kernel(ygrid, *w, nt=1201) for w in WINDOWS}


def G_interp(y, w):
    g = Ggrid[w]
    return np.interp(y, ygrid, g.real) + 1j * np.interp(y, ygrid, g.imag)


edges = np.arange(0.0, YMAX + 1e-9, 0.25)
nb = len(edges) - 1
S = {k: np.zeros(nb) for k in ("T", "dim", "zub", "prime", "T2", "cnt")}
trans = {w: 0j for w in WINDOWS}
wsum = 0.0
t0 = time.time()
BUDGET = float(os.environ.get("BUDGET", "450"))   # wall-clock seconds; stop cleanly and save what is done
done = 0


def save():
    np.savez(out, edges=edges, nlev=done, wsum=wsum, **S,
             trans=np.array([trans[w] for w in WINDOWS]), windows=np.array(WINDOWS), last=last)


last = 0
for N in levels:
    if time.time() - t0 > BUDGET:
        break
    dim = int(pari('mfdim([%d, 2], 0)' % N))
    if dim == 0:
        continue
    P, Tv = pari('V(%d, %s)' % (N, YMAX))
    P = np.array([int(q) for q in P], dtype=float)
    Tv = np.array([float(v) for v in Tv])
    y = P / N
    b = np.digitize(y, edges) - 1
    ok = (b >= 0) & (b < nb)
    np.add.at(S["T"], b[ok], Tv[ok])
    np.add.at(S["T2"], b[ok], (Tv[ok] / dim) ** 2 * dim)
    np.add.at(S["cnt"], b[ok], 1.0)
    np.add.at(S["dim"], b[ok], dim)
    np.add.at(S["zub"], b[ok], dim * M_zub(y[ok]))
    np.add.at(S["prime"], b[ok], dim * M_prime(y[ok]))
    for w in WINDOWS:
        G = G_interp(y, w)
        trans[w] += dim * (-(1 / (2 * np.pi)) * np.sum(np.log(P) / P * (Tv / dim) * G))
    wsum += dim
    done += 1
    last = N
    save()
print(f"{fam} [{N0},{N1}]: {done} of {len(levels)} levels done (last N = {last}), total dim {wsum:.0f}, {time.time() - t0:.0f}s", flush=True)
