"""Compare gen_data output with exact PARI eigenform coefficients: sum_f a_f(p)^2 = trace over the Hecke fields of a(p)^2."""
import sys, numpy as np, cypari2
pari = cypari2.Pari(); pari.allocatemem(900 * 10 ** 6)
pari('''S2(N,k,B) = my(mf=mfinit([N,k],0), L=mfeigenbasis(mf), P=primes([2,B]), S=vector(#P), K=mffields(mf));
  for(j=1,#L, my(v=mfcoefs(L[j],B), pol=K[j]);
     for(i=1,#P, my(a=v[P[i]+1]); S[i] += if(poldegree(pol)>1, trace(Mod(lift(a)^2,pol)), a^2)));
  [P, S]''')
outdir = sys.argv[1]
for k, N in ((2, 101), (4, 37), (2, 53)):
    d = np.load("%s/data_k%d_N%d.npz" % (outdir, k, N))
    B = int(d["p"].max())
    r = pari('S2(%d,%d,%d)' % (N, k, B))
    ref = {int(p): int(s) for p, s in zip(r[0], r[1])}
    rel = [abs(s * p ** (k - 1) - ref[int(p)]) / ref[int(p)] for p, s in zip(d["p"], d["s"])]
    print("k=%d N=%d primes=%d max rel err %.2e" % (k, N, len(rel), max(rel)), flush=True)
