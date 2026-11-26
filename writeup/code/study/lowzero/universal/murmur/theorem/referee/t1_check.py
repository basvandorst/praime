# Referee check of T1 (trace formula) with root numbers taken from the L-function (lfunmf + lfuncheckfeq),
# not from a_f(N). Levels chosen outside the note's list.
import cypari2
pari = cypari2.Pari()
pari.allocatemem(400 * 10**6)
pari.set_real_precision(38)

pari("Ucheb(n, x) = polchebyshev(n, 2, x);")
pari("""rhs(N, k, p) = {
  my(s = qfbhclassno(4*p*N)/2, r = 1);
  while (r^2*N < 4*p,
    my(U = Ucheb(k-2, 'X), V = 0);
    for (j = 0, (k-2)/2, V += polcoef(U, 2*j, 'X) * (r^2*N/4)^j * p^((k-2)/2 - j));
    s += (-1)^(k/2-1) * V / p^((k-2)/2) * qfbhclassno(4*p*N - r^2*N^2);
    r++);
  if (k == 2, s -= p + 1);
  s
};""")
pari("""lhs_data(N, k, B) = {
  my(mf = mfinit([N, k], 0), Lb = mfeigenbasis(mf), res = vector(B), feqworst = -1000, nfo = 0, epslist = List());
  for (i = 1, #Lb,
    my(LL = lfunmf(mf, Lb[i]), dg = poldegree(mffields(mf)[i]));
    if (dg == 1, LL = [LL]);
    for (j = 1, #LL,
      my(Lj = LL[j], ch = lfuncheckfeq(Lj), w = round(real(lfunrootres(Lj)[3])), an = lfunan(Lj, B));
      feqworst = max(feqworst, ch);
      listput(epslist, w);
      nfo++;
      for (n = 1, B, res[n] += w * an[n])));
  [res, feqworst, nfo, Vec(epslist)]
};""")

cases = [(167, 2, 700), (173, 2, 700), (197, 2, 800), (29, 4, 150), (31, 4, 150), (17, 6, 100),
         (19, 8, 100), (13, 12, 70), (11, 14, 60), (17, 12, 70)]
worst_all = 0.0
for N, k, B in cases:
    out = pari(f"lhs_data({N},{k},{B})")
    res, feq, nf, eps = out[0], out[1], out[2], out[3]
    dim_cusp = int(pari(f"mfdim([{N},{k}],1)"))
    worst = 0.0
    npairs = 0
    for p in pari(f"primes([2,{B}])"):
        p = int(p)
        if p == N:
            continue
        lhs = complex(res[p - 1]) / p ** ((k - 2) / 2)
        r = float(pari(f"rhs({N},{k},{p})"))
        worst = max(worst, abs(lhs - r))
        npairs += 1
    worst_all = max(worst_all, worst)
    print(f"N={N} k={k}: newforms {int(nf)} (dim S_k(G0(N)) = {dim_cusp}), eps sum {sum(int(e) for e in eps)}, "
          f"worst lfuncheckfeq {int(feq)}, primes {npairs}, max |LHS-RHS| = {worst:.2e}")
print("overall max deviation", worst_all)
