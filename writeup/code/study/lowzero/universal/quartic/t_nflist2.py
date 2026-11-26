import cypari2, time
pari = cypari2.Pari(); pari.allocatemem(300_000_000)
for (a,b) in ((5000,9999),(10000,19999)):
  for s in (0,1,2):
    t=time.time(); L = pari('nflist("S4", [%d, %d], %d)'%(a,b,s)); print((a,b), s, len(L), '%.1fs'%(time.time()-t), flush=True)
