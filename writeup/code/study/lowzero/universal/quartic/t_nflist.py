import cypari2, time
pari = cypari2.Pari(); pari.allocatemem(300_000_000)
print(pari('default(datadir)'), pari.version())
for s in (0,1,2):
    t=time.time(); L = pari('nflist("S4", [1, 3000], %d)'%s); print(s, len(L), '%.1fs'%(time.time()-t), L[0] if len(L) else None, flush=True)
