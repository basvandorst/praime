import time, cypari2
pari = cypari2.Pari(); pari.allocatemem(300*10**6)
t=time.time()
print(pari('qfbhclassno(400000003*4-0)'), time.time()-t)
t=time.time()
tf = pari('mftraceform([3*3989,2],0)')
print(pari.mfcoef(tf, 3*3989), pari('mfdim([3*3989,2],0)'), time.time()-t)
t=time.time()
mf = pari('mfinit([3*13,2],0)'); B = pari.mfeigenbasis(mf); print(len(B), pari.mfatkineigenvalues(mf,3), time.time()-t)
