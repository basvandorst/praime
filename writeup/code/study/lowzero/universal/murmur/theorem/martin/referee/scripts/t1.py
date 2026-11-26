import time, cypari2
pari = cypari2.Pari(); pari.allocatemem(300*10**6)
tf = pari('mftraceform([3*3989,2],0)')
for n in [3989*1201, 3*3989*1201, 3989*9001, 3*3989*9001, 1201, 3*1201, 9001]:
    t=time.time(); c = pari.mfcoef(tf, n); print(n, c, round(time.time()-t,3))
