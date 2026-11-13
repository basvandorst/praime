import mpmath as mp, sys, time
mp.mp.dps = 30
lo, hi, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
t0 = time.time()
with open(out, 'w') as f:
    for n in range(lo, hi):
        f.write(f'{n} {mp.nstr(mp.zetazero(n).imag, 28)}\n'); f.flush()
print(f'zeros {lo}..{hi-1} done in {time.time()-t0:.0f} s')
