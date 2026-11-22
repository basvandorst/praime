import flint, sys
flint.ctx.prec = 64
N = int(sys.argv[1]); out=[]
step=500
for s in range(1, N+1, step):
    zs = flint.acb.zeta_zeros(s, min(step, N-s+1))
    out += [float(z.imag.mid()) for z in zs]
open('zeros_%d.txt'%N,'w').write('\n'.join('%.15f'%g for g in out))
print(len(out), out[-1])
