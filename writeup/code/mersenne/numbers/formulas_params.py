import numpy as np, math
src = open('formulas.py').read()
exec(src[src.index('def strips'):src.index('def exact_check')])
D = np.load('seqs.npz'); real = D['real']
b, (u, v) = best_hits(real, 'F1')
c, r = math.exp(u), math.exp(v)
print(b, repr(c), repr(r))
open('floor_cr_params.txt', 'w').write(f'{c!r} {r!r}\n')
