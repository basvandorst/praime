import numpy as np
src = open('verify.py').read()
exec(src[:src.index('def best_exact')].replace("assert all(sympy.isprime(v) for v in real)", "").replace("ll_small = ", "ll_small = []#").replace("print('1.", "0 and print('1."))
from grammar import build, apply
P = {nm: fn for nm, fn, _ in build()}
s = rows[188]; X = np.array(s[:-1], np.int64)[None, :]; Y = np.array(s[1:], np.int64)[None, :]; N = np.arange(1, 52)[None, :]
for nm, f in G0.items():
    a = sum(1 for i in range(51) if f(s[i]) == s[i + 1])
    b = int((apply(P[nm], X, N) == Y).sum())
    if a != b:
        F = apply(P[nm], X, N)[0]
        bad = [(s[i], f(s[i]), int(F[i]), s[i + 1]) for i in range(51) if f(s[i]) != F[i]]
        print(nm, 'plain', a, 'numpy', b, bad[:5])
