# sup over L in [0,1] of the worst exponent A for kappa = 0.24 (a stand-in for Michel-Venkatesh 1/4 - delta),
# with Stucky's estimate set (Coleman LS, Ricci C as stated, Coleman/Huxley large values), and check A(kappa) = 2 + 8 kappa.
import numpy as np
ns = {}
src = open("stucky_exponent.py").read()
exec(src[:src.index("import sys\ncases")], ns)
worstA = ns["worstA"]
res = []
for L in np.linspace(0.0, 1.0, 21):
    A, b, s = worstA(float(L), 0.24, True, True)
    res.append((A, L, b, s))
    print(f"kappa=0.24 L={L:.2f}: A={A:.4f} (beta={b:.3f}, sigma={s:.4f})", flush=True)
A, L, b, s = max(res)
print(f"sup over L: A={A:.4f} at L={L:.2f}; 2+8kappa = {2+8*0.24:.4f}; theta = 1/A = {1/A:.4f} > 1/4: {1/A > 0.25}")
