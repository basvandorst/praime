# Worst exponent A(kappa) for the twisted family at L near 1 (full norm range, T <= x^eps),
# using the estimates of Stucky's paper; kappa = exponent of the pointwise bound on Re s = 1/2.
import numpy as np, sys
sys.path.insert(0, '.')
import importlib.util
spec = importlib.util.spec_from_file_location("se", "stucky_exponent.py")
src = open("stucky_exponent.py").read()
src = src[:src.index("import sys\ncases")]
ns = {}
exec(src, ns)
worstA = ns["worstA"]
for L in [1.0, 0.98, 0.9]:
    for kappa in [0.25, 0.249, 0.245, 0.24, 0.23, 0.2, 1/6]:
        for C in [True, False]:
            A, b, s = worstA(L, kappa, C, True)
            print(f"L={L:.2f} kappa={kappa:.4f} RicciC={str(C):5s} A={A:.4f} 1/A={1/A:.4f} at beta={b:.3f} sigma={s:.4f}", flush=True)
