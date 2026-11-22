import math, sys
import region2 as R
for g0, beta in [(14.135, 0.99), (14.135, 0.9), (14.135, 0.75), (1e3, 0.9), (1e3, 0.6), (1e6, 0.9)]:
    yn = R.yneeded(beta, g0, rig=True)
    print(f"rigorous (Turan) one-sided: t={g0}, beta={beta}: needs Y1 = 10^{yn if yn is None else f'{yn:.1f}'}", flush=True)
