# Check the model positions against the page's extracted canvas positions.
import numpy as np, common as C
for name in C.LAYOUTS:
    L = C.layout(name); fx, fy = C.load_xy(name)
    d = np.hypot(L['x'][1:]-fx[1:], L['y'][1:]-fy[1:])
    print(f"{name:11s} max|model-page| = {np.nanmax(d):.4f} px   median {np.nanmedian(d):.2e}")
