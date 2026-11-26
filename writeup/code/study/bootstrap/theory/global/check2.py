import numpy as np
# corrected Lemma H' identity for zeta: sum_{n != 0} c sinc^2(pi c n) = 1 - c for 0 < c <= 1
n = np.arange(1, 2_000_001, dtype=np.float64)
for c in [0.3, 0.77, 1.0]:
    S = 2*np.sum(c*(np.sin(np.pi*c*n)/(np.pi*c*n))**2)
    tail = 2*c/(2*np.pi**2*c**2*n[-1]) if c < 1 else 0.0   # mean of sin^2 is 1/2
    print("c=%.2f  partial sum + tail = %.8f   1-c = %.8f" % (c, S + tail, 1 - c))
