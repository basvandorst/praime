"""pi(n) - li(n) up to 10^7: sign, extremes, and the log-average of (pi - li)/(-li(sqrt n)/2) (Riemann: pi ~ li(x) - li(sqrt x)/2 - ...)."""
import numpy as np, math
from scipy.special import expi
for N in (500000, 10**7):
    s = np.ones(N + 1, bool); s[:2] = False
    for i in range(2, int(N ** .5) + 1):
        if s[i]: s[i * i::i] = False
    pi = np.cumsum(s)
    n = np.arange(2, N + 1); li = expi(np.log(n))
    d = pi[2:] - li
    g = np.unique(np.exp(np.linspace(math.log(1000), math.log(N), 100000)).astype(int))
    r = (pi[g] - expi(np.log(g))) / (-0.5 * expi(0.5 * np.log(g)))
    print(f'N={N:,}: max of pi - li {d.max():.2f} at {int(n[d.argmax()])}, min {d.min():.2f}; positive anywhere: {bool((d > 0).any())}; '
          f'pi(N) - li(N) = {d[-1]:.2f}, -li(sqrt N)/2 = {-0.5 * expi(0.5 * math.log(N)):.2f}; log-average ratio over [1e3, N]: {r.mean():.3f}')
mu = [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1, 0, -1, 0, -1, 0, 1, 1, -1, 0]
def R_minus_li(x):
    return sum(mu[k] / k * expi(math.log(x) / k) for k in range(2, 25) if mu[k] and x ** (1 / k) > 1.0001)
for N in (500000, 10**7):
    s = np.ones(N + 1, bool); s[:2] = False
    for i in range(2, int(N ** .5) + 1):
        if s[i]: s[i * i::i] = False
    pi = np.cumsum(s)
    g = np.unique(np.exp(np.linspace(math.log(1000), math.log(N), 20000)).astype(int))
    r = np.array([(pi[x] - expi(math.log(x))) / R_minus_li(x) for x in g])
    print(f'N={N:,}: log-average of (pi - li)/(R - li) over [1e3, N]: {r.mean():.3f}; R(N) - li(N) = {R_minus_li(N):.2f}, pi(N) - li(N) = {pi[N] - expi(math.log(N)):.2f}')
