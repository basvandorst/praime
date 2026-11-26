import mpmath as mp, time
mp.mp.dps = 15
for h in (10, 100, 1000, 5000):
    t = time.time(); n = 200
    for k in range(n): mp.zeta(mp.mpc(0.6, h + 0.01 * k))
    print(h, '%.3f ms' % (1e3 * (time.time() - t) / n))
