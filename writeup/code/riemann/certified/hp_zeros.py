import mpmath as mp, numpy as np, time
mp.mp.dps = 30
t0 = time.time(); Z = [mp.zetazero(n).imag for n in range(1, 401)]
np.save('hp_zeros_400.npy', np.array([float(z) for z in Z]))
with open('hp_zeros_400.txt', 'w') as f:
    for z in Z: f.write(mp.nstr(z, 28) + '\n')
G = np.fromfile('z2m.f64')[:400]
d = np.abs(G - np.array([float(z) for z in Z]))
print(f'400 zeros to 28 digits in {time.time()-t0:.0f} s; last {mp.nstr(Z[-1], 20)}; page zeros vs mpmath: max {d.max():.2e} at #{d.argmax()+1}, rms {np.sqrt((d**2).mean()):.2e}')
