import numpy as np, mpmath as mp
z = list(np.load('hp_zeros_400.npy'))
txt = [l.split() for l in open('hpA.txt')] + [l.split() for l in open('hpB.txt')]
assert [int(t[0]) for t in txt] == list(range(401, 2420))
strs = [l.strip() for l in open('hp_zeros_400.txt')] + [t[1] for t in txt]
Z = np.array([float(s) for s in strs])
assert np.all(np.diff(Z) > 0)
np.save('hp2419.npy', Z)
open('hp2419.txt', 'w').write('\n'.join(strs) + '\n')
G = np.fromfile('z2m.f64')[:2419]
d = np.abs(G - Z)
print(f'{len(Z)} zeros, last {strs[-1]}; page zeros vs mpmath: max {d.max():.2e} at #{d.argmax()+1}, rms {np.sqrt((d**2).mean()):.2e}, above 1e-6: {int((d>1e-6).sum())}')
# sign-change check of every zero with Z(t) (mpmath siegelz, 40 digits): Z(g - 1e-20) and Z(g + 1e-20) differ in sign
mp.mp.dps = 40
bad = 0
for i, s in enumerate(strs):
    g = mp.mpf(s); a = mp.siegelz(g - mp.mpf('1e-20')); b = mp.siegelz(g + mp.mpf('1e-20'))
    if not (a * b < 0): bad += 1
print(f'sign change of Z within 1e-20 of every listed zero: {len(strs) - bad} of {len(strs)}')
