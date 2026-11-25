# Round 157: the test set of the conductor spectrometer. Seeded draw, no a(p) beyond p <= 100 is computed here.
# Box: a1, a3 in {0,1}, a2 in {-1,0,1}, |a4| <= 20, |a6| <= 30. One curve per (conductor, a_p for p <= 100) (an isogeny
# class has one list a(p)); classes of the ten discovery curves excluded; conductor <= 3000 (the search runs over N' <= 6000,
# so that 2N is among the rivals).  Strata: 20 curves with N <= 500, 20 with 500 < N <= 3000; numpy default_rng(SEED).
import json, numpy as np, cypari2
pari = cypari2.Pari(); pari.allocatemem(200_000_000)
SEED = 20261126 + 157
DISC = {'11a1': [0, -1, 1, -10, -20], '15a1': [1, 1, 1, -10, -10], '27a1': [0, 0, 1, 0, -7], '37a1': [0, 0, 1, -1, 0],
        '43a1': [0, 1, 1, 0, 0], '89a1': [1, 1, 1, -1, 0], '389a1': [0, 1, 1, -2, 0], 'c_A': [1, 0, 1, -3, 1],
        'c_B': [0, 0, 1, -3, 1], '5077a1': [0, 0, 1, -7, 6]}
PR = [int(p) for p in pari.primes(25)]
def key(ai):
    E = pari.ellinit(ai)
    if len(E) == 0: return None
    try: N = int(pari.ellglobalred(E)[0])
    except Exception: return None
    return (N,) + tuple(int(pari.ellap(E, p)) for p in PR)
used = {key(v) for v in DISC.values()}
classes = {}
for a1 in (0, 1):
    for a3 in (0, 1):
        for a2 in (-1, 0, 1):
            for a4 in range(-20, 21):
                for a6 in range(-30, 31):
                    ai = [a1, a2, a3, a4, a6]
                    try: k = key(ai)
                    except Exception: k = None
                    if k is None or k[0] > 3000 or k in used: continue
                    if k not in classes: classes[k] = ai          # first model in the loop order
ks = sorted(classes)
lo = [k for k in ks if k[0] <= 500]; hi = [k for k in ks if k[0] > 500]
rng = np.random.default_rng(SEED)
pick = [lo[i] for i in rng.permutation(len(lo))[:20]] + [hi[i] for i in rng.permutation(len(hi))[:20]]
out = {('T%02d_N%d' % (j, k[0])): classes[k] for j, k in enumerate(pick)}
json.dump(out, open('testset157.json', 'w'), indent=0)
print('classes with N <= 3000 in the box:', len(ks), ' N <= 500:', len(lo), ' 500 < N <= 3000:', len(hi))
print('drawn:', ' '.join(sorted(out, key=lambda s: int(s.split('_N')[1]))))
