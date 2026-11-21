"""Run every probe on one window for one generator and append one JSON line per draw.

usage: engine.py ROUND WINDOW GEN DRAWS [SEED0]   (GEN = primes, cramer, sieve13, sieve47, sieve100, class, shuffle,
                                                  or theory: writes the theory values once)
seed of draw j = SEED0 + 1000 * hash-free index of GEN + j, so every draw can be replayed.
"""
import json, sys, time, os
import numpy as np
import common as C
import probes as PR

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data', 'study')


def main():
    rnd, win, gen, draws = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    seed0 = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    L, W = C.WINDOWS[win]
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, 'results_%s.jsonl' % rnd)
    if gen == 'theory':
        rec = {'round': rnd, 'window': win, 'gen': 'theory', 'stats': PR.theory_values(L, W)}
        with open(path, 'a') as f:
            f.write(json.dumps(rec) + '\n')
        return
    if gen == 'primes':
        # whole window and K blocks; each block is compared with the theory of that block
        K = 10
        M = C.primes_mask(L, W)
        with open(path, 'a') as f:
            for b in range(-1, K):
                t = time.time()
                Lb, Wb = (L, W) if b < 0 else (L + b * (W // K), W // K)
                Mb = M[Lb - L:Lb - L + Wb]
                rec = {'round': rnd, 'window': win, 'gen': 'primes', 'block': b, 'L': Lb, 'W': Wb,
                       'members': int(Mb.sum()), 'stats': PR.all_stats(Mb, Lb, Wb),
                       'theory': PR.theory_values(Lb, Wb), 'secs': None}
                rec['secs'] = round(time.time() - t, 1)
                f.write(json.dumps(rec) + '\n')
                print('primes block', b, rec['secs'], flush=True)
        return
    Mp = C.primes_mask(L, W) if gen == 'shuffle' else None
    gi = C.GENS.index(gen)
    for j in range(draws):
        t = time.time()
        seed = seed0 + 1000 * gi + j
        rng = np.random.default_rng(seed)
        M = C.make_mask(gen, L, W, rng, Mp)
        st = PR.all_stats(M, L, W)
        rec = {'round': rnd, 'window': win, 'gen': gen, 'draw': j, 'seed': seed, 'members': int(M.sum()),
               'secs': round(time.time() - t, 1), 'stats': st}
        del M
        with open(path, 'a') as f:
            f.write(json.dumps(rec) + '\n')
        print(gen, j, rec['secs'], flush=True)


if __name__ == '__main__':
    main()
