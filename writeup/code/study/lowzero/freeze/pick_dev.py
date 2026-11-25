"""pick_dev.py: choose the 24 extra development characters, 1300 < |D| <= 1500, fundamental, distinct |D|,
seeded and balanced in p = 2, 3, 5, 7 (greedy on a seeded permutation: accept D if it brings a prime group
closer to its target, else fill).  Writes dev_extra.json."""
import json, sys
import numpy as np


def squarefree(n):
    n = abs(n)
    d = 2
    while d * d <= n:
        if n % (d * d) == 0:
            return False
        d += 1
    return True


def fundamental(D):
    if D in (0, 1):
        return False
    if D % 4 == 1 and squarefree(D):
        return True
    if D % 4 == 0:
        m = D // 4
        if m % 4 in (2, 3) and squarefree(m):
            return True
    return False


def main():
    cand = [D for a in range(1301, 1501) for D in (a, -a) if fundamental(D)]
    rng = np.random.default_rng(20261125 + 149)
    perm = [cand[i] for i in rng.permutation(len(cand))]
    target = {2: 8, 3: 8, 5: 7, 7: 7}
    have = {p: 0 for p in target}
    kept, absk = [], set()
    # pass 1: greedy for prime groups
    for D in perm:
        if abs(D) in absk:
            continue
        bad = [p for p in target if D % p == 0]
        if any(have[p] < target[p] for p in bad):
            kept.append(D); absk.add(abs(D))
            for p in bad:
                have[p] += 1
        if len(kept) >= 24 or all(have[p] >= target[p] for p in target):
            break
    # pass 2: fill with D having none of 2, 3, 5, 7
    for D in perm:
        if len(kept) >= 24:
            break
        if abs(D) in absk or any(D % p == 0 for p in target):
            continue
        kept.append(D); absk.add(abs(D))
    kept.sort()
    have = {p: sum(1 for D in kept if D % p == 0) for p in target}
    print('candidates', len(cand), 'kept', len(kept), kept)
    print('groups', have)
    json.dump({'dev_extra': kept, 'rule': 'pick_dev.py: fundamental 1300<|D|<=1500, default_rng(20261274) permutation, '
               'greedy to targets 2:8 3:8 5:7 7:7, fill with D coprime to 210'}, open(sys.argv[1], 'w'), indent=1)


if __name__ == '__main__':
    main()
