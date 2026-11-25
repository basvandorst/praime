"""test_set.py: the seeded test set.  Candidates: all fundamental discriminants D with 1500 < |D| <= 2400 (both signs),
distinct conductors (walk keeps the first D of each |D|).  This range is disjoint from round 149's reserved 56 and
development 8 (all |D| <= 1300) and from the 24 development characters of this proposal (1300 < |D| <= 1500).
Selection: numpy default_rng(20261125 + 1500) permutation of the sorted candidate list; four passes over it:
  pass 1 keep D if 7 | D and fewer than 56 kept with 7 | D;  pass 2 the same for 5 (target 64);
  pass 3 for 3 (target 80);  pass 4 for 2 (target 96);  then fill with D coprime to 210 up to 240.
(Group counts include D kept in earlier passes.)  Writes test_set.json."""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pick_dev import fundamental


def main():
    cand = sorted(D for a in range(1501, 2401) for D in (a, -a) if fundamental(D))
    rng = np.random.default_rng(20261125 + 1500)
    perm = [cand[i] for i in rng.permutation(len(cand))]
    kept, absk = [], set()
    count = lambda p: sum(1 for D in kept if D % p == 0)
    for p, tgt in ((7, 56), (5, 64), (3, 80), (2, 96)):
        for D in perm:
            if count(p) >= tgt:
                break
            if abs(D) in absk or D % p:
                continue
            kept.append(D); absk.add(abs(D))
    for D in perm:
        if len(kept) >= 240:
            break
        if abs(D) in absk or any(D % p == 0 for p in (2, 3, 5, 7)):
            continue
        kept.append(D); absk.add(abs(D))
    kept.sort()
    groups = {p: count(p) for p in (2, 3, 5, 7, 11, 13)}
    print('candidates %d (distinct |D| %d); kept %d; groups %s; none of 2,3,5,7: %d; D<0: %d' % (
        len(cand), len(set(abs(D) for D in cand)), len(kept), groups,
        sum(1 for D in kept if all(D % p for p in (2, 3, 5, 7))), sum(1 for D in kept if D < 0)))
    json.dump({'test': kept, 'groups': groups,
               'rule': 'test_set.py: fundamental D, 1500 < |D| <= 2400, distinct |D|; default_rng(20261125 + 1500) '
                       'permutation of the sorted list; passes for 7 (56), 5 (64), 3 (80), 2 (96); fill with D '
                       'coprime to 210 to 240'}, open(sys.argv[1], 'w'), indent=1)


if __name__ == '__main__':
    main()
