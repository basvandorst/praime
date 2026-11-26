import cypari2
pari = cypari2.Pari(); pari.allocatemem(300_000_000)
for r2 in (0, 2):
    full = sorted(int(pari.nfdisc(f)) for f in pari('nflist("S4", [1, 20000], %d)' % r2))
    ch = []
    for c0 in range(1, 20000, 2500):
        L = [int(pari.nfdisc(f)) for f in pari('nflist("S4", [%d, %d], %d)' % (c0, min(20000, c0 + 2499), r2))]
        print(r2, c0, len(L), 'outside:', sum(1 for d in L if not c0 <= abs(d) <= c0 + 2499), 'dups:', len(L) - len(set(L)))
        ch += [d for d in L if c0 <= abs(d) <= c0 + 2499]
    print('r2', r2, 'full', len(full), len(set(full)), 'chunks filtered', len(ch), len(set(ch)), 'same set:', set(full) == set(ch))
    print(' full within [1,20000]:', all(1 <= abs(d) <= 20000 for d in full), ' discs with >1 field:', len(full) - len(set(full)))
