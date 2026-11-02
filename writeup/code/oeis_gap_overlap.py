# Which OEIS sequences share long runs of consecutive differences with the prime gaps? (FINDINGS.txt)
# Scans references/oeis/stripped.gz (first terms of every sequence): a sequence whose differences contain any 8-gram
# of the prime-gap sequence is extended to the longest common run. Run from the repo root (about 10 seconds):
#   .venv/bin/python writeup/code/oeis_gap_overlap.py
import gzip, os
from sympy import primerange
here = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'references', 'oeis')
P = list(primerange(2, 300000)); gaps = [P[i + 1] - P[i] for i in range(len(P) - 1)]
K = 8
grams = {}
for i in range(len(gaps) - K + 1): grams.setdefault(tuple(gaps[i:i + K]), i)
names = {}
with gzip.open(os.path.join(here, 'names.gz'), 'rt') as f:
    for line in f:
        if line.startswith('A'): a, _, nm = line.partition(' '); names[a] = nm.strip()
hits = []
with gzip.open(os.path.join(here, 'stripped.gz'), 'rt') as f:
    for line in f:
        if not line.startswith('A'): continue
        a, _, terms = line.partition(' ')
        t = [int(x) for x in terms.strip().strip(',').split(',') if x and len(x) < 12]
        if len(t) < K + 2: continue
        d = [t[i + 1] - t[i] for i in range(len(t) - 1)]
        if min(d) <= 0 or max(d) > 1500: continue
        best = 0
        for s in range(len(d) - K + 1):
            pos = grams.get(tuple(d[s:s + K]))
            if pos is None: continue
            L = K
            while s + L < len(d) and pos + L < len(gaps) and d[s + L] == gaps[pos + L]: L += 1
            best = max(best, L)
        if best: hits.append((best, len(d), a, names.get(a, '')))
hits.sort(reverse=True)
non = [h for h in hits if 'prime' not in h[3].lower()]
print('%d sequences share >= %d consecutive differences with the prime gaps; %d without "prime" in the name:' % (len(hits), K, len(non)))
for best, n, a, nm in non[:30]: print('  run %2d of %2d  %s %s' % (best, n, a, nm[:100]))
