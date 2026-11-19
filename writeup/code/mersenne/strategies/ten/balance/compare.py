import json, sys
A = json.load(open(sys.argv[1])); B = json.load(open(sys.argv[2]))
bad = []; n = 0
for key in sorted(A):
    a, b = A[key], B[key]; n += 1
    exact = key.endswith('|one') or key.endswith('|val')
    for f in ('k', 'eq', 'aL', 'aR', 'ch', 'bL', 'bR', 'fL', 'fR', 'eqsteps'):
        if a[f] != b[f]: bad.append((key, f, a[f], b[f]))
    for f in ('L', 'R'):
        if exact and a[f] != b[f]: bad.append((key, f, a[f], b[f]))
        if not exact and abs(a[f] - b[f]) > 1e-9 * max(1, abs(a[f])): bad.append((key, f, a[f], b[f]))
print('%d runs compared (count, totals, equal steps, time ahead, lead changes, longest leads, first leads): %d differences' % (n, len(bad)))
for x in bad[:20]: print('  ', x)
