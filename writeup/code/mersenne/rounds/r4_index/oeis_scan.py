"""OEIS (local stripped.gz / names.gz, read-only) for every derived sequence, real and null.

Relatives: sequences whose name mentions Mersenne, perfect numbers, Lucas-Lehmer, '2^p-1 is prime', 2^prime(n)-1,
or the A-numbers A000043, A001348, A016027, A000668, A000396;
they are excluded from the statistics (they restate the definition) and listed apart for the real sequences.
  prefix : largest k such that b_1..b_k occur as consecutive terms of a non-relative OEIS sequence
           (all queries at once: a trie of the query prefixes, walked from every position of every sequence).
           Values below 4 are reported as found (any k); the statistic is k.
  block  : (real only) longest b_i..b_j (any i) consecutive in a non-relative sequence (trie of all suffixes).
  inorder: (real only, increasing sequences) longest common subsequence with an OEIS sequence, and its purity
           prec = inorder / #terms of s between b_1 and the last common term (oeis.py of section 6).
Output: oeis.json {seq: {'prefix': {'side': 'hi', 'real': k, 'null': [...]}}}, oeis_real.txt
"""
import gzip, re, json, bisect, time, sys
import numpy as np

OE = '/home/bas/praime/references/oeis/'
D = np.load('data.npz')
NAMES = [str(x) for x in D['names']]
REL = re.compile(r'mersenne|perfect numbers?|even perfect|lucas.?lehmer|2\^[pnkm] *- *1 is (a )?prime|2\^prime\(n\) *- *1|2\^\(prime|A000043|A001348|A016027|A000668|A000396', re.I)
t0 = time.time()
names = {}
relset = set()
with gzip.open(OE + 'names.gz', 'rt', encoding='utf-8', errors='replace') as f:
    for line in f:
        if line.startswith('A'):
            an = line[:7]; nm = line[8:].strip(); names[an] = nm
            if REL.search(nm): relset.add(an)
print('names', len(names), 'relatives', len(relset), '%.1f s' % (time.time() - t0), flush=True)

# queries
Q = []   # (seq name, row, tuple of terms)
for nm in NAMES:
    M = D[nm]; Ln = D[nm + '_len']
    for r in range(M.shape[0]):
        Q.append((nm, r, tuple(int(x) for x in M[r, :Ln[r]])))
# prefix trie: node id -> dict(value -> child id); depth per node; visited flag
child = [dict()]
depth = [0]
qpath = []
for nm, r, t in Q:
    node = 0; path = [0]
    for v in t:
        nx = child[node].get(v)
        if nx is None:
            nx = len(child); child.append(dict()); depth.append(depth[node] + 1); child[node][v] = nx
        node = nx; path.append(node)
    qpath.append(path)
visited = np.zeros(len(child), bool)
visitedA = {}  # node -> first A-number (for real query nodes)
realnodes = set()
for (nm, r, t), path in zip(Q, qpath):
    if r == 0: realnodes.update(path)
# suffix trie for real queries (block statistic)
schild = [dict()]; sdepth = [0]; sfirst = {}; sowner = {}
for nm in NAMES:
    t = [int(x) for x in D[nm][0, :D[nm + '_len'][0]]]
    for i in range(len(t)):
        node = 0
        for j in range(i, len(t)):
            v = t[j]
            nx = schild[node].get(v)
            if nx is None:
                nx = len(schild); schild.append(dict()); sdepth.append(sdepth[node] + 1); schild[node][v] = nx
            node = nx
            sowner.setdefault(node, set()).add((nm, i))
svis = {}
# in-order for increasing real sequences
INC = [nm for nm in NAMES if all(np.diff(D[nm][0, :D[nm + '_len'][0]]) > 0)]
IDX = {nm: {int(v): i for i, v in enumerate(D[nm][0, :D[nm + '_len'][0]])} for nm in INC}
REALT = {nm: [int(v) for v in D[nm][0, :D[nm + '_len'][0]]] for nm in NAMES}
inorder = {nm: [] for nm in INC}
relhits = {nm: [] for nm in NAMES}  # relatives that contain the real prefix
nseq = 0
with gzip.open(OE + 'stripped.gz', 'rt') as f:
    for line in f:
        if not line.startswith('A'): continue
        an, _, tt = line.partition(' ')
        try:
            s = [int(x) for x in tt.strip().strip(',').split(',') if x]
        except ValueError:
            continue
        if len(s) < 3: continue
        nseq += 1
        isrel = an in relset
        n = len(s)
        for p in range(n):
            node = 0
            for q in range(p, n):
                nx = child[node].get(s[q])
                if nx is None: break
                node = nx
                if isrel:
                    continue
                if not visited[node]:
                    visited[node] = True
                    if node in realnodes: visitedA[node] = an
            # suffix trie (real block)
            node = 0
            for q in range(p, n):
                nx = schild[node].get(s[q])
                if nx is None: break
                node = nx
            if node and sdepth[node] >= 4:
                key = 'rel' if isrel else 'non'
                cur = svis.get((node, key))
                if cur is None: svis[(node, key)] = an
        for nm in INC:
            ix = IDX[nm]
            hit = [ix.get(v, -1) for v in s]
            if sum(h >= 0 for h in hit) < 4: continue
            tails = []
            for h in hit:
                if h >= 0:
                    k = bisect.bisect_left(tails, h)
                    if k == len(tails): tails.append(h)
                    else: tails[k] = h
            last = REALT[nm][tails[-1]]
            nin = sum(1 for v in s if REALT[nm][0] <= v <= last)
            inorder[nm].append((len(tails), len(tails) / max(nin, 1), an, isrel))
print('scanned', nseq, 'sequences %.1f s' % (time.time() - t0), flush=True)

# prefix statistic per query
res = {}
pref = {}
for (nm, r, t), path in zip(Q, qpath):
    k = 0
    for d, node in enumerate(path):
        if d and visited[node]: k = d
    pref.setdefault(nm, {})[r] = k
for nm in NAMES:
    v = pref[nm]
    res[nm] = {'prefix': {'side': 'hi', 'real': v[0], 'null': [v[r] for r in sorted(v) if r > 0]}}
out = []
for nm in NAMES:
    t = REALT[nm]
    path = qpath[[i for i, q in enumerate(Q) if q[0] == nm and q[1] == 0][0]]
    k = res[nm]['prefix']['real']
    an = visitedA.get(path[k]) if k else None
    nul = np.array(res[nm]['prefix']['null'])
    p = (1 + (nul >= k).sum()) / (1 + len(nul))
    out.append(f'== {nm}: real {t[:15]}')
    out.append(f'   prefix (non-relatives): k = {k} (first in {an}: {names.get(an, "")[:90]}); '
               f'null median {np.median(nul):.0f}, 95% {np.quantile(nul, .95):.0f}, max {nul.max()}, p = {p:.3f}')
    # block: deepest suffix node visited by a non-relative / by a relative
    for key in ('non', 'rel'):
        bestb = (0, None, None)
        for (node, kk), a in svis.items():
            if kk != key: continue
            own = [o for o in sowner[node] if o[0] == nm]
            if own and sdepth[node] > bestb[0]:
                bestb = (sdepth[node], min(o[1] for o in own), a)
        if bestb[0]:
            out.append(f'   block ({"relatives" if key == "rel" else "non-relatives"}): {bestb[0]} consecutive terms from b_{bestb[1] + 1} = {t[bestb[1]]} in {bestb[2]}: {names.get(bestb[2], "")[:90]}')
    if nm in inorder:
        L = sorted(inorder[nm], key=lambda x: (x[0], x[1]), reverse=True)
        top = [x for x in L if not x[3]][:3]
        topr = [x for x in L if x[3]][:3]
        Lp = sorted([x for x in L if x[1] >= 0.5 and not x[3]], key=lambda x: (x[0], x[1]), reverse=True)[:3]
        out.append('   in-order (non-relatives): ' + '; '.join(f'{a} {c} terms prec {pp:.2f}' for c, pp, a, _ in top))
        out.append('   in-order, prec >= 0.5 (non-relatives): ' + '; '.join(f'{a} {c} terms prec {pp:.2f} {names.get(a, "")[:60]}' for c, pp, a, _ in Lp))
        out.append('   in-order (relatives): ' + '; '.join(f'{a} {c} terms prec {pp:.2f} {names.get(a, "")[:70]}' for c, pp, a, _ in topr))
json.dump(res, open('oeis.json', 'w'))
open('oeis_real.txt', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
print('time %.1f s' % (time.time() - t0))
