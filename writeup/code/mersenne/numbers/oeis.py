"""OEIS (local stripped.gz, read-only): sequences that agree longest with A000043.

For each sequence s (terms as stored in stripped.gz, about 260 characters):
  prefix  = largest k with a_1..a_k consecutive somewhere in s (pos = where the run starts in s);
  block   = longest a_i..a_j (any i) that occurs consecutively in s;
  inorder = longest common subsequence of s and a (s may hold other terms between),
            computed as a longest strictly increasing run of a-indices (patience sorting);
  prec    = inorder / #terms of s in [2, last common term]  (how pure the agreement is).
Relatives of A000043 (name mentions Mersenne / perfect / Lucas-Lehmer / 2^p-1 style) are listed apart.
"""
import gzip, re, bisect, sys, os
sys.path.insert(0, '.')
from seqs import REAL

OE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../../references/oeis/')
IDX = {v: i for i, v in enumerate(REAL)}
rows = []
with gzip.open(OE + 'stripped.gz', 'rt') as f:
    for line in f:
        if not line.startswith('A'):
            continue
        an, _, t = line.partition(' ')
        terms = t.strip().strip(',').split(',')
        try:
            s = [int(x) for x in terms if x]
        except ValueError:
            continue
        if len(s) < 3:
            continue
        hit = [IDX.get(v, -1) for v in s]
        if sum(h >= 0 for h in hit) < 3:
            continue
        # prefix and block
        pre, prepos, blk, blki = 0, -1, 0, -1
        for p, h in enumerate(hit):
            if h < 0:
                continue
            r = 0
            while p + r < len(s) and h + r < len(REAL) and s[p + r] == REAL[h + r]:
                r += 1
            if h == 0 and r > pre:
                pre, prepos = r, p
            if r > blk:
                blk, blki = r, h
        # LIS on a-indices
        tails = []
        for h in hit:
            if h >= 0:
                k = bisect.bisect_left(tails, h)
                if k == len(tails):
                    tails.append(h)
                else:
                    tails[k] = h
        lcs = len(tails)
        last = REAL[tails[-1]] if tails else 0
        nin = sum(1 for v in s if 2 <= v <= last)
        rows.append((an, pre, prepos, blk, blki, lcs, lcs / max(nin, 1), len(s), s[:]))

wanted = set()
def top(key, k):
    return sorted(rows, key=key, reverse=True)[:k]
T1 = top(lambda r: (r[1], -r[2]), 60)
T2 = top(lambda r: (r[3], r[6]), 60)
T3 = top(lambda r: (r[5], r[6]), 60)
T4 = top(lambda r: (r[5] if r[6] >= 0.5 else 0, r[6]), 60)
for T in (T1, T2, T3, T4):
    wanted |= {r[0] for r in T}
names = {}
with gzip.open(OE + 'names.gz', 'rt', encoding='utf-8', errors='replace') as f:
    for line in f:
        if line[:7] in wanted:
            names[line[:7]] = line[8:].strip()
REL = re.compile(r'mersenne|perfect|lucas.?lehmer|2\^p *- *1|2\^n *- *1 is prime|2\^k *- *1 is prime|M\(|M_', re.I)


def rel(an):
    return an == 'A000043' or bool(REL.search(names.get(an, '')))


def diverge(r):
    an, pre, pos, *_ , s = r
    nxt = s[pos + pre] if pos + pre < len(s) else None
    return (f'agrees on a_1..a_{pre} (= ..{REAL[pre - 1]}); next term {nxt} vs a_{pre + 1} = {REAL[pre]}'
            if nxt is not None else f'agrees on a_1..a_{pre}; stored terms end there')


out = []
def section(title, T, fmt, n=15):
    out.append('\n== ' + title)
    k = 0
    relat = []
    for r in T:
        if rel(r[0]):
            relat.append(r); continue
        out.append(fmt(r) + '  ' + names.get(r[0], '')[:110]); k += 1
        if k >= n: break
    out.append('  relatives excluded: ' + '; '.join(f'{r[0]} ({names.get(r[0], "")[:50]})' for r in relat[:12]))

out.append(f'sequences scanned with >= 3 terms of a: {len(rows)}')
section('longest prefix a_1..a_k as consecutive terms', T1,
        lambda r: f'{r[0]} k={r[1]:2d} (at position {r[2]}): ' + diverge(r))
section('longest consecutive block a_i..a_j anywhere', T2,
        lambda r: f'{r[0]} block={r[3]:2d} from a_{r[4] + 1}={REAL[r[4]]}  inorder={r[5]} prec={r[6]:.2f}')
section('most terms of a in order (LCS)', T3,
        lambda r: f'{r[0]} inorder={r[5]:2d} prec={r[6]:.2f} nterms={r[7]}')
section('most terms of a in order, among sequences with prec >= 0.5', T4,
        lambda r: f'{r[0]} inorder={r[5]:2d} prec={r[6]:.2f} prefix={r[1]}')
open('oeis.txt', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
