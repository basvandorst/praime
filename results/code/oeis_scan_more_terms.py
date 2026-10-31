# Scan: puzzle pages that ask for more terms and whose list is a prefix of an OEIS sequence (references/oeis/stripped.gz).
# Print the page, the A-number and the next OEIS terms. Filter the hits by hand: keep only those whose next term is not on the page.
# Pages that ask for more terms and whose list of terms matches an OEIS sequence: does OEIS have terms beyond the page?
import gzip, re, glob, os
seqs = {}
with gzip.open('/home/bas/praime/references/oeis/stripped.gz', 'rt') as f:
    for line in f:
        if line.startswith('A'):
            a, s = line.split(' ', 1); seqs[a] = [int(x) for x in s.strip().strip(',').split(',') if x.lstrip('-').isdigit()]
names = {}
with gzip.open('/home/bas/praime/references/oeis/names.gz', 'rt') as f:
    for line in f:
        if line.startswith('A'): a, _, n = line.rstrip().partition(' '); names[a] = n
done = {os.path.splitext(f)[0] for f in os.listdir('/home/bas/praime/results')} | {os.path.splitext(f)[0] for f in os.listdir('/home/bas/praime/writeup')}
for fn in sorted(glob.glob('/home/bas/praime/data/puzzles/puzzle-*.txt')):
    num = str(int(re.search(r'(\d{4})', fn).group(1)))
    if num in done: continue
    txt = open(fn, encoding='utf-8', errors='replace').read().split('=' * 20, 1)[-1]
    if not re.search(r'more terms|next term|extend|find more|the next', txt, re.I): continue
    best = None
    for m in re.finditer(r'(?<![\d.])(\d{1,15}(?:\s*,\s*\d{1,15}){4,})', txt):
        t = [int(x) for x in re.split(r'\s*,\s*', m.group(1))]
        if any(b <= a for a, b in zip(t, t[1:])) or t[-1] < 100: continue
        for a, s in seqs.items():
            if len(s) > len(t) and s[:len(t)] == t:          # page list is a prefix of the OEIS terms
                pagemax = max(int(x.replace(',', '')) for x in re.findall(r'\d[\d,]{0,40}\d', txt) if len(x) < 45 and x.replace(',', '').isdigit()) if False else None
                extra = s[len(t):]
                if best is None or len(extra) > len(best[2]): best = (a, t, extra)
    if best:
        a, t, extra = best
        print(f"puzzle {num}: page list {t[:4]}..{t[-1]} ({len(t)} terms) = prefix of {a} ({names.get(a,'')[:70]}); OEIS has {len(extra)} more, next {extra[:3]}")
