# Rebuild writeup/INDEX.md: every puzzle with its status and a link to our write-up.
# Run from the repo root: python3 writeup/make_index.py
import os, re

titles = {}
last = None
for line in open('data/INDEX.txt', encoding='utf-8'):
    line = line.rstrip('\n')
    if line.startswith('puzzle-'):
        _, url, title = line.split('\t', 2)
        last = int(re.search(r'puzz_(\d+)', url).group(1))
        titles[last] = title.strip()
    elif re.match(r'(problem|conjecture)-', line) or line.lstrip().startswith('->'):
        last = None
    elif last is not None and line.strip():
        titles[last] += ' ' + line.strip()   # title wrapped onto the next line

# Status from the README tables. Each dated Results section has its own tables;
# a table's header row says whether it lists beaten records or answers.
readme = open('README.md', encoding='utf-8').read()
beaten, answered, kind = set(), set(), None
for row in readme.splitlines():
    if row.lstrip().startswith('| # |'):
        kind = {'Puzzle': beaten, 'Open question': answered}.get(row.split('|')[2].strip())
    for n in re.findall(r'\| \[(\d+)[^\]]*\]\(https://www\.primepuzzles\.net/puzzles/', row):
        if kind is not None:
            kind.add(int(n))
answered -= beaten

def name(n):
    # drop the "Puzzle 123.- " prefix, but only when it is this puzzle's own number (or off by one)
    t = titles[n]
    m = re.match(r'^P?uzzle\s*(?:No\.\s*)?(\d+)\s*[.:-]*\s*', t)
    if m and abs(int(m.group(1)) - n) <= 1 and t[m.end():].strip():
        t = t[m.end():]
    return t.strip().replace('|', '\\|')

rows = []
for n in sorted(titles):
    page = 'https://www.primepuzzles.net/puzzles/puzz_%03d.htm' % n
    if n in beaten or n in answered:
        status = '[%s](../results/%d.md)' % ('beaten' if n in beaten else 'answered', n)
    elif os.path.exists('writeup/%d.md' % n):
        status = '[open](%d.md)' % n
    else:
        status = 'open'
    rows.append('| [%d](%s) | %s | %s |' % (n, page, name(n), status))

attempts = sum(os.path.exists('writeup/%d.md' % n) for n in titles if n not in beaten | answered)
with open('writeup/INDEX.md', 'w', encoding='utf-8') as f:
    f.write('# All puzzles\n\n')
    f.write('%d puzzles: %d beaten, %d answered, %d open (%d with a write-up of our attempt). '
            'Status links go to our write-up. Rebuild with `python3 writeup/make_index.py`.\n\n'
            % (len(titles), len(beaten), len(answered), len(titles) - len(beaten) - len(answered), attempts))
    f.write('| # | Puzzle | Status |\n|---|---|---|\n')
    f.write('\n'.join(rows) + '\n')
