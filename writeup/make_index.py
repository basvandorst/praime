# Rebuild writeup/INDEX.md: every puzzle and every conjecture with its status and a link to our write-up.
# Run from the repo root: python3 writeup/make_index.py
import os, re

titles, ctitles = {}, {}
last = None
for line in open('data/INDEX.txt', encoding='utf-8'):
    line = line.rstrip('\n')
    if line.startswith('puzzle-'):
        _, url, title = line.split('\t', 2)
        last = (titles, int(re.search(r'puzz_(\d+)', url).group(1)))
        titles[last[1]] = title.strip()
    elif line.startswith('conjecture-'):
        _, url, title = line.split('\t', 2)
        last = (ctitles, int(re.search(r'conj_(\d+)', url).group(1)))
        ctitles[last[1]] = title.strip()
    elif line.startswith('problem-') or line.lstrip().startswith('->'):
        last = None
    elif last is not None and line.strip():
        last[0][last[1]] += ' ' + line.strip()   # title wrapped onto the next line

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
# Where each README row links to (usually writeup/N.md; puzzle 900 points to its conjecture write-up).
target = {}
for row in readme.splitlines():
    m = re.match(r'\| \[(\d+)[^\]]*\]\(https://www\.primepuzzles\.net/puzzles/', row)
    t = re.findall(r'\]\(writeup/([^)]+\.md)\)', row)
    if m and t:
        target.setdefault(int(m.group(1)), t[-1])

def name(n):
    # drop the "Puzzle 123.- " prefix, but only when it is this puzzle's own number (or off by one)
    t = titles[n]
    m = re.match(r'^P?uzzle\s*(?:No\.\s*)?(\d+)', t)
    if m and int(m.group(1)) != n:   # stale HTML title (copied from the previous page): use the page heading
        page = open('data/puzzles/puzzle-%04d.txt' % n, encoding='utf-8').read().split('=' * 20, 1)[-1]
        h = re.search(r'^Puzzle\s*%d\b.*$' % n, page, re.M)
        t = h.group(0).strip() if h else t
    m = re.match(r'^P?uzzle\s*(?:No\.\s*)?(\d+)\s*[.:-]*\s*', t)
    if m and abs(int(m.group(1)) - n) <= 1 and t[m.end():].strip():
        t = t[m.end():]
    return t.strip().replace('|', '\\|')

rows = []
for n in sorted(titles):
    page = 'https://www.primepuzzles.net/puzzles/puzz_%03d.htm' % n
    if n in beaten or n in answered:
        status = '[%s](%s)' % ('beaten' if n in beaten else 'answered', target.get(n, '%d.md' % n))
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


# Conjectures: status is the Status line of writeup/conj-N.md up to the first '(' or '·'.
# A Status line that starts with 'open' marks an attempt; any other status is a result.
def cname(n):
    t = ctitles[n]
    m = re.match(r'^Conjecture\s*(\d+)', t)
    if m and int(m.group(1)) != n:   # the site's title is wrong for some pages (conj_101 repeats 100): use the page heading
        page = open('data/puzzles/conjecture-%04d.txt' % n, encoding='utf-8').read().split('=' * 20, 1)[-1]
        h = re.search(r'^Conjecture\s*%d\b.*$' % n, page, re.M)
        t = h.group(0).strip() if h else t
    t2 = re.sub(r'^Conjecture\s*\d+\s*[.:-]*\s*', '', t).strip()
    return (t2 or t).replace('|', '\\|')

crows, settled, cattempts = [], 0, 0
for n in sorted(ctitles):
    page = 'https://www.primepuzzles.net/conjectures/conj_%03d.htm' % n
    path = 'writeup/conj-%d.md' % n
    status = 'open'
    if os.path.exists(path):
        m = re.search(r'^\*\*Status:\*\*\s*([a-z][a-z ]*?)\s*(?:\(|·|$)', open(path, encoding='utf-8').read(), re.M)
        word = m.group(1) if m else 'answered'
        if word.startswith('open'):
            status, cattempts = '[open](conj-%d.md)' % n, cattempts + 1
        else:
            status, settled = '[%s](conj-%d.md)' % (word, n), settled + 1
    crows.append('| [%d](%s) | %s | %s |' % (n, page, cname(n), status))

with open('writeup/INDEX.md', 'a', encoding='utf-8') as f:
    f.write('\n# All conjectures\n\n')
    f.write('%d conjectures: %d with a result of ours (proved, disproved, explained, answered, or a variant settled), '
            '%d open (%d with a write-up of our attempt). Status links go to our write-up.\n\n'
            % (len(ctitles), settled, len(ctitles) - settled, cattempts))
    f.write('| # | Conjecture | Status |\n|---|---|---|\n')
    f.write('\n'.join(crows) + '\n')
