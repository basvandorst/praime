# Rebuild the "## Status" section of README.md: every puzzle, conjecture and problem by status.
# Beaten, solved and partially solved numbers link to our write-up, unsolved ones to the page on primepuzzles.net.
# Run from the repo root: python3 writeup/code/make_status.py
import os, re

pages = {'puzzle': [], 'conjecture': [], 'problem': []}
for line in open('data/INDEX.txt', encoding='utf-8'):
    m = re.match(r'^(puzzle|conjecture|problem)-\d+\t(\S+)', line)
    if m:
        pages[m.group(1)].append((int(re.search(r'_(\d+)\.htm', m.group(2)).group(1)), m.group(2)))

# Puzzle status comes from the README tables; for a write-up without a row there, its Status line decides. A table's header row says whether it lists beaten records or answers.
readme = open('README.md', encoding='utf-8').read()
beaten, answered, target, kind = set(), set(), {}, None
for row in readme.splitlines():
    if row.lstrip().startswith('| # |'):
        kind = {'Puzzle': beaten, 'Open question': answered}.get(row.split('|')[2].strip())
    m = re.match(r'\| \[(\d+)[^\]]*\]\(https://www\.primepuzzles\.net/puzzles/', row.lstrip())
    if m and kind is not None:
        kind.add(int(m.group(1)))
    t = re.findall(r'\]\((writeup/[^)]+\.md)\)', row)
    if m and t:
        target.setdefault(int(m.group(1)), t[-1])
answered -= beaten


def status_word(path):
    # The Status line of a write-up, up to the first '(' or '·'. 'open' marks an attempt.
    m = re.search(r'^\*\*Status:\*\*\s*([a-z][a-z ]*?)\s*(?:\(|·|$)', open(path, encoding='utf-8').read(), re.M)
    return m.group(1) if m else 'answered'


def classify(kind, n):
    if kind == 'puzzle':
        if n in beaten:
            return 'Beaten', target.get(n, 'writeup/%d.md' % n)
        if n in answered:
            return 'Solved', target.get(n, 'writeup/%d.md' % n)
        path = 'writeup/%d.md' % n
    else:
        path = 'writeup/%s-%d.md' % ({'conjecture': 'conj', 'problem': 'prob'}[kind], n)
    if not os.path.exists(path):
        return 'Unsolved', None
    return ('Partially solved' if status_word(path).startswith('open') else 'Solved'), path


section = ['## Status', '']
for kind, title in (('puzzle', 'Puzzles'), ('conjecture', 'Conjectures'), ('problem', 'Problems')):
    groups = {'Beaten': [], 'Solved': [], 'Partially solved': [], 'Unsolved': []}
    for n, url in sorted(pages[kind]):
        group, link = classify(kind, n)
        groups[group].append('[%d](%s)' % (n, link or url))
    counts = ', '.join('%d %s' % (len(v), k.lower()) for k, v in groups.items() if v)
    section += ['### %s (%d): %s' % (title, len(pages[kind]), counts), '']
    for group, links in groups.items():
        if links:
            section += ['**%s:** %s' % (group, ', '.join(links)), '']
section.append('Unsolved numbers link to the page on primepuzzles.net, all others to our write-up. '
               'Rebuild with `python3 writeup/code/make_status.py`.')
block = '\n'.join(section) + '\n\n'

if re.search(r'^## Status\n', readme, re.M):
    readme = re.sub(r'^## Status\n.*?(?=^## )', lambda _: block, readme, count=1, flags=re.M | re.S)
else:
    readme = re.sub(r'^## Results ', lambda m: block + m.group(0), readme, count=1, flags=re.M)
open('README.md', 'w', encoding='utf-8').write(readme)
print('\n'.join(line for line in section if line.startswith('### ')))
