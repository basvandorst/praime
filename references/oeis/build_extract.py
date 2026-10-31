"""Build references/oeis/puzzle-sequences.txt and -index.txt from the .seq files in references/oeis/oeisdata.
Run clone_oeisdata.sh first. Usage: CAP=1000 python3 references/oeis/build_extract.py"""
import os, re, subprocess, sys

S = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(S, '..', '..'))
OUT = os.path.join(REPO, 'references/oeis')
SEQ = os.path.join(OUT, 'oeisdata', 'seq')
cited = sorted(set(re.findall(r'A\d{6}', ''.join(open(os.path.join(REPO, 'data/puzzles', f), encoding='utf-8', errors='replace').read()
    for f in os.listdir(os.path.join(REPO, 'data/puzzles')) if f.endswith('.txt')))))
ids = [a for a in cited if os.path.exists(os.path.join(SEQ, a[:4], a + '.seq'))]
notseq = [a for a in cited if a not in ids]

def cites(aid):
    r = subprocess.run(['grep', '-lE', aid + '([^0-9]|$)'] + sorted(
        os.path.join(REPO, 'data/puzzles', f) for f in os.listdir(os.path.join(REPO, 'data/puzzles'))
        if f.endswith('.txt')), capture_output=True, text=True)
    return [os.path.basename(x) for x in r.stdout.split()]

# real A-numbers that only occur as the head of a longer digit string (e.g. A0000058 in problem 83)
prefixonly = [a for a in ids if not cites(a)]
ids = [a for a in ids if a not in prefixonly]

def parse(aid):
    p = os.path.join(SEQ, aid[:4], aid + '.seq')
    if not os.path.exists(p):
        return None
    d = {'S': [], 'F': [], 'C': [], 'N': None, 'O': None, 'K': None, 'I': None}
    for line in open(p, encoding='utf-8', errors='replace'):
        line = line.rstrip('\n')
        m = re.match(r'^%([A-Za-z]) (A\d{6}) ?(.*)$', line)
        if not m or m.group(2) != aid:
            continue
        t, rest = m.group(1), m.group(3)
        if t in 'STU':
            d['S'].append(rest.strip())
        elif t == 'F':
            d['F'].append(rest)
        elif t == 'C':
            d['C'].append(rest)
        elif t in 'NOKI':
            d[t] = rest
    if d['N'] is None:
        return None
    terms = ','.join(x.strip(',') for x in d['S'] if x.strip(','))
    d['terms'] = [x for x in terms.split(',') if x != '']
    return d

entries, failed = {}, []
for a in ids:
    d = parse(a)
    if d is None:
        failed.append(a)
    else:
        entries[a] = d
nF = sum(1 for d in entries.values() if d['F'])

CAP = int(os.environ.get('CAP', '0'))      # cap term line at CAP chars (0 = no cap)
NOC = os.environ.get('NOCOMMENTS') == '1'
ncut = 0
blocks = []
for a in sorted(entries):
    d = entries[a]
    b = [a, d['N']]
    tl = ', '.join(d['terms'])
    if CAP and len(tl) > CAP:
        tl = tl[:CAP].rsplit(', ', 1)[0] + ', ...'
        ncut += 1
    b.append(tl)
    b += ['F: ' + f for f in d['F']]
    b.append('K: ' + (d['K'] or ''))
    b.append('offset: ' + (d['O'] or ''))
    for c in ([] if NOC else d['C'][:3]):
        b.append('C: ' + (c if len(c) <= 300 else c[:297] + '...'))
    blocks.append('\n'.join(b))

COMMIT = subprocess.run(['git', '-C', os.path.join(OUT, 'oeisdata'), 'log', '-1', '--format=%h, %cd'],
                        capture_output=True, text=True).stdout.strip()
NOTICE = open(os.path.join(S, "notice.txt")).read()
header = NOTICE + f"""# OEIS sequences cited in the primepuzzles.net pages (data/puzzles/*.txt)
# Made from the .seq files (OEIS internal format) in references/oeis/oeisdata,
#   a sparse clone of https://github.com/oeis/oeisdata (commit {COMMIT}).
# List made with: grep -ohE "A[0-9]{{6}}" data/puzzles/*.txt | sort -u  ({len(cited)} hits);
#   {len(notseq)} hits are digit strings in puzzle text, not sequences, and are left out:
#   {' '.join(notseq)}
#   {len(prefixonly)} more are real A-numbers but occur only as the head of a longer digit string
#   (such as A0000058 in problem-0083), so they are not citations and are left out too:
#   {' '.join(prefixonly)}
# Kept per entry: %N name, the whole term list (%S/%T/%U), all %F formulas, %K keywords,
#   %O offset, at most the first 3 %C comments (each cut to 300 characters).
#   Dropped: programs, references, links, examples, cross-references.
# Block layout: A-number / name / terms (one line, comma-separated) / F: lines / K: / offset: / C: lines
# Counts: {len(entries)} sequences ({len(failed)} unreadable{': ' + ', '.join(failed) if failed else ''}),
#   {nF} entries with at least one %F formula.
"""
with open(os.path.join(OUT, 'puzzle-sequences.txt'), 'w') as f:
    f.write(header + (f"# Term lines longer than {CAP} characters are cut and end with ', ...' ({ncut} entries).\n" if CAP else '') + '\n' + '\n\n'.join(blocks) + '\n')

with open(os.path.join(OUT, 'puzzle-sequences-index.txt'), 'w') as f:
    f.write(NOTICE + '# Columns: A-number<TAB>name<TAB>puzzle files citing it (grep -l)\n')
    for a in ids:
        name = entries[a]['N'] if a in entries else '(not fetched)'
        f.write(f"{a}\t{name}\t{' '.join(cites(a))}\n")

print('cut', ncut, 'ok', len(entries), 'failed', failed, 'withF', nF)
