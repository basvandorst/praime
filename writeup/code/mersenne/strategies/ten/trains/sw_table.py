"""(The HTML is not stored in the repository.) Read Table 1 of Sorenson and Waiss (arXiv:2509.00236v1, saved HTML): h(k, x) = #{n <= x : f(n) = k} for each column
x, and compute the mean of f, sum_k k h(k, x) / sum_k h(k, x).   usage: sw_table.py sw.html"""
import sys, re, html
t = open(sys.argv[1], encoding='utf-8', errors='replace').read()
i = t.find('Table 1'); k2 = t.find('Table 2', i + 10); seg = t[i:k2 if k2 > 0 else i + 400000]
rows = re.findall(r'<tr[^>]*>(.*?)</tr>', seg, re.S)
cells = [[html.unescape(re.sub(r'<[^>]+>', ' ', c)).strip() for c in re.findall(r'<t[dh][^>]*>(.*?)</t[dh]>', r, re.S)] for r in rows]
head = []; data = {}
for c in cells:
    c = [re.sub(r'\s+', '', v) for v in c]
    if c and c[0] == 'Count': head += c[1:]; continue
    if c and re.fullmatch(r'\d+', c[0] or 'x'):
        data.setdefault(int(c[0]), []).extend(int(v) if re.fullmatch(r'\d+', v) else 0 for v in c[1:])
ncol = max(len(v) for v in data.values())
for j in range(ncol):
    col = [data[k][j] if j < len(data[k]) else 0 for k in sorted(data)]
    tot = sum(col); s = sum(k * v for k, v in zip(sorted(data), col))
    m = re.search(r'\^\{(\d+)\}', head[j]) if j < len(head) else None
    if not (tot and m): continue
    e = int(m.group(1))
    if abs(tot - 10 ** e) > 1e-5 * 10 ** e: continue          # parts of other tables mixed into the column: skip
    print('x = 10^%-2d  sum of h = %d  mean f = %.7f' % (e, tot, s / tot))
# the first block of the table (10^2 ... 10^9) on its own, to compare with our own counts
seg1 = t[t.find('Table 1'):]; seg1 = seg1[:seg1.find('</table>')]
rows1 = [[re.sub(r'\s+', '', html.unescape(re.sub(r'<[^>]+>', ' ', c))) for c in re.findall(r'<t[dh][^>]*>(.*?)</t[dh]>', r, re.S)] for r in re.findall(r'<tr[^>]*>(.*?)</tr>', seg1, re.S)]
hd = rows1[0][1:]; dat = {int(r[0]): [int(v) if v.isdigit() else 0 for v in r[1:]] for r in rows1[1:] if r and r[0].isdigit()}
for j, h in enumerate(hd):
    m = re.search(r'\^\{(\d+)\}', h)
    col = [dat[k][j] if j < len(dat[k]) else 0 for k in sorted(dat)]; tot = sum(col)
    if m and tot == 10 ** int(m.group(1)): print('first block, x = 10^%-2s mean f = %.7f' % (m.group(1), sum(k * v for k, v in zip(sorted(dat), col)) / tot))

