import re, zlib, sys
data = open(sys.argv[1], 'rb').read()
out = []
for m in re.finditer(rb'stream\r?\n(.*?)\r?\nendstream', data, re.S):
    try:
        s = zlib.decompress(m.group(1))
    except Exception:
        continue
    # extract text from TJ/Tj
    txt = []
    for t in re.finditer(rb'\[(.*?)\]\s*TJ|\((.*?)\)\s*Tj', s, re.S):
        if t.group(1) is not None:
            parts = re.findall(rb'\(((?:\\.|[^\\)])*)\)|(-?\d+\.?\d*)', t.group(1))
            w = b''
            for p, num in parts:
                if p: w += p
                elif num and float(num) < -200: w += b' '
            txt.append(w)
        else:
            txt.append(t.group(2))
    out.append(b' '.join(txt))
full = b'\n'.join(out).decode('latin1')
for kw in sys.argv[2:]:
    for mm in re.finditer(kw, full):
        print('...', full[max(0, mm.start()-700):mm.end()+900].replace('\n', ' '), '\n=====')
