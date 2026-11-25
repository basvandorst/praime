import sys, pypdf
for f in sys.argv[1:]:
    r = pypdf.PdfReader(f)
    with open(f[:-4] + '.txt', 'w') as o:
        for i, p in enumerate(r.pages):
            o.write('\n=== page %d ===\n' % (i+1))
            try: o.write(p.extract_text())
            except Exception as e: o.write('[err %s]' % e)
