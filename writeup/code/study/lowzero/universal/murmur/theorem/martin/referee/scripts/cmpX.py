import json, math
for fn in ['rows_1500.json', 'rows.json']:
    rows = json.load(open('/tmp/claude-1000/referee_martin/' + fn))
    out = []
    for xt in sorted(set(r[0] for r in rows)):
        sel = [r for r in rows if r[0] == xt]
        out.append('x=%.1f fixed rms %.3f, rms(A-M) %.3f, chi2/n %.2f' % (xt, math.sqrt(sum(r[5]**2 for r in sel)/len(sel)),
                   math.sqrt(sum((r[4]-r[7])**2 for r in sel)/len(sel)), sum(((r[4]-r[7])/r[6])**2 for r in sel)/len(sel)))
    allf = math.sqrt(sum(r[5]**2 for r in rows)/len(rows))
    print(fn, '; '.join(out), '| overall fixed rms %.3f' % allf)
