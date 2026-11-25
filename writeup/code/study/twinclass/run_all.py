"""Decomposition for all seen windows: whole windows and their nine sub-windows.  -> all_decomp.json"""
import sys; sys.path.insert(0, '/tmp/claude-1000/-home-bas-praime/86d93f3a-2189-4d59-981e-63b2fd6b59ef/scratchpad/twinclass')
import json, time
import decomp
out = {}
t = time.time()
for name, path, subs in (('1e9-1e10', 'win_1e9.npz', True), ('1e8-1e9', 'win_1e8.npz', True), ('1e7-1e8', 'win_1e7.npz', True),
                         ('A', 'win_A.npz', False), ('B', 'win_B.npz', False), ('C', 'win_C.npz', False)):
    out[name] = decomp.analyse(path, (2,))
    print(name, 'whole done %.0fs' % (time.time() - t), flush=True)
    if subs:
        for s in range(9):
            out['%s/sub%d' % (name, s)] = decomp.analyse(path, (2,), sub=s)
        print(name, 'subs done %.0fs' % (time.time() - t), flush=True)
json.dump(out, open('all_decomp.json', 'w'), indent=1)
