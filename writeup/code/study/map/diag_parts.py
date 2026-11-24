"""diagnostic: LR slope (sum over directions, q=5 k=1, y=17) for different part schemes of the response stretch 5e9"""
import sys, math, numpy as np
sys.path.insert(0, '.')
import map as M
S = M.MSet('new=val_new'); q, k = 5, 1
G = M.group(q); C = G.C
Q = 510510; RU, _ = M.rough_index(Q, 17); NR = len(RU)
for parts, CH, cont in ((2, 5*10**6, True), (4, 5*10**6, False), (4, 5*10**6, True), (10, 10**7, False), (20, 5*10**6, False), (8, 10**6, False)):
    acc = {'n': np.zeros((parts, NR), np.int64), 1: np.zeros((parts, NR*C*C), np.int32)}
    M.count_stretch(5*10**9, q, [1], parts, acc, CH=CH, contiguous=cont)
    r = M.lr_project(q, 1, acc, math.log(5.05e9), S.basis(q, 1), S.lam, parts, ys=(17,), rev=(parts == 2))
    W = np.array(r['y17']['W']); Wjk = np.array(r['y17']['Wjk'])
    sl = (W[2] - W[0]).mean() / math.log(100); sljk = (Wjk[:, 2] - Wjk[:, 0]).mean(1) / math.log(100)
    se = math.sqrt((parts - 1) / parts * ((sljk - sljk.mean()) ** 2).sum())
    msg = 'parts %2d chunk %8d contiguous %d: slope/dof %+.5f +- %.5f  W(1e4,1e5,1e6) %s  V %.4f' % (parts, CH, cont, sl, se, np.round(W.mean(1), 5), np.mean(r['y17']['V']))
    if 'rev' in r:
        Wr = np.array(r['rev']['W']); msg += ' | pooled: slope %+.5f V %.4f' % ((Wr[2] - Wr[0]).mean() / math.log(100), np.mean(r['rev']['V']))
    print(msg, flush=True)
