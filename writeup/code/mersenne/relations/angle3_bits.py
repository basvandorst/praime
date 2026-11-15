"""Angle 3: the terms as bit strings.

lz_cat     LZ76 complexity (Kaspar-Schuster) of the concatenated MSB-first binary
           strings, normalised c*log2(n)/n                                  (lo = structure)
lz_nolead  same with the forced leading 1 of every term removed           (lo)
dfa_var    minimal DFA of the finite language {bin(a) : a < 2^k}, MSB first,
           summed over k = 8..20, divided by the summed trie sizes          (lo)
dfa_lsb    same with LSB-first strings                                      (lo)
dfa_fix    minimal DFA for {a < 2^k} as k-bit strings (zero padded), ratio   (lo)
suf_chi_L  chi-square of the last L bits of the odd terms >= 2^L against
           uniform on odd residues, L = 2..6                                (hi)
pre_chi_L  chi-square of the first L bits (leading 1 included) of terms with
           >= L+2 bits against binary Benford, L = 2..5                     (hi)
lcs2,lcs3  longest binary substring found in >= 2 (>= 3) distinct terms     (hi)
lcp2,lsf2  longest common prefix / suffix shared by two distinct terms      (hi)
"""
import json, sys, time
import numpy as np


def lz76(s):
    n = len(s)
    i, k, l = 0, 1, 1
    c, kmax = 1, 1
    while True:
        if s[i + k - 1] == s[l + k - 1]:
            k += 1
            if l + k > n:
                c += 1
                break
        else:
            if k > kmax:
                kmax = k
            i += 1
            if i == l:
                c += 1
                l += kmax
                if l + 1 > n:
                    break
                i, k, kmax = 0, 1, 1
            else:
                k = 1
    return c


def min_dfa(strings):
    """States of the minimal DFA of a finite language (dead state excluded) and trie size."""
    trie = [{}]
    acc = [False]
    for s in strings:
        v = 0
        for ch in s:
            nxt = trie[v].get(ch)
            if nxt is None:
                trie.append({})
                acc.append(False)
                nxt = len(trie) - 1
                trie[v][ch] = nxt
            v = nxt
        acc[v] = True
    # post-order canonical ids
    canon = {}
    cid = [None] * len(trie)
    order = []
    stack = [(0, False)]
    while stack:
        v, done = stack.pop()
        if done:
            order.append(v)
        else:
            stack.append((v, True))
            for ch in trie[v]:
                stack.append((trie[v][ch], False))
    for v in order:
        sig = (acc[v], cid[trie[v]['0']] if '0' in trie[v] else -1, cid[trie[v]['1']] if '1' in trie[v] else -1)
        if sig not in canon:
            canon[sig] = len(canon)
        cid[v] = canon[sig]
    return len(canon), len(trie)


def bits(a):
    return [bin(x)[2:] for x in a]


def chi_suffix(a, L):
    vals = [x % (1 << L) for x in a if x % 2 == 1 and x >= (1 << L)]
    m = len(vals)
    cats = list(range(1, 1 << L, 2))
    exp = m / len(cats)
    cnt = np.array([vals.count(c) for c in cats])
    return float(((cnt - exp) ** 2 / exp).sum()), dict(zip(cats, cnt.tolist()))


def chi_prefix(a, L):
    vals = [x >> (x.bit_length() - L) for x in a if x.bit_length() >= L + 2]
    m = len(vals)
    cats = list(range(1 << (L - 1), 1 << L))
    pr = np.array([np.log2((c + 1) / c) for c in cats])
    exp = m * pr
    cnt = np.array([vals.count(c) for c in cats])
    return float(((cnt - exp) ** 2 / exp).sum()), dict(zip(cats, cnt.tolist()))


def lcs_multi(strs, need):
    best = 0
    maxlen = max(len(s) for s in strs)
    for L in range(1, maxlen + 1):
        seen = {}
        found = False
        for idx, s in enumerate(strs):
            subs = {s[i:i + L] for i in range(len(s) - L + 1)}
            for u in subs:
                seen[u] = seen.get(u, 0) + 1
                if seen[u] >= need:
                    found = True
        if found:
            best = L
        else:
            break
    return best


def common_prefix_suffix(strs):
    lcp = lsf = 0
    for i in range(len(strs)):
        for j in range(i + 1, len(strs)):
            s, t = strs[i], strs[j]
            p = 0
            while p < min(len(s), len(t)) and s[p] == t[p]:
                p += 1
            q = 0
            while q < min(len(s), len(t)) and s[-1 - q] == t[-1 - q]:
                q += 1
            lcp = max(lcp, p)
            lsf = max(lsf, q)
    return lcp, lsf


def all_stats(a):
    out = {}
    b = bits(a)
    cat = ''.join(b)
    n = len(cat)
    out['lz_cat'] = lz76(cat) * np.log2(n) / n
    nl = ''.join(s[1:] for s in b)
    out['lz_nolead'] = lz76(nl) * np.log2(len(nl)) / len(nl)
    sv = tv = sl = tl = sf = tf = 0
    for k in range(8, 21):
        sub = [x for x in a if x < (1 << k)]
        s1, t1 = min_dfa(bits(sub))
        s2, t2 = min_dfa([s[::-1] for s in bits(sub)])
        s3, t3 = min_dfa([format(x, '0%db' % k) for x in sub])
        sv += s1; tv += t1; sl += s2; tl += t2; sf += s3; tf += t3
    out['dfa_var'] = sv / tv
    out['dfa_lsb'] = sl / tl
    out['dfa_fix'] = sf / tf
    for L in range(2, 7):
        out['suf_chi_%d' % L] = chi_suffix(a, L)[0]
    for L in range(2, 6):
        out['pre_chi_%d' % L] = chi_prefix(a, L)[0]
    out['lcs2'] = lcs_multi(b, 2)
    out['lcs3'] = lcs_multi(b, 3)
    out['lcp2'], out['lsf2'] = common_prefix_suffix(b)
    return out


SIDES = {'lz_cat': 'lo', 'lz_nolead': 'lo', 'dfa_var': 'lo', 'dfa_lsb': 'lo', 'dfa_fix': 'lo',
         'lcs2': 'hi', 'lcs3': 'hi', 'lcp2': 'hi', 'lsf2': 'hi'}
SIDES.update({'suf_chi_%d' % L: 'hi' for L in range(2, 7)})
SIDES.update({'pre_chi_%d' % L: 'hi' for L in range(2, 6)})

if __name__ == '__main__':
    D = json.load(open('nulls.json'))
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    t0 = time.time()
    A = D['real']
    real = all_stats(A)
    with open('angle3_real.txt', 'w') as rep:
        rep.write('real stats: %s\n' % real)
        for k in range(8, 21):
            sub = [x for x in A if x < (1 << k)]
            rep.write('k=%2d terms=%2d minDFA(var,MSB)=%s minDFA(LSB)=%s minDFA(fixed k bits)=%s\n'
                      % (k, len(sub), min_dfa(bits(sub)), min_dfa([s[::-1] for s in bits(sub)]),
                         min_dfa([format(x, '0%db' % k) for x in sub])))
        for L in range(2, 7):
            rep.write('suffix L=%d counts %s\n' % (L, chi_suffix(A, L)[1]))
        for L in range(2, 6):
            rep.write('prefix L=%d counts %s\n' % (L, chi_prefix(A, L)[1]))
    out = {s: {'side': SIDES[s], 'real': real[s], 'a': [], 'b': []} for s in SIDES}
    for fam in ('a', 'b'):
        for seq in D[fam][:nmax]:
            st = all_stats(seq)
            for s in SIDES:
                out[s][fam].append(st[s])
    json.dump(out, open('angle3.json', 'w'))
    print('done', time.time() - t0)
