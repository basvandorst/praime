# Puzzle 896 (writeup/896.md): is there a digit string of length L containing every word of the set as a subsequence?
# Same encoding as p896_verify.py (subsequence automaton per prefix, shared over a trie), solved with OR-Tools CP-SAT.
# Usage: python3 p896_cpsat.py A|B|C|P4|P5 L [time_s]
import sys, time
from sympy import primerange
from ortools.sat.python import cp_model

which, L = sys.argv[1], int(sys.argv[2]); tl = float(sys.argv[3]) if len(sys.argv) > 3 else 1800
A = '4242064951 5519385407 6886477529 3816749731 4700006089 9650560979 4371455527 5419320281 8761823531 3454277693'.split()
words = {'A': A, 'B': [str(p) for p in primerange(2, 100)], 'C': [str(p) for p in primerange(2, 1000)],
         'P4': [str(p) for p in primerange(2, 10**4)], 'P5': [str(p) for p in primerange(2, 10**5)]}[which]
M = cp_model.CpModel()
x = [[M.NewBoolVar('') for d in range(10)] for i in range(L)]
for i in range(L): M.AddExactlyOne(x[i])
T = M.NewConstant(1); F = M.NewConstant(0)
memo = {(): [T] * (L + 1)}
def m(pre):
    if pre in memo: return memo[pre]
    prev = m(pre[:-1]); d = pre[-1]; cur = [F]
    for i in range(1, L + 1):
        v = M.NewBoolVar(''); a, b, xi = cur[i - 1], prev[i - 1], x[i - 1][d]
        # v <-> a or (b and xi)
        M.AddBoolOr([a.Not(), v]); M.AddBoolOr([b.Not(), xi.Not(), v]); M.AddBoolOr([v.Not(), a, b]); M.AddBoolOr([v.Not(), a, xi])
        cur.append(v)
    memo[pre] = cur; return cur
for w in words: M.Add(m(tuple(map(int, w)))[L] == 1)
# symmetry: digit counts must be at least the maximum count in any single word (valid cut)
for d in range(10):
    need = max(w.count(str(d)) for w in words)
    if need: M.Add(sum(x[i][d] for i in range(L)) >= need)
s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl; s.parameters.num_workers = 2
t0 = time.time(); st = s.Solve(M)
print('%s L=%d: %s (%.0f s)' % (which, L, s.StatusName(st), time.time() - t0))
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    out = ''.join(str(next(d for d in range(10) if s.Value(x[i][d]))) for i in range(L))
    def subseq(w, t):
        it = iter(t); return all(c in it for c in w)
    assert all(subseq(w, out) for w in words)
    print(out)
