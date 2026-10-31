# Proof for every length d: finite automaton over symmetric digit pairs (outside in) for S = M + rev(M), M ending in 9.
# Pair i: right column i, left column d-1-i, pair sum c_i = m_i + m_{d-1-i} (c_0 = 9 + f, f = first digit 1..9).
# r_i = carry into right column i (r_0 = 0); l_i = carry into left column d-1-i (from the inside, guessed).
# Right: R = c + r_i, digit a_i = R % 10, r_{i+1} = R // 10.  Left: Lv = c + l_i, digit b_i = Lv % 10, carry out Lv // 10.
# Palindrome S_i = S_{d-i}:  i = 0: a_0 = carry out of column d-1;  i >= 1: a_i = b_{i-1}, and carry out of left
# column d-1-i must equal l_{i-1}. State after pair i: (r_{i+1}, l_i, b_i).
from itertools import product
def pair_options(i0):
    return range(10, 19) if i0 else range(0, 19)      # c_0 = 9 + f with f >= 1
# initial pairs
start = {}
for f in range(1, 10):
    c = 9 + f
    for l0 in (0, 1):
        R = c; a0 = R % 10; r1 = R // 10; Lv = c + l0; b0 = Lv % 10
        if a0 == Lv // 10:                              # S_0 = S_d
            start.setdefault((r1, l0, b0), set()).add(('f', f))
print("start states:", start)
# explore reachable states and allowed inner pair sums
reach = set(start); frontier = list(start); allowed = {}
while frontier:
    st = frontier.pop(); r, lprev, bprev = st
    for c, l in product(range(0, 19), (0, 1)):
        R = c + r; a = R % 10; rn = R // 10; Lv = c + l; b = Lv % 10
        if a == bprev and Lv // 10 == lprev:
            allowed.setdefault(st, set()).add(c)
            ns = (rn, l, b)
            if ns not in reach: reach.add(ns); frontier.append(ns)
print("reachable states:", sorted(reach))
print("allowed inner pair sums per state:", {s: sorted(v) for s, v in allowed.items()})
# termination: even d (innermost pair): l_{h-1} == r_h ; odd d: middle digit m with 2m + r_h -> digit == b_{h-1}, carry == l_{h-1}
even_ok = [s for s in reach if s[1] == s[0]]
odd_mid = {s: [m for m in range(10) if (2 * m + s[0]) % 10 == s[2] and (2 * m + s[0]) // 10 == s[1]] for s in reach}
print("even end accepted in states:", even_ok)
print("odd end: allowed middle digits per state:", odd_mid)
