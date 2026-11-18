import json, random, subprocess, sys
def ringpow(n, r, kind):
    if kind == 2:
        h = r
        cur = [0]*h; cur[0] = 1; cur[1] = 1
        def mul(a, b):
            c = [0]*h
            for i in range(h):
                if a[i]:
                    for j in range(h):
                        k = i+j
                        if k >= h: c[k-h] -= a[i]*b[j]
                        else: c[k] += a[i]*b[j]
            return [v % n for v in c]
    else:
        cur = [0]*r
        if kind == 1: cur[0] = 2; cur[1] += 1; cur[r-1] += 1
        else: cur[0] = 1; cur[1] = 1
        def mul(a, b):
            c = [0]*r
            for i in range(r):
                if a[i]:
                    for j in range(r): c[(i+j) % r] += a[i]*b[j]
            return [v % n for v in c]
    base = cur[:]; res = None
    for bit in bin(n)[2:]:
        if res is None: res = base[:]; continue
        res = mul(res, res)
        if bit == '1': res = mul(res, base)
    return res
def liar_full(n, r):
    v = ringpow(n, r, 0); w = [0]*r; w[0] += 1; w[n % r] += 1
    return v == [x % n for x in w]
if __name__ == "__main__":
    random.seed(20261118)
    psp = [int(x) for x in open('psp2_to1e10.txt')]
    known = json.load(open('ring_psp2_to1e10.json'))
    rings = [5,7,9,11,13,15,16,17,19,21,23,25,27,29,31,32,33,35,37,39]
    out = subprocess.run(['./ring3','-fit']+[str(r) for r in rings], input='\n'.join(map(str,psp))+'\n', capture_output=True, text=True)
    got = {}
    for l in out.stdout.split('\n'):
        if l.startswith('ring'): _, r, _, n = l.split(); got.setdefault(int(r), []).append(int(n))
    bad = 0; tot = 0
    for r in rings:
        a = sorted(got.get(r, [])); b = sorted(x for x in known[str(r)] if x % 2 and __import__('math').gcd(x, r) == 1)
        tot += len(b)
        if a != b: bad += 1; print('ring', r, 'differs', len(a), len(b))
    print('known liars below 1e10 on', len(rings), 'rings (numbers that fit):', tot, '| rings where the fast program gives a different list:', bad)
    out2 = subprocess.run(['./ring3','-plain','-fit']+[str(r) for r in rings], input='\n'.join(map(str,psp))+'\n', capture_output=True, text=True)
    print('same list without the first step:', out2.stdout == out.stdout)
    big = [int(x) for x in open('sample64.txt')]
    pick = random.sample(big, 60) + [n for n in big if n > 2**63][:20] + random.sample(psp, 20) + [18446744073709551557, 18446744073709551533, 18446744073709551615, 3*5*7*11*13*17*19*23*29*31*37*41*43*47]
    wrong = 0; cmpd = 0
    for kind_r in [(17,0),(17,1),(19,0),(19,1),(25,0),(25,1),(32,0),(16,2)]:
        r, kind = kind_r
    rr = [17, 19, 25, 32]
    out3 = subprocess.run(['./ring3','-dump','-fit']+[str(r) for r in rr], input='\n'.join(map(str,pick))+'\n', capture_output=True, text=True)
    for l in out3.stdout.split('\n'):
        f = l.split()
        if len(f) < 4 or f[0] == 'ring': continue
        n, r, kind = int(f[0]), int(f[1]), int(f[2]); cells = [int(x) for x in f[3:]]
        ref = ringpow(n, r, kind)
        cmpd += 1
        if cells != ref: wrong += 1
    print('final rings compared with plain Python (rings 17, 19, 25, 32 and their first steps, numbers up to 2^64):', cmpd, '| different:', wrong)
