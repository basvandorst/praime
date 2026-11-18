"""Independent check of the tournament in plain Python (integers and lists only; ref.py holds the ring test).
Recomputes from scratch, without the C programs:
 a. the numbers that fool the doubling wheel below 1e7, and their count below 1e8
 b. Pascal's wheel on the smallest misfit ring: all candidates below 1e10, a sample of the list up to 2^64, a sample of primes
 c. each ring of 5, 7, 9, 11, 13, 16 cells on the candidates below 1e10 that do not fit it
 d. the doubling wheel with the Fibonacci wheel for n = 2, 3 (mod 5), below 1e10
 e. first mistakes and counts of the single wheels below 1e6
 f. the two lemmas on the data (sum of the cells; a ring folds onto its divisors)
Reads: psp2_to1e10.txt (our own list), sample64.txt (every 4,000th line of the published list), primes_sample.txt,
       ring_psp2_to1e10.json and the single_*.txt logs for the comparison."""
import json, sys, time
from math import gcd, isqrt
from ref import ring_test
t0=time.time()
def misfit(n):
    for r in range(2,400):
        if gcd(n,r)==1 and (n*n)%r!=1: return r
    return 0
def fibmod(n,m):
    a,b=0,1
    for bit in bin(n)[2:]:
        c=a*((2*b-a)%m)%m; d=(a*a+b*b)%m
        a,b=(d,(c+d)%m) if bit=='1' else (c,d)
    return a
def perrin_home(n):
    a,b,c=3,0,2
    M=[[0,1,0],[0,0,1],[1,1,0]]
    def mul(A,B): return [[sum(A[i][k]*B[k][j] for k in range(3))%n for j in range(3)] for i in range(3)]
    R=[[1,0,0],[0,1,0],[0,0,1]]; e=n; A=M
    while e:
        if e&1: R=mul(R,A)
        A=mul(A,A); e>>=1
    return (R[0][0]*3+R[0][1]*0+R[0][2]*2)%n==0
N7=10**7; sieve=bytearray([1])*(N7+1); sieve[0]=sieve[1]=0
for i in range(2,isqrt(N7)+1):
    if sieve[i]: sieve[i*i::i]=bytearray(len(sieve[i*i::i]))
mine=[n for n in range(9,N7,2) if not sieve[n] and pow(2,n,n)==2]
L=[int(x) for x in open('psp2_to1e10.txt')]
print("a. odd composites below 1e7 that fool the doubling wheel: %d; equal to the start of the C list: %s"%(len(mine),mine==[n for n in L if n<N7]))
print("   C list: %d below 1e8, %d below 1e10"%(sum(1 for n in L if n<10**8),len(L)))
bad=[n for n in L if (lambda r:r and ring_test(n,r))(misfit(n))]
print("b. Pascal's wheel on the smallest misfit ring, all %d candidates below 1e10: came home: %d"%(len(L),len(bad)),flush=True)
S=[int(x) for x in open('sample64.txt')]
use={}
bad=[]
for n in S:
    r=misfit(n); use[r]=use.get(r,0)+1
    if r==0 or ring_test(n,r): bad.append(n)
print("   sample of the published list (every 4,000th line): %d numbers from %d to %d; came home: %d; rings used: %s"%(len(S),S[0],S[-1],len(bad),sorted(use.items())),flush=True)
P=[int(x) for x in open('primes_sample.txt')]; P=P[::57]
print("   primes of 8 to 64 bits (every 57th of the sample): %d; came home: %d"%(len(P),sum(1 for p in P if ring_test(p,misfit(p)))),flush=True)
for r in (5,7,9,11,13,16):
    c=[n for n in L if gcd(n,r)==1 and (n*n)%r!=1]
    print("c. ring of %2d cells: %d candidates below 1e10 do not fit; came home: %d"%(r,len(c),sum(1 for n in c if ring_test(n,r))),flush=True)
c=[n for n in L if n%5 in (2,3)]
print("d. doubling and Fibonacci wheel, n = 2 or 3 (mod 5), below 1e10: %d candidates; Fibonacci wheel home (n | F(n+1)): %d"%(len(c),sum(1 for n in c if fibmod(n+1,n)==0)))
N6=10**6
comp=[n for n in range(9,N6,2) if not sieve[n]]
fib=[n for n in comp if n%5 and fibmod(n-(1 if n%5 in (1,4) else -1),n)==0]
per=[n for n in comp if perrin_home(n)]
def clist(fn): return [int(x) for x in open(fn).read().split('\n') if x.strip().isdigit()]
print("e. below 1e6: Fibonacci wheel %d mistakes, first %d; Perrin's wheel %d mistakes %s; doubling wheel %d"%(len(fib),fib[0],len(per),per,sum(1 for n in mine if n<N6)))
try:
    cp=clist('single_perrin_1e8.txt'); cr=clist('single_ring5_1e8.txt')
    print("   C logs agree: Perrin below 1e6 %s; ring of 5 cells among ALL odd composites below 1e8 = the %d found on the list: %s"%([n for n in cp if n<N6]==per,len(cr),cr==[n for n in json.load(open('ring_psp2_to1e10.json'))['5'] if n<10**8]))
except FileNotFoundError: print("   (C logs not present)")
d=json.load(open('ring_psp2_to1e10.json')); Ls=set(L)
print("f. every ring mistake below 1e10 fools the doubling wheel: %s; a mistake on ring r is a mistake on every ring d | r (d >= 2): %s"%(all(n in Ls for v in d.values() for n in v),all(set(d[str(r)])<=set(d[str(q)]) for r in range(2,41) for q in range(2,r) if r%q==0)))
print("   ring of 5 cells, first mistake by plain Python among all odd composites below 3e5 (no shortcut):",[n for n in range(9,300000,2) if not sieve[n] and ring_test(n,5)][:2])
print("time %.0f s"%(time.time()-t0))
