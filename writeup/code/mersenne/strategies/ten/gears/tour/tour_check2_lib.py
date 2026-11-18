"""The two-cell wheel in plain Python (shared by tour_check2.py and duo_row.py)."""
from sympy import jacobi_symbol, primerange
SP=list(primerange(3,200))
def step(n):
    for c in SP:
        if n==c: continue
        j=jacobi_symbol(c,n)
        if j==0: return -c
        if j==-1: return c
    return 0
def home(n):
    c=step(n)
    if c<=0: return False
    a,b=1,1
    for bit in bin(n)[3:]:
        a,b=(a*a+c*b*b)%n,(2*a*b)%n
        if bit=='1': a,b=(a+c*b)%n,(a+b)%n
    return a==1 and b==n-1
