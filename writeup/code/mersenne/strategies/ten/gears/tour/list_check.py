"""Compare our own list of odd composites with 2^n = 2 (mod n) (psp2.c: segmented sieve and a Fermat test)
with the start of the published list of Feitsma and Galway (psps-below-2-to-64.txt.bz2).
usage: list_check.py PATH_TO_BZ2"""
import struct, subprocess, sys, hashlib
def rd(fn):
    d=open(fn,'rb').read(); n=len(d)//8; return list(struct.unpack('<%dQ'%n,d[:8*n]))
a=rd('psp2_1e8.bin'); b=rd('psp2_1e10.bin'); c=rd('psp2_1e11.bin')
own=a+b+c
p=subprocess.Popen('bzcat %s | head -%d'%(sys.argv[1],len(own)),shell=True,stdout=subprocess.PIPE,text=True)
pub=[int(x) for x in p.stdout.read().split()]
print("our list: %d below 1e8, %d below 1e10, and %d more up to %d (the run to 1e11 was stopped there)"%(len(a),len(a)+len(b),len(c),c[-1]))
print("identical to the first %d lines of the published list: %s"%(len(own),own==pub))
h=hashlib.sha256()
with open(sys.argv[1],'rb') as f:
    for blk in iter(lambda:f.read(1<<22),b''): h.update(blk)
print("published file: https://www.cecm.sfu.ca/Pseudoprimes/psps-below-2-to-64.txt.bz2, sha256",h.hexdigest())
