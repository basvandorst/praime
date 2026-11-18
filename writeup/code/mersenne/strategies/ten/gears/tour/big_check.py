"""64-bit check: numbers near 2^64 through tour.c and the Python reference must agree."""
import subprocess, random
from ref import ring_test, poly_test
random.seed(5)
nums=[random.randrange(2**63,2**64)|1 for _ in range(60)]+[18446744073709551557,18446744073709551533,2**64-1,(2**32-5)*(2**32-17),4294967291*4294967279]
open('big.txt','w').write('\n'.join(map(str,nums))+'\n')
for mach,fn in (( ['ring','5','1'], lambda n:ring_test(n,5)),(['ring','11','1'],lambda n:ring_test(n,11)),(['poly','3','1','0','1','1'],lambda n:poly_test(n,[1,1,0])),(['exp','2'],lambda n:pow(2,n,n)==2%n)):
    out=subprocess.run(['./tour','stdin']+mach,stdin=open('big.txt'),capture_output=True,text=True).stdout.split()
    mine=[n for n in nums if fn(n)]
    print(mach,'C passes:',len(out),'Python passes:',len(mine),'equal:',sorted(map(int,out))==sorted(mine))
