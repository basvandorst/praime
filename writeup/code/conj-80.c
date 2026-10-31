// Conjecture 80: running 5-adic balance W over all primes. Build: gcc -O3 -o conj80 conj-80.c -lm
// Run: ./conj80 7 8000000001   (about 12 s; prints the total and the minimum prefix from p = 7 on; W(3) = 2)
// Conjecture 80: for primes p in [lo, hi) (lo >= 7, odd), accumulate d(p) = v5(p^2+1) if p = 2,3 mod 5,
// -v5(p-1) if p = 1 mod 5, 0 if p = 4 mod 5. Prints the total, the minimum prefix sum (and where), the
// number of primes. Segmented sieve of Eratosthenes over odd numbers.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
static int v5_sq1(uint64_t p){ int k=0; unsigned __int128 m=5; while(1){ unsigned __int128 r=p%m; if((r*r+1)%m) break; k++; m*=5; } return k; }
static int v5(uint64_t x){ int k=0; while(x%5==0){ x/=5; k++; } return k; }
int main(int argc,char**argv){
  uint64_t lo=strtoull(argv[1],0,10), hi=strtoull(argv[2],0,10);
  uint64_t R=(uint64_t)sqrt((double)hi)+2; char *bs=calloc(R+1,1); uint64_t *bp=malloc(sizeof(uint64_t)*R); long nb=0;
  for(uint64_t i=2;i<=R;i++) if(!bs[i]){ if(i>2) bp[nb++]=i; for(uint64_t j=i*i;j<=R;j+=i) bs[j]=1; }
  const uint64_t S=1<<22;                 // odd numbers per segment
  unsigned char *seg=malloc(S);
  long long sum=0, minpre=0; uint64_t minat=0, nprimes=0;
  if(lo%2==0) lo++;
  for(uint64_t base=lo; base<hi; base+=2*S){          // segment covers odd numbers base, base+2, ..., base+2(S-1)
    memset(seg,0,S);
    uint64_t top=base+2*S;
    for(long t=0;t<nb;t++){ uint64_t p=bp[t]; if(p*p>=top) break;
      uint64_t st=((base+p-1)/p)*p; if(st<p*p) st=p*p; if(st%2==0) st+=p;   // first odd multiple >= base
      for(uint64_t j=st;j<top;j+=2*p) seg[(j-base)/2]=1; }
    for(uint64_t i=0;i<S;i++){ uint64_t p=base+2*i; if(p>=hi) break; if(seg[i] || p<2) continue;
      if(p<=R && bs[p]) continue;          // (only matters if segment overlaps base range; p composite)
      nprimes++;
      int r=p%5; if(r==1) sum-=v5(p-1); else if(r==2||r==3) sum+=v5_sq1(p);
      if(sum<minpre){ minpre=sum; minat=p; }
    }
  }
  printf("lo=%llu hi=%llu primes=%llu total=%lld minprefix=%lld at p=%llu\n",(unsigned long long)lo,(unsigned long long)hi,(unsigned long long)nprimes,sum,minpre,(unsigned long long)minat);
  return 0;
}
