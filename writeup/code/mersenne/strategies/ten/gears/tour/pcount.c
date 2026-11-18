/* pcount LO HI: the number of odd primes in [LO, HI), by a segmented sieve of Eratosthenes (no connection with the wheels). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
typedef uint64_t u64;
int main(int argc,char**argv){u64 lo=strtoull(argv[1],0,10),hi=strtoull(argv[2],0,10);if(lo<3)lo=3;if(lo%2==0)lo++;
  u64 r=(u64)sqrtl((long double)hi)+2;char*sm=calloc(r+1,1);u64*pr=malloc(sizeof(u64)*(r/2+10));int np=0;
  for(u64 i=3;i<=r;i+=2){if(!sm[i]){pr[np++]=i;for(u64 j=i*i;j<=r;j+=2*i)sm[j]=1;}}
  const u64 SEG=1<<22;char*sv=malloc(SEG);u64 cnt=0;
  for(u64 a=lo;a<hi;a+=2*SEG){u64 b=a+2*SEG;if(b>hi)b=hi;u64 len=(b-a+1)/2;memset(sv,0,len);
    for(int k=0;k<np;k++){u64 p=pr[k];if(p*p>=b)break;u64 s=p*p;if(s<a){s=((a+p-1)/p)*p;if(s%2==0)s+=p;}for(u64 j=(s-a)/2;j<len;j+=p)sv[j]=1;}
    for(u64 j=0;j<len;j++)cnt+=!sv[j];}
  printf("odd primes in [%llu, %llu): %llu\n",(unsigned long long)lo,(unsigned long long)hi,(unsigned long long)cnt);return 0;}
