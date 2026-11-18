#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
typedef uint64_t u64;
int main(int argc,char**argv){u64 lo=strtoull(argv[1],0,10),hi=strtoull(argv[2],0,10);u64 sq=(u64)sqrt((double)hi)+2;char*small=calloc(sq+1,1);u64 np=0;u64*pr=malloc(sizeof(u64)*(sq/2+10));
  for(u64 i=2;i<=sq;i++){if(!small[i]){pr[np++]=i;for(u64 j=i*i;j<=sq;j+=i)small[j]=1;}}
  const u64 BL=1<<20;char*seg=malloc(BL);if(lo%2==0)lo++;if(lo<3)lo=3;static char ob[1<<16];setvbuf(stdout,ob,_IOFBF,sizeof ob);
  for(u64 base=lo;base<hi;base+=2*BL){memset(seg,0,BL);u64 end=base+2*BL;if(end>hi)end=hi;
    for(u64 k=1;k<np;k++){u64 p=pr[k];if(p*p>=end)break;u64 st=p*p;if(st<base)st=((base+p-1)/p)*p;if(st%2==0)st+=p;for(u64 j=st;j<end;j+=2*p)seg[(j-base)>>1]=1;}
    for(u64 n=base;n<end;n+=2)if(seg[(n-base)>>1])printf("%llu\n",(unsigned long long)n);}
  return 0;}
