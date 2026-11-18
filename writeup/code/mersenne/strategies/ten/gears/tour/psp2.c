#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef unsigned __int128 u128; typedef uint64_t u64;
static inline u64 mont_inv(u64 n){u64 x=n;for(int i=0;i<6;i++)x*=2-n*x;return x;}
static inline u64 redc(u128 t,u64 n,u64 ni){u64 m=(u64)t*ni;u128 mn=(u128)m*n;u64 r=(u64)((t-mn)>>64);return (t<mn)?r+n:r;}
int main(int argc,char**argv){
  u64 lo=strtoull(argv[1],0,10),hi=strtoull(argv[2],0,10);FILE*out=fopen(argv[3],"wb");
  u64 sq=(u64)sqrt((double)hi)+2;char*small=calloc(sq+1,1);u64 np=0;u64*pr=malloc(sizeof(u64)*(sq/2+10));
  for(u64 i=2;i<=sq;i++){if(!small[i]){pr[np++]=i;for(u64 j=i*i;j<=sq;j+=i)small[j]=1;}}
  const u64 BL=1<<21;char*seg=malloc(BL);u64 cnt=0;if(lo%2==0)lo++;if(lo<3)lo=3;
  for(u64 base=lo;base<hi;base+=2*BL){
    memset(seg,0,BL);u64 end=base+2*BL;if(end>hi)end=hi;
    for(u64 k=1;k<np;k++){u64 p=pr[k];if(p*p>=end)break;u64 st=p*p;if(st<base){st=((base+p-1)/p)*p;}if(st%2==0)st+=p;for(u64 j=st;j<end;j+=2*p)seg[(j-base)>>1]=1;}
    for(u64 n=base;n<end;n+=2){if(!seg[(n-base)>>1])continue;
      u64 ni=mont_inv(n);u64 one=(u64)(((u128)1<<64)%n);u64 two=one+one;if(two>=n)two-=n;
      u64 e=n-1;int top=63-__builtin_clzll(e);u64 x=two;
      for(int b=top-1;b>=0;b--){x=redc((u128)x*x,n,ni);if((e>>b)&1){x+=x;if(x>=n)x-=n;}}
      if(x==one){fwrite(&n,8,1,out);cnt++;}}
  }
  fclose(out);printf("%llu %llu %llu\n",(unsigned long long)lo,(unsigned long long)hi,(unsigned long long)cnt);return 0;}
