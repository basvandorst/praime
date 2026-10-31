// all n in [lo, hi) with 1 + tau(n)*sigma(n) == n*rad(n); also report n | 1 + tau*sigma
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef unsigned __int128 u128; typedef uint64_t u64;
int main(int c,char**v){ u64 lo=strtoull(v[1],0,10), hi=strtoull(v[2],0,10); const u64 S=1<<20;
  u64 P=(u64)sqrt((double)hi)+2; char*cmp=calloc(P+1,1); u64*pr=malloc(sizeof(u64)*P); u64 np=0;
  for(u64 i=2;i<=P;i++) if(!cmp[i]){pr[np++]=i; for(u64 k=i*i;k<=P;k+=i)cmp[k]=1;}
  u64 *rem=malloc(8*S), *sg=malloc(8*S), *rd=malloc(8*S); uint32_t *tu=malloc(4*S); u64 hits=0;
  for(u64 b=lo;b<hi;b+=S){ u64 n0=b, m=(hi-b<S)?hi-b:S;
    for(u64 i=0;i<m;i++){rem[i]=n0+i; sg[i]=1; rd[i]=1; tu[i]=1;}
    for(u64 j=0;j<np;j++){ u64 p=pr[j]; if(p*p>n0+m) break; u64 st=(n0+p-1)/p*p;
      for(u64 x=st;x<n0+m;x+=p){ u64 i=x-n0; u64 r=rem[i]; int e=0; u64 pp=1, s=1; while(r%p==0){r/=p;e++;pp*=p;s+=pp;} rem[i]=r; tu[i]*=e+1; sg[i]*=s; rd[i]*=p; } }
    for(u64 i=0;i<m;i++){ u64 n=n0+i; if(n<2)continue; if(rem[i]>1){tu[i]*=2; sg[i]*=rem[i]+1; rd[i]*=rem[i];}
      u128 L=(u128)tu[i]*sg[i]+1, R=(u128)n*rd[i];
      if(L==R){printf("SOL %llu\n",(unsigned long long)n);fflush(stdout);}
      else if(L%n==0){printf("div %llu\n",(unsigned long long)n);fflush(stdout);} } }
  return 0; }
