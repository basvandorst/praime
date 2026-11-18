/* Second, independent implementation of the ring test (x+1)^n = x^n + 1 mod (n, x^r - 1).
   Plain 128-bit remainders (no Montgomery), right-to-left binary exponentiation, full products.
   usage: ring2 R [all] < numbers   ; tests the numbers with gcd(n,R)=1 and n^2 != 1 mod R (or all coprime ones with "all") */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef unsigned __int128 u128;
static u64 n; static int r;
static void mul(const u64*a,const u64*b,u64*out){
  u64 t[128]; for(int k=0;k<r;k++)t[k]=0;
  for(int i=0;i<r;i++){ if(!a[i])continue;
    for(int j=0;j<r;j++){ if(!b[j])continue; int k=(i+j)%r;
      u64 p=(u64)(((u128)a[i]*b[j])%n); u128 s=(u128)t[k]+p; t[k]=(u64)(s%n);} }
  memcpy(out,t,sizeof(u64)*r); }
static u64 gcd(u64 a,u64 b){while(b){u64 t=a%b;a=b;b=t;}return a;}
int main(int argc,char**argv){ r=atoi(argv[1]); int all=argc>2; char buf[64]; u64 tested=0,skipped=0,liars=0,last=0;
  while(fgets(buf,64,stdin)){ n=strtoull(buf,0,10); if(n<3||n%2==0)continue; last=n;
    u64 m=n%r; if(gcd(m,r)!=1){skipped++;continue;} if(!all&&(m*m)%r==1){skipped++;continue;}
    u64 res[128],b[128]; for(int k=0;k<r;k++){res[k]=0;b[k]=0;} res[0]=1; b[0]=1; b[1]=1; u64 e=n;
    while(e){ if(e&1)mul(res,b,res); e>>=1; if(e)mul(b,b,b); }
    int ok=1; for(int k=0;k<r;k++){ u64 w=(k==0?1:0)+(k==(int)m?1:0); if(res[k]!=w%n){ok=0;break;} }
    tested++; if(ok){liars++;printf("%llu\n",n);fflush(stdout);}
    if((tested&0x3FFFFF)==0)fprintf(stderr,"  .. tested %llu at %llu liars %llu\n",tested,n,liars); }
  fprintf(stderr,"ring %d: tested %llu skipped %llu liars %llu last %llu\n",r,tested,skipped,liars,last); return 0; }
