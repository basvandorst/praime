// List every prime power p^a <= X (p = 2: any a >= 1; p odd: a even >= 2) whose sigma(p^a) = q^j, q odd prime.
// Output: q j p a
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef unsigned __int128 u128; typedef uint64_t u64;
static u64 mulm(u64 a,u64 b,u64 m){return (u128)a*b%m;}
static u64 powm(u64 a,u64 e,u64 m){u64 r=1;a%=m;while(e){if(e&1)r=mulm(r,a,m);a=mulm(a,a,m);e>>=1;}return r;}
static int mr(u64 n){ if(n<2)return 0; static const u64 sp[]={2,3,5,7,11,13,17,19,23,29,31,37};
  for(int i=0;i<12;i++){if(n%sp[i]==0)return n==sp[i];}
  u64 d=n-1;int s=0;while(!(d&1)){d>>=1;s++;}
  static const u64 B[]={2,325,9375,28178,450775,9780504,1795265022};
  for(int i=0;i<7;i++){u64 a=B[i]%n; if(!a)continue; u64 x=powm(a,d,n); if(x==1||x==n-1)continue; int ok=0;
    for(int r=1;r<s;r++){x=mulm(x,x,n); if(x==n-1){ok=1;break;}} if(!ok)return 0;} return 1;}
static u64 iroot(u64 s,int k){ u64 r=(u64)llround(pow((double)s,1.0/k)); if(r<1)r=1;
  for(u64 c=(r>2?r-2:1);c<=r+2;c++){u128 t=1;int over=0;for(int i=0;i<k;i++){t*=c;if(t>s){over=1;break;}} if(!over&&t==s)return c;} return 0;}
// q^j? returns q, sets *j
static u64 ppow(u64 s,int*j){ if(s%2==0) return 0; if(mr(s)){*j=1;return s;}
  static const int L[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61};
  for(int i=0;i<18;i++){u64 r=iroot(s,L[i]); if(r>1){int jj; u64 q=ppow(r,&jj); if(q){*j=jj*L[i];return q;}}} return 0;}
int main(int c,char**v){ double Xd=atof(v[1]); u64 X=(u64)Xd; u64 P=(u64)sqrtl((long double)X);
  // p = 2
  for(int a=1;a<=62 && (1ULL<<a)<=X;a++){int j; u64 q=ppow((1ULL<<(a+1))-1,&j); if(q) printf("%llu %d 2 %d\n",(unsigned long long)q,j,a);}
  // odd p <= P, segmented sieve
  u64 r=(u64)sqrt((double)P)+2; char*small=calloc(r+1,1); u64*pr=malloc(sizeof(u64)*r); int np=0;
  for(u64 i=2;i<=r;i++) if(!small[i]){pr[np++]=i; for(u64 k=i*i;k<=r;k+=i) small[k]=1;}
  const u64 SEG=1<<24; char*seg=malloc(SEG); u64 cnt=0;
  for(u64 lo=3;lo<=P;lo+=SEG){ u64 hi=lo+SEG-1; if(hi>P)hi=P; memset(seg,0,SEG);
    for(int i=0;i<np;i++){u64 p=pr[i]; if(p*p>hi)break; u64 st=(lo+p-1)/p*p; if(st<p*p)st=p*p; for(u64 k=st;k<=hi;k+=p) seg[k-lo]=1;}
    for(u64 n=lo;n<=hi;n++){ if(seg[n-lo]||!(n&1)) continue; u64 p=n; cnt++;
      // sigma(p^a) for even a>=2 with p^a<=X
      u128 pa=(u128)p*p, s=1+(u128)p+pa; int a=2;
      while(pa<=X){ int j; u64 q=ppow((u64)s,&j); if(q) printf("%llu %d %llu %d\n",(unsigned long long)q,j,(unsigned long long)p,a);
        pa*= (u128)p*p; s += (pa/p) + pa; a+=2; if(s>(u128)UINT64_MAX)break; } } }
  fprintf(stderr,"odd primes <= %llu: %llu\n",(unsigned long long)P,(unsigned long long)cnt); }
