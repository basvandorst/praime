// Puzzle 1104, Q1: primes p whose sums S_n of n-th powers of primitive roots are prime for n = 1..5.
// Only p with p-1 squarefree can work (else p | S_1). This program lists p <= LIMIT with S_1 and S_2 prime;
// S_3..S_5 are tested in Python. Usage: q1 LIMIT
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m){ return (u128)a*b % m; }
static uint64_t powmod(uint64_t a, uint64_t e, uint64_t m){ uint64_t r=1; a%=m; while(e){ if(e&1) r=mulmod(r,a,m); a=mulmod(a,a,m); e>>=1;} return r; }
static int isprime64(uint64_t n){
  if(n<2) return 0; static const uint64_t sp[]={2,3,5,7,11,13,17,19,23,29,31,37};
  for(int i=0;i<12;i++){ if(n%sp[i]==0) return n==sp[i]; }
  uint64_t d=n-1; int s=0; while(!(d&1)){d>>=1;s++;}
  for(int i=0;i<12;i++){ uint64_t x=powmod(sp[i],d,n); if(x==1||x==n-1) continue; int ok=0;
    for(int r=1;r<s;r++){ x=mulmod(x,x,n); if(x==n-1){ok=1;break;} } if(!ok) return 0; }
  return 1;
}
int main(int argc,char**argv){
  uint64_t LIMIT=strtoull(argv[1],0,10);
  char *comp=calloc(LIMIT+1,1); for(uint64_t i=2;i*i<=LIMIT;i++) if(!comp[i]) for(uint64_t j=i*i;j<=LIMIT;j+=i) comp[j]=1;
  char *mark=malloc(LIMIT+1);
  long tested=0, sqf=0, s1=0, s2=0;
  for(uint64_t p=3;p<=LIMIT;p++){ if(comp[p]) continue; tested++;
    uint64_t m=p-1, t=m; int nf=0; uint64_t f[20]; int squarefree=1;
    for(uint64_t q=2;q*q<=t;q++) if(t%q==0){ f[nf++]=q; int e=0; while(t%q==0){t/=q;e++;} if(e>1) squarefree=0; }
    if(t>1) f[nf++]=t;
    if(!squarefree) continue; sqf++;
    // find a generator
    uint64_t g=2; for(;;g++){ int ok=1; for(int i=0;i<nf;i++) if(powmod(g,m/f[i],p)==1){ok=0;break;} if(ok) break; }
    // mark exponents i with gcd(i, m) > 1
    memset(mark,0,m); for(int i=0;i<nf;i++) for(uint64_t j=0;j<m;j+=f[i]) mark[j]=1;
    uint64_t S1=0; u128 S2=0; uint64_t x=1;
    for(uint64_t i=0;i<m;i++){ if(!mark[i]){ S1+=x; S2+=(u128)x*x; } x=mulmod(x,g,p); }
    if(!isprime64(S1)) continue; s1++;
    if(S2>>63) { fprintf(stderr,"S2 overflow at p=%llu\n",(unsigned long long)p); return 1; }
    if(!isprime64((uint64_t)S2)) continue; s2++;
    printf("%llu %llu %llu\n",(unsigned long long)p,(unsigned long long)S1,(unsigned long long)S2);
  }
  fprintf(stderr,"primes %ld, p-1 squarefree %ld, S1 prime %ld, S1 and S2 prime %ld\n",tested,sqf,s1,s2);
  return 0;
}
