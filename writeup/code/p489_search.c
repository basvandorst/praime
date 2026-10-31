// all n <= N with sigma(sigma(phi(n))) == n ; spf for odd numbers stored as uint16 (0 = prime)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
static uint16_t *spf; static uint64_t N;
static inline uint32_t sp(uint64_t m){ /* m odd > 1 */ uint16_t s=spf[m>>1]; return s? s : 0; }
// factor m (<= N) ; returns k, primes in P, exps in E
static int fac(uint64_t m, uint64_t *P, int *E){
  int k=0; if(m%2==0){int e=0; while(m%2==0){m/=2;e++;} P[k]=2;E[k++]=e;}
  while(m>1){ uint64_t p=sp(m); if(!p) p=m; int e=0; while(m%p==0){m/=p;e++;} P[k]=p;E[k++]=e; }
  return k;}
static uint64_t sigma(uint64_t m){ uint64_t P[40]; int E[40]; int k=fac(m,P,E); uint64_t s=1;
  for(int i=0;i<k;i++){uint64_t t=1,pp=1; for(int j=0;j<E[i];j++){pp*=P[i];t+=pp;} s*=t;} return s;}
int main(int c,char**v){ N=strtoull(v[1],0,10); uint64_t H=N/2+1;
  spf=calloc(H,2); uint64_t r=(uint64_t)sqrt((double)N)+1;
  for(uint64_t p=3;p<=r;p+=2) if(!spf[p>>1]) for(uint64_t j=p*p;j<=N;j+=2*p) if(!spf[j>>1]) spf[j>>1]=p;
  uint64_t cnt=0;
  for(uint64_t n=1;n<=N;n++){ uint64_t P[40]; int E[40]; int k=fac(n,P,E); uint64_t ph=1;
    for(int i=0;i<k;i++){ ph*=P[i]-1; for(int j=1;j<E[i];j++) ph*=P[i]; }
    uint64_t M=sigma(ph); if(M>=n && !(n==1)) continue;   // sigma(M)>M except M=1
    if(sigma(M)==n){ printf("%llu\n",(unsigned long long)n); fflush(stdout);} cnt++; }
  fprintf(stderr,"N=%llu checked, %llu passed M<n\n",(unsigned long long)N,(unsigned long long)cnt); }
