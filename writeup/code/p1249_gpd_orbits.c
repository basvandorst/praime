// Puzzles 1249/1250. Build: gcc -O2 -o p1249_gpd_orbits p1249_gpd_orbits.c
// Run: ./p1249_gpd_orbits 2 3 200000000 5 13 29 61 43 89 181 73 149
//      ./p1249_gpd_orbits 3 2 200000000 5 17 53 23 71 43 131 79 239 719 127 383 1151 691 83 251 151 13 41 167 503 1511 907 389
// For A(p) = GPD(a p + b): every prime 5 <= p < X (or 3 <= p for 1250) reaches a value < p or a known cycle.
// Processing primes in increasing order, this proves that every prime < X ends in a known cycle.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
static uint32_t *pr; static int np;
static uint64_t gpd(uint64_t n){ uint64_t g=1; for(int i=0;i<np;i++){ uint64_t q=pr[i]; if(q*q>n) break; if(n%q==0){ g=q; while(n%q==0) n/=q; } } if(n>1) g=n; return g; }
int main(int argc,char**argv){
  uint64_t a=atoll(argv[1]), b=atoll(argv[2]), X=atoll(argv[3]); int ncyc=argc-4; uint64_t cyc[64]; for(int i=0;i<ncyc;i++) cyc[i]=atoll(argv[4+i]);
  uint64_t B=2000000; char *s=calloc(B,1); np=0; pr=malloc(4*B);
  for(uint64_t i=2;i<B;i++){ if(!s[i]){ pr[np++]=i; for(uint64_t j=i*i;j<B;j+=i) s[j]=1; } }
  char *isp=calloc(X,1); for(uint64_t i=2;i<X;i++) isp[i]=1; for(uint64_t i=2;i*i<X;i++) if(isp[i]) for(uint64_t j=i*i;j<X;j+=i) isp[j]=0;
  uint64_t maxval=0, maxstart=0; long long steps=0; int pmin = (a==2)?5:3;
  for(uint64_t p=pmin;p<X;p++){ if(!isp[p]) continue;
    uint64_t x=p; int st=0;
    while(1){ int known=0; for(int i=0;i<ncyc;i++) if(x==cyc[i]) known=1; if(known) break; if(x<p && x>=(uint64_t)pmin) break;
      if(x<(uint64_t)pmin){ printf("orbit of %llu reaches %llu < %d\n",(unsigned long long)p,(unsigned long long)x,pmin); return 1; }
      if(a*x+b >= 4000000000000ULL){ printf("value too large for trial division\n"); return 1; }
      x=gpd(a*x+b); st++; if(x>maxval){maxval=x;maxstart=p;} if(st>100000){ printf("long orbit at %llu\n",(unsigned long long)p); return 1; } }
    steps+=st; }
  printf("a=%llu b=%llu: every prime %d <= p < %llu ends in a known cycle; %lld steps; largest value %llu (from %llu)\n",(unsigned long long)a,(unsigned long long)b,pmin,(unsigned long long)X,steps,(unsigned long long)maxval,(unsigned long long)maxstart);
  return 0; }
