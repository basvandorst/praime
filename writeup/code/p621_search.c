// Puzzle 621. Build: gcc -O2 -o p621_search p621_search.c -lm
// Run: ./p621_search C X PRIME   (PRIME = 1: entries prime, offsets multiples of 6; PRIME = 0: any centre, step 1)
//      ./p621_search 30000 2000 0  finds Resta's square (c = 25818, a = 381, b = 1); ./p621_search 20000000 9000 1 finds none
// 3x3 magic squares M of primes (centre c, offsets +-a, +-b, +-(a+b), +-(a-b)) such that p(M) is magic, p(k) = k-th prime.
// p(M) magic <=> p(c+x) + p(c-x) = 2 p(c) for x in {a, b, a+b, a-b} and d(a+b) = d(a) + d(b), d(a-b) = d(a) - d(b),
// where d(x) = p(c+x) - p(c).  Entries prime forces c prime and a, b = 0 (mod 6).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
int main(int argc,char**argv){ long C=atol(argv[1]), X=atol(argv[2]); int PR=atoi(argv[3]); long STEP=PR?6:1; long K=C+X+10;
  // sieve to get the first K primes
  long N=(long)(K*(__builtin_log(K)+__builtin_log(__builtin_log(K)))*1.1)+100;
  char *sv=calloc(N+1,1); uint32_t *p=malloc(sizeof(uint32_t)*(K+1)); long np=1;
  for(long i=2;i<=N && np<=K;i++){ if(!sv[i]){ p[np++]=i; if((long)i*i<=N) for(long j=(long)i*i;j<=N;j+=i) sv[j]=1; } }
  if(np<=K){ fprintf(stderr,"sieve too small\n"); return 1; }
  // isprime for indices up to K: sv covers that range (K < N)
  long found=0; static long S[100000];
  for(long c=5;c<=C;c++){ if(PR && (sv[c]||c%2==0)) continue;
    long ns=0; for(long x=STEP;x<=X && x<c;x+=STEP) if((long)p[c+x]+(long)p[c-x]==2L*p[c]) S[ns++]=x;
    if(ns<4) continue;
    for(long i=0;i<ns;i++) for(long j=0;j<i;j++){ long a=S[i], b=S[j]; // a > b
      long s=a+b, d=a-b; if(s>X||s>=c) continue;
      // s and d must be in S
      int hs=0, hd=0; for(long t=0;t<ns;t++){ if(S[t]==s) hs=1; if(S[t]==d) hd=1; } if(!hs||!hd) continue;
      long da=p[c+a]-p[c], db=p[c+b]-p[c], ds=p[c+s]-p[c], dd=p[c+d]-p[c];
      if(ds!=da+db || dd!=da-db) continue;
      long e[9]={c, c+a, c-a, c+b, c-b, c+s, c-s, c+d, c-d}; int allp=1; for(int t=0;t<9;t++) if(sv[e[t]]){allp=0;break;}
      if(PR && !allp) continue;
      printf("c=%ld a=%ld b=%ld  magic sum %ld; p(M) magic sum %ld\n", c, a, b, 3*c, 3L*p[c]); found++; fflush(stdout);
      if(found>=5) return 0; } }
  fprintf(stderr,"searched centres <= %ld, offsets <= %ld: %ld found\n",C,X,found); return 0; }
