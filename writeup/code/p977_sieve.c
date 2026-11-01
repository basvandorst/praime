/* Puzzle 977 sieve: q = 3k*2^n - 1, r = k*2^(n+2) - 1, s = 3k*2^(n+1) - 1 (= 2q + 1), odd k <= K.
   Remove k if a prime p <= P divides q, r or s. Prints surviving k. Usage: ./sieve n K P */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned long long u64; typedef unsigned __int128 u128;
static u64 pw(u64 b, u64 e, u64 m){ u64 r=1%m; b%=m; while(e){ if(e&1) r=(u128)r*b%m; b=(u128)b*b%m; e>>=1;} return r; }
static u64 inv(u64 a, u64 m){ return pw(a, m-2, m); }
int main(int argc, char**argv){
  u64 n = atoll(argv[1]), K = atoll(argv[2]), P = atoll(argv[3]);
  u64 nk = (K+1)/2;                      /* k = 2i+1, i in [0, nk) */
  uint8_t *bits = calloc(nk/8+1, 1);
  char *comp = calloc(P+1, 1);
  for(u64 i=2;i*i<=P;i++) if(!comp[i]) for(u64 j=i*i;j<=P;j+=i) comp[j]=1;
  for(u64 p=3;p<=P;p++){ if(comp[p]) continue;
    u64 t = pw(2, n, p), cls[3]; int nc = 0;
    /* q: 3k t = 1 ; r: 4k t = 1 ; s: 6k t = 1 */
    if(p != 3){ cls[nc++] = inv(3*t % p, p); cls[nc++] = inv(6*t % p, p); }
    cls[nc++] = inv(4*t % p, p);
    for(int c=0;c<nc;c++){
      u64 k0 = cls[c];                   /* k = k0 mod p; need k odd: k = k0 or k0+p */
      if(k0 % 2 == 0) k0 += p;
      for(u64 k=k0; k<=K; k+=2*p){ u64 i=(k-1)/2; bits[i>>3] |= (uint8_t)(1u<<(i&7)); }
    }
  }
  u64 cnt=0;
  for(u64 i=0;i<nk;i++) if(!(bits[i>>3] & (1u<<(i&7)))){ printf("%llu\n", 2*i+1); cnt++; }
  fprintf(stderr, "n=%llu K=%llu P=%llu survivors %llu\n", n, K, P, cnt);
  return 0;
}
