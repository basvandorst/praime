// Conjecture 115 scanner. Build: gcc -O3 -o conj115 conj-115.c   Run: ./conj115 5 30000  (prints zeros "Z" and failures "FAIL")
// Faster scan: G(m) mod p for all m (using G(m) = -sum_{j>m} j^m mod p when m > (p-1)/2, m != p-1),
// then G(m) mod p^2 only where p | G(m). Prime p fails iff some k <= (p-1)/2 has v(G(k)) + v(G(p-k)) >= 2.
// Valid for p < 65536 (p^2 < 2^32).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static int isprime(uint64_t n){ if(n<2) return 0; for(uint64_t d=2; d*d<=n; d++) if(n%d==0) return 0; return 1; }
static uint64_t pw_mod(uint64_t b, uint64_t e, uint64_t m){ uint64_t r=1%m; b%=m; while(e){ if(e&1) r=r*b%m; b=b*b%m; e>>=1;} return r; }
static int vcap2(uint64_t p, uint64_t m){ // v_p(G(m)) capped at 2, assuming p | G(m)
  uint64_t p2=p*p, s=0; for(uint64_t j=1;j<=m;j++){ s+=pw_mod(j,m,p2); if(s>=p2) s-=p2; } return s==0?2:1; }
int main(int argc,char**argv){
  uint64_t lo=atoll(argv[1]), hi=atoll(argv[2]);
  uint32_t *pw=malloc(4*70000); unsigned char *v=malloc(70000);
  for(uint64_t p=(lo<5?5:lo); p<hi; p++){
    if(!isprime(p)) continue;
    uint64_t h=(p-1)/2;
    for(uint64_t m=1;m<p;m++) v[m]=0;
    // m = 1..h : pw[j] = j^m for j<=m
    for(uint64_t m=1;m<=h;m++){
      uint64_t s=0;
      for(uint64_t j=1;j<m;j++){ pw[j]=(uint32_t)((uint64_t)pw[j]*j%p); s+=pw[j]; }
      pw[m]=(uint32_t)pw_mod(m,m,p); s+=pw[m];
      if(s%p==0) v[m]=1;
    }
    // m = h+1..p-2 : pw[j] = j^m for j in (m, p-1]
    uint64_t m0=h+1;
    for(uint64_t j=m0+1;j<p;j++) pw[j]=(uint32_t)pw_mod(j,m0,p);
    for(uint64_t m=m0;m<=p-2;m++){
      if(m>m0) for(uint64_t j=m+1;j<p;j++) pw[j]=(uint32_t)((uint64_t)pw[j]*j%p);
      uint64_t s=0; for(uint64_t j=m+1;j<p;j++) s+=pw[j];
      if(s%p==0) v[m]=1;
    }
    // m = p-1: G(p-1) = sum j^(p-1) = p-1 mod p, never 0.
    for(uint64_t m=1;m<p-1;m++) if(v[m]){ v[m]=vcap2(p,m); printf("Z p=%lu m=%lu v>=%d\n",(unsigned long)p,(unsigned long)m,v[m]); }
    for(uint64_t k=1;k<=h;k++) if(v[k]+v[p-k]>=2) printf("FAIL p=%lu k=%lu v(G(k))=%d v(G(p-k))=%d\n",(unsigned long)p,(unsigned long)k,v[k],v[p-k]);
    fflush(stdout);
  }
  return 0;
}
