/* Search F(y) = A y^2 + B y + C (A >= 1, 0 < B < A, A + B even, C odd) for long runs of consecutive y in [-Y, Y]
   with |F(y)| prime. Filter (reasoning): for every odd prime p <= PS, either p | A and p | B (then F = C mod p, need
   p !| C), or the discriminant D = B^2 - 4AC is a non-residue mod p (then p never divides F). Runs >= MINRUN printed.
   usage: qs AMIN AMAX CMAX PS Y MINRUN */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned __int128 u128;
static uint8_t *bits; static uint64_t LIM;
static uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m){ return (u128)a*b % m; }
static uint64_t powmod(uint64_t a, uint64_t e, uint64_t m){ uint64_t r=1; a%=m; while(e){ if(e&1) r=mulmod(r,a,m); a=mulmod(a,a,m); e>>=1;} return r; }
static int mr(uint64_t n){ if(n<2) return 0; static const uint64_t P[]={2,3,5,7,11,13,17,19,23,29,31,37};
  for(int i=0;i<12;i++){ if(n%P[i]==0) return n==P[i]; }
  uint64_t d=n-1; int s=0; while(!(d&1)){d>>=1;s++;}
  for(int i=0;i<12;i++){ uint64_t x=powmod(P[i],d,n); if(x==1||x==n-1) continue; int ok=0; for(int r=1;r<s;r++){ x=mulmod(x,x,n); if(x==n-1){ok=1;break;} } if(!ok) return 0; }
  return 1; }
static inline int isp(int64_t v){ uint64_t n = v<0 ? (uint64_t)(-v) : (uint64_t)v; if(n<LIM){ if(n==2) return 1; if(!(n&1)) return 0; return (bits[n>>4]>>((n>>1)&7))&1; } return mr(n); }
int main(int argc, char **argv){
  int64_t AMIN=atoll(argv[1]), AMAX=atoll(argv[2]), CMAX=atoll(argv[3]); int PS=atoi(argv[4]), Y=atoi(argv[5]), MINRUN=atoi(argv[6]);
  LIM = 1ULL<<29; bits = malloc(LIM>>4); memset(bits,0xff,LIM>>4); bits[0]&=~2; /* 1 not prime */
  for(uint64_t i=3;i*i<LIM;i+=2) if((bits[i>>4]>>((i>>1)&7))&1) for(uint64_t j=i*i;j<LIM;j+=2*i) bits[j>>4]&=~(1<<((j>>1)&7));
  int pr[64], np=0; for(int p=3;p<=PS;p+=2){ int q=1; for(int d=3;d*d<=p;d+=2) if(p%d==0) q=0; if(q) pr[np++]=p; }
  static uint8_t nonres[64][128];
  for(int i=0;i<np;i++){ int p=pr[i]; memset(nonres[i],1,128); for(int x=0;x<p;x++) nonres[i][(x*x)%p]=0; }
  int W = 2*Y+1; int64_t *vals = malloc(sizeof(int64_t)*W); uint8_t *ok = malloc(W);
  long long tested=0;
  for(int64_t A=AMIN; A<=AMAX; A++) for(int64_t B=1; B<A; B++){
    if((A+B)&1) continue;
    /* per prime: allowed C mod p */
    static uint8_t allow[64][128]; int bad=0;
    for(int i=0;i<np;i++){ int p=pr[i];
      if(A%p==0){ if(B%p){ bad=1; break; } for(int c=0;c<p;c++) allow[i][c]= (c!=0); }
      else for(int c=0;c<p;c++){ int64_t D = ((B*B - 4*A*c) % p + p) % p; allow[i][c]=nonres[i][D]; } }
    if(bad) continue;
    for(int64_t C=-CMAX|1; C<=CMAX; C+=2){
      int pass=1; for(int i=0;i<np;i++){ int p=pr[i]; if(!allow[i][((C%p)+p)%p]){ pass=0; break; } }
      if(!pass) continue; tested++;
      /* checkpoints at multiples of h = ceil(MINRUN/2) */
      int h=(MINRUN+1)/2; int any=0;
      for(int y=-(Y/h)*h; y+h<=Y; y+=h){ if(isp(A*y*y+B*y+C) && isp(A*(y+h)*(y+h)+B*(y+h)+C)){ any=1; break; } }
      if(!any) continue;
      for(int y=-Y;y<=Y;y++){ ok[y+Y]=isp(A*y*y+B*y+C); }
      int best=0,bs=0,cur=0; for(int k=0;k<W;k++){ cur= ok[k]? cur+1:0; if(cur>best){best=cur; bs=k-cur+1-Y;} }
      if(best>=MINRUN) printf("%lld %lld %lld %d %d\n",(long long)A,(long long)B,(long long)C,bs,best);
    }
    fflush(stdout);
  }
  fprintf(stderr,"tested %lld\n",tested); return 0; }
