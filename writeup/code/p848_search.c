// Puzzle 848 Q1, k >= 5: one of p, q is 7.
// mode A: N_n = 7*10^(d+n) + x  (p = 7, q = x with d digits)     mode B: N_n = x*10^(n+1) + 7  (p = x, q = 7)
// usage: s848 mode KS xlo xhi   -> prints x with x prime and N_0..N_(KS-1) PRP; reports full run length (up to 16)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef struct { int a; int s; void *d; } mpz_s; typedef mpz_s mpz_t[1];
void __gmpz_init(mpz_t); void __gmpz_import(mpz_t, size_t, int, size_t, int, size_t, const void*); int __gmpz_probab_prime_p(const mpz_t, int);
typedef unsigned __int128 u128; typedef uint64_t u64;
static mpz_t Zt; static int prp(u128 n){ u64 w[2]={(u64)n,(u64)(n>>64)}; __gmpz_import(Zt,2,-1,8,0,0,w); return __gmpz_probab_prime_p(Zt,2)>0; }
static u64 powm(u64 a,u64 e,u64 m){u64 r=1;a%=m;while(e){if(e&1)r=r*a%m;a=a*a%m;e>>=1;}return r;}
#define NP 400
int main(int c,char**v){ char mode=v[1][0]; int KS=atoi(v[2]); u64 xlo=strtoull(v[3],0,10), xhi=strtoull(v[4],0,10); __gmpz_init(Zt);
  int pr[NP],np=0; for(int q=3;np<NP;q+=2){int ok=1;for(int d=3;d*d<=q;d+=2)if(q%d==0){ok=0;break;} if(ok&&q!=5&&q!=7)pr[np++]=q;}
  u128 p10[40]; p10[0]=1; for(int i=1;i<40;i++) p10[i]=p10[i-1]*10;
  const u64 BL=1<<20; unsigned char*sv=malloc(BL); static int forb[NP][20]; static int nf[NP]; u64 tested=0;
  for(u64 lo=xlo; lo<xhi; lo+=BL){ u64 n=BL; if(lo+n>xhi) n=xhi-lo;
    int d=0; { u64 t=lo; while(t){d++;t/=10;} } u64 dlim=(u64)p10[d]; if(mode=='A' && lo+n>dlim) n=dlim-lo; // keep one digit length per block
    for(int j=0;j<np;j++){ int q=pr[j]; nf[j]=0; forb[j][nf[j]++]=0;          // x itself prime
      for(int i=0;i<KS;i++){ int r; if(mode=='A') r=(int)((q-(7*powm(10,d+i,q))%q)%q); else r=(int)((q-7*powm(powm(10,i+1,q),q-2,q)%q)%q);
        int dup=0; for(int k=0;k<nf[j];k++) if(forb[j][k]==r) dup=1; if(!dup) forb[j][nf[j]++]=r; } }
    memset(sv,1,n);
    for(int j=0;j<np;j++){ u64 q=pr[j]; for(int k=0;k<nf[j];k++){ u64 st=(forb[j][k]+q-lo%q)%q; for(u64 x=st;x<n;x+=q) sv[x]=0; } }
    for(u64 i=0;i<n;i++){ u64 x=lo+i; if(!sv[i]||!(x&1)||x%5==0) continue; tested++;
      if(!prp(x)) continue; int run=0;
      for(int t=0;t<16;t++){ u128 N = (mode=='A') ? 7*p10[d+t]+x : (u128)x*p10[t+1]+7; if(!prp(N)) break; run++; }
      if(run>=KS){ printf("%c x=%llu run=%d (k=%d) p+q=%llu\n",mode,(unsigned long long)x,run,run-1,(unsigned long long)(x+7)); fflush(stdout);} }
    lo += n - BL; }  // adjust when block was shortened
  fprintf(stderr,"%c KS=%d [%llu,%llu) tested=%llu\n",mode,KS,(unsigned long long)xlo,(unsigned long long)xhi,(unsigned long long)tested); }
