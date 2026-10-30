// Puzzle 868: 3x3 magic squares of L-digit palindromic primes.
// Entries: d0*(10^(L-1)+1) + sum_j v_j*w_j, j=0 center (10^h), j>=1: 10^(h-j)+10^(h+j), h=(L-1)/2, m free half-positions.
// Square = c + alpha*a + beta*b (alpha,beta in {-1,0,1}) digitwise without carries.
// usage: pal L m d0 nprimes
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef struct { int a; int s; void *d; } mpz_s; typedef mpz_s mpz_t[1];
void __gmpz_init(mpz_t); void __gmpz_set(mpz_t, const mpz_t); void __gmpz_set_ui(mpz_t, unsigned long);
void __gmpz_ui_pow_ui(mpz_t, unsigned long, unsigned long); void __gmpz_add(mpz_t, const mpz_t, const mpz_t);
void __gmpz_addmul_ui(mpz_t, const mpz_t, unsigned long); void __gmpz_mul_ui(mpz_t, const mpz_t, unsigned long);
unsigned long __gmpz_fdiv_ui(const mpz_t, unsigned long); int __gmpz_probab_prime_p(const mpz_t, int);
typedef uint64_t u64;
static int L,m,d0; static u64 M; static u64 *S;
static inline int getS(int64_t i){ return (S[i>>6]>>(i&63))&1; }
static int dig[12]; static u64 P10[12];
static inline u64 bits(int64_t p,int w){ u64 wd=p>>6; int off=p&63; u64 v=S[wd]>>off; if(off) v|=S[wd+1]<<(64-off); return v&((1ULL<<w)-1); }
static inline u64 rev(u64 v,int w){ u64 r=0; for(int i=0;i<w;i++) if(v>>i&1) r|=1ULL<<(w-1-i); return r; }
int main(int argc,char**argv){ L=atoi(argv[1]); m=atoi(argv[2]); d0=atoi(argv[3]); int NQ=atoi(argv[4]); int h=(L-1)/2;
  P10[0]=1; for(int i=1;i<12;i++) P10[i]=P10[i-1]*10; M=P10[m];
  mpz_t base,w[12],t,N; __gmpz_init(base); __gmpz_init(t); __gmpz_init(N);
  __gmpz_ui_pow_ui(base,10,L-1); __gmpz_set_ui(t,1); __gmpz_add(base,base,t); __gmpz_mul_ui(base,base,d0);
  for(int j=0;j<m;j++){ __gmpz_init(w[j]); __gmpz_ui_pow_ui(w[j],10,h-j); if(j){ __gmpz_ui_pow_ui(t,10,h+j); __gmpz_add(w[j],w[j],t);} }
  // primes for the sieve
  int *Q=malloc(sizeof(int)*NQ); int nq=0; for(int q=3;nq<NQ;q+=2){int ok=1;for(int d=3;d*d<=q;d+=2)if(q%d==0){ok=0;break;} if(ok&&q!=5)Q[nq++]=q;}
  int ml=m/2, mh=m-ml; u64 NL=P10[ml], NH=P10[mh];
  // lo residues bucketed per q
  uint32_t **ord=malloc(sizeof(void*)*nq), **st=malloc(sizeof(void*)*nq); unsigned long *rb=malloc(sizeof(long)*nq), (*wr)[12]=malloc(sizeof(*wr)*nq);
  for(int k=0;k<nq;k++){ int q=Q[k]; rb[k]=__gmpz_fdiv_ui(base,q); for(int j=0;j<m;j++) wr[k][j]=__gmpz_fdiv_ui(w[j],q);
    uint32_t *res=malloc(4*NL); st[k]=calloc(q+1,4); ord[k]=malloc(4*NL);
    for(u64 lo=0;lo<NL;lo++){ unsigned long r=0; u64 x=lo; for(int j=0;j<ml;j++){ r+=(x%10)*wr[k][j]; x/=10;} res[lo]=r%q; st[k][res[lo]+1]++; }
    for(int r=0;r<q;r++) st[k][r+1]+=st[k][r]; uint32_t *pos=malloc(4*(q+1)); memcpy(pos,st[k],4*(q+1));
    for(u64 lo=0;lo<NL;lo++) ord[k][pos[res[lo]]++]=lo; free(pos); free(res); }
  S=calloc(M/64+2,8); unsigned char *mark=malloc(NL); u64 tested=0, nprime=0;
  for(u64 hi=0;hi<NH;hi++){ memset(mark,1,NL);
    for(int k=0;k<nq;k++){ int q=Q[k]; unsigned long r=rb[k]; u64 x=hi; for(int j=ml;j<m;j++){ r+=(x%10)*wr[k][j]; x/=10;} r%=q;
      unsigned long need=(q-r)%q; for(uint32_t i=st[k][need];i<st[k][need+1];i++) mark[ord[k][i]]=0; }
    for(u64 lo=0;lo<NL;lo++) if(mark[lo]){ u64 idx=hi*NL+lo; __gmpz_set(N,base); u64 x=idx; for(int j=0;j<m;j++){ unsigned long dj=x%10; x/=10; if(dj) __gmpz_addmul_ui(N,w[j],dj);} tested++;
      if(__gmpz_probab_prime_p(N,1)){ S[idx>>6]|=1ULL<<(idx&63); nprime++; } } }
  fprintf(stderr,"L=%d m=%d d0=%d: tested %llu, PRP %llu (density %.5f)\n",L,m,d0,(unsigned long long)tested,(unsigned long long)nprime,(double)nprime/M);
  // phase 2
  long squares=0; int maxD=1<<16; int64_t *DA=malloc(8*maxD); signed char (*Dv)[12]=malloc(12*maxD);
  for(u64 c=0;c<M;c++){ if(!getS(c)) continue; int cd[12], tt[12]; u64 x=c; for(int j=0;j<m;j++){cd[j]=x%10;x/=10; tt[j]=cd[j]<9-cd[j]?cd[j]:9-cd[j];}
    int nD=0; int a[12]={0}; for(int j=1;j<m;j++) a[j]=-tt[j]; int w0=2*tt[0]+1;
    while(1){ int64_t O=0; for(int j=1;j<m;j++) O+=a[j]*(int64_t)P10[j];
      u64 x=bits((int64_t)c+O-tt[0],w0), y=bits((int64_t)c-O-tt[0],w0);   // y bit k <-> c-O-t0+k  = c-(O+a0) with a0=t0-k
      u64 msk=x & rev(y,w0);
      while(msk){ int k=__builtin_ctzll(msk); msk&=msk-1; int a0=k-tt[0]; int64_t A=O+a0; if(A<=0||nD>=maxD) continue;
        DA[nD]=A; Dv[nD][0]=a0; for(int j=1;j<m;j++) Dv[nD][j]=a[j]; nD++; }
      int j=1; while(j<m){ if(a[j]<tt[j]){a[j]++;break;} a[j]=-tt[j]; j++; } if(j==m) break; }
    for(int i=0;i<nD;i++) for(int k=i+1;k<nD;k++){ int ok=1; for(int j=0;j<m;j++){ int s=abs(Dv[i][j])+abs(Dv[k][j]); if(s>tt[j]){ok=0;break;} } if(!ok) continue;
      int64_t A=DA[i],B=DA[k]; if(A==B||A==-B||A==2*B||A==-2*B||2*A==B||2*A==-B) continue;
      if(getS(c+A+B)&&getS(c-A-B)&&getS(c+A-B)&&getS(c-A+B)){ squares++; printf("SQ c=%llu A=%lld B=%lld\n",(unsigned long long)c,(long long)A,(long long)B); fflush(stdout);} } }
  fprintf(stderr,"squares found: %ld\n",squares); return 0; }
