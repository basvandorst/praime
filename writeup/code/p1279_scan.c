// Puzzle 1279. Build: gcc -O3 -march=native -o p1279_scan p1279_scan.c
// Run:   ./p1279_scan ../data/p1279_primes.txt 10080 SIGN K0 K1 REPORT   (SIGN 1 Sierpinski, -1 Riesel)
// Prints K lines (k with a minimal cover of >= REPORT primes), N lines (k with >= 2 minimal covers),
// Z lines (covered k without a minimal cover of <= 20 primes) and totals on stderr.
// Scan odd k in [K0, K1) for k*2^n + s (s = +1 Sierpinski, -1 Riesel) covered by primes of a list
// (p, m_p = ord_p(2), all m_p | L).  Stage 1: block sieve of a 64-point necessary test (n = 0..63).
// Stage 2: full test on Z/L.  Stage 3: all minimal covering subsets (size <= MAXC), by a canonical DFS.
// Output: "K k sizes..." for every covered k with a minimal cover of size >= REPORT, plus totals.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned __int128 u128;
#define MAXP 256
#define MAXL 20160
#define W ((MAXL+63)/64)
static int NP, L, LW, S, MAXC=20, REPORT=11; static long long nzero=0; static long long multi[16];
static uint64_t P[MAXP]; static int Mo[MAXP];
static uint32_t *res[MAXP];      // useful residues, sorted, with class index
static int *cls[MAXP];
static uint64_t pw(uint64_t b, uint64_t e, uint64_t m){uint64_t r=1;b%=m;while(e){if(e&1)r=(u128)r*b%m;b=(u128)b*b%m;e>>=1;}return r;}
static int findclass(int j, uint64_t r){ // binary search r in res[j]
  int lo=0, hi=Mo[j]-1; while(lo<=hi){int mid=(lo+hi)/2; if(res[j][mid]==r) return cls[j][mid]; if(res[j][mid]<r) lo=mid+1; else hi=mid-1;} return -1; }
static uint64_t full[W];
static int nused; static int up[MAXP], ua[MAXP]; static uint64_t umask[MAXP][W];
static int sizecount[64]; static uint64_t found[4096]; static int nfound; static int best_size;
static long long ncov_total=0, nsample=0;
static int popc_cover(uint64_t *m){ for(int w=0;w<LW;w++) if(m[w]!=full[w]) return 0; return 1; }
static int chosen[64], nch;
static void dfs(uint64_t *cov, uint64_t chmask){
  if(popc_cover(cov)){
    // minimality: each chosen prime has a point not covered by the others
    for(int c=0;c<nch;c++){ uint64_t o[W]; memset(o,0,sizeof(uint64_t)*LW);
      for(int d=0;d<nch;d++) if(d!=c) for(int w=0;w<LW;w++) o[w]|=umask[chosen[d]][w];
      if(popc_cover(o)) return; }
    for(int f=0;f<nfound;f++) if(found[f]==chmask) return;
    if(nfound<4096) found[nfound++]=chmask; return; }
  if(nch>=MAXC) return;
  // smallest uncovered point
  int x=-1; for(int w=0;w<LW;w++){ uint64_t u=~cov[w]&full[w]; if(u){ x=w*64+__builtin_ctzll(u); break; } }
  for(int i=0;i<nused;i++){ if(chmask>>i&1) continue; if(!(umask[i][x>>6]>>(x&63)&1)) continue;
    uint64_t nc[W]; for(int w=0;w<LW;w++) nc[w]=cov[w]|umask[i][w];
    chosen[nch++]=i; dfs(nc, chmask|(1ULL<<i)); nch--; }
}
static void analyse(uint64_t k){
  uint64_t cov[W]; memset(cov,0,sizeof cov); nused=0;
  for(int j=0;j<NP;j++){ int a=findclass(j, k%P[j]); if(a<0) continue;
    up[nused]=j; ua[nused]=a; memset(umask[nused],0,sizeof(uint64_t)*LW);
    for(int n=a;n<L;n+=Mo[j]) umask[nused][n>>6]|=1ULL<<(n&63);
    for(int w=0;w<LW;w++) cov[w]|=umask[nused][w];
    nused++; if(nused>=64) break; }
  if(!popc_cover(cov)) return;
  ncov_total++;
  uint64_t z[W]; memset(z,0,sizeof z); nfound=0; nch=0; dfs(z,0);
  int mx=0, has[64]={0};
  for(int f=0;f<nfound;f++){ int c=__builtin_popcountll(found[f]); has[c]=1; if(c>mx) mx=c; }
  for(int c=0;c<64;c++) if(has[c]) sizecount[c]++;
  if(nfound>=2){ printf("N %llu %d", (unsigned long long)k, nfound); if(nfound>=3) for(int f=0;f<nfound;f++){ printf(" |"); for(int i=0;i<nused;i++) if(found[f]>>i&1) printf(" %llu",(unsigned long long)P[up[i]]); } printf("\n"); }
  multi[nfound<16?nfound:15]++;
  if(nfound==0){ nzero++; printf("Z %llu covered, no minimal cover with <= %d primes\n",(unsigned long long)k,MAXC); }
  if(mx>=REPORT){ printf("K %llu sizes", (unsigned long long)k); for(int c=0;c<64;c++) if(has[c]) printf(" %d",c);
    for(int f=0;f<nfound;f++) if(__builtin_popcountll(found[f])>=REPORT){ printf(" |"); for(int i=0;i<nused;i++) if(found[f]>>i&1) printf(" %llu",(unsigned long long)P[up[i]]); }
    printf("\n"); fflush(stdout); }
}
int main(int argc, char **argv){
  // args: listfile L sign K0 K1 [REPORT]
  FILE *f=fopen(argv[1],"r"); L=atoi(argv[2]); S=atoi(argv[3]);
  uint64_t K0=strtoull(argv[4],0,10), K1=strtoull(argv[5],0,10); if(argc>6) REPORT=atoi(argv[6]);
  LW=(L+63)/64; memset(full,0,sizeof full); for(int n=0;n<L;n++) full[n>>6]|=1ULL<<(n&63);
  unsigned long long p; int m; NP=0;
  while(fscanf(f,"%llu %d",&p,&m)==2){ if(L%m) continue; P[NP]=p; Mo[NP]=m; NP++; }
  // useful residues: k*2^a + s == 0 mod p  <=>  k == -s * 2^{-a} mod p
  for(int j=0;j<NP;j++){ res[j]=malloc(4*Mo[j]); cls[j]=malloc(4*Mo[j]);
    uint64_t inv2=(P[j]+1)/2, t=1; // t = 2^{-a}
    for(int a=0;a<Mo[j];a++){ uint64_t r = (S>0)? (P[j]-t)%P[j] : t; res[j][a]=r; cls[j][a]=a; t=(u128)t*inv2%P[j]; }
    for(int a=1;a<Mo[j];a++){ uint32_t r=res[j][a]; int c=cls[j][a]; int b=a-1; while(b>=0 && res[j][b]>r){res[j][b+1]=res[j][b]; cls[j][b+1]=cls[j][b]; b--;} res[j][b+1]=r; cls[j][b+1]=c; } }
  // sample masks: point n in 0..63 covered by class a mod m  <=>  n % m == a
  const int BS=1<<20; uint64_t *arr=malloc(8*(size_t)BS);
  if(K0%2==0) K0++;
  // next index per (prime, class)
  int64_t **nxt=malloc(sizeof(int64_t*)*NP); uint64_t **smask=malloc(sizeof(uint64_t*)*NP);
  for(int j=0;j<NP;j++){ nxt[j]=malloc(8*Mo[j]); smask[j]=malloc(8*Mo[j]);
    uint64_t inv2=(P[j]+1)/2;
    for(int a=0;a<Mo[j];a++){ uint64_t r=res[j][a]; // k0 + 2i == r mod p  -> i == (r - k0) * inv2
      uint64_t d=(r + P[j] - K0%P[j])%P[j]; nxt[j][a]=(int64_t)((u128)d*inv2%P[j]);
      uint64_t sm=0; for(int n=cls[j][a]; n<64; n+=Mo[j]) sm|=1ULL<<n; smask[j][a]=sm; } }
  uint64_t nk=(K1-K0+1)/2; 
  for(uint64_t base=0; base<nk; base+=BS){
    int64_t len = (nk-base<BS)? (int64_t)(nk-base) : BS;
    memset(arr,0,8*(size_t)len);
    for(int j=0;j<NP;j++){ int64_t pj=P[j]; for(int a=0;a<Mo[j];a++){ int64_t i=nxt[j][a]; uint64_t sm=smask[j][a];
        if(sm) for(; i<len; i+=pj) arr[i]|=sm; else if(i<len) i+=((len-1-i)/pj+1)*pj;
        nxt[j][a]=i-len; } }
    for(int64_t i=0;i<len;i++) if(arr[i]==~0ULL){ nsample++; analyse(K0+2*(base+i)); }
  }
  fprintf(stderr,"range [%llu,%llu) sign %d L %d primes %d: sample-pass %lld covered %lld\n",(unsigned long long)K0,(unsigned long long)K1,S,L,NP,nsample,ncov_total);
  fprintf(stderr,"covered with no minimal cover of <= %d primes: %lld\n",MAXC,nzero);
  fprintf(stderr,"k by number of minimal covers:"); for(int c=0;c<16;c++) if(multi[c]) fprintf(stderr," %d:%lld",c,multi[c]); fprintf(stderr,"\n");
  fprintf(stderr,"k with a minimal cover of size c:"); for(int c=0;c<64;c++) if(sizecount[c]) fprintf(stderr," %d:%d",c,sizecount[c]); fprintf(stderr,"\n");
  return 0; }
