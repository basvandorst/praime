// Puzzle 1279, second implementation (no sieve). Build: gcc -O3 -march=native -o p1279_scan_direct p1279_scan_direct.c
// Run: ./p1279_scan_direct ../data/p1279_primes.txt 10080 SIGN K0 K1
// Independent check: for each odd k in [K0,K1) compute k mod p directly for every prime of the list,
// look the residue up in a hash table of useful residues, mark the covered n in a byte array over Z/L,
// and if all of Z/L is covered enumerate the minimal covers by include/exclude recursion.
// Prints "C k count sizes s1 s2 ..." (number of minimal covers, their distinct sizes) for every covered k.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned __int128 u128;
static int NP, L, S; static uint64_t P[256]; static int Mo[256];
static uint64_t *hk[256]; static int *hv[256]; static int hsz[256];
static int lookup(int j, uint64_t r){ int h=(int)((r*0x9E3779B97F4A7C15ULL)>>40)%hsz[j]; while(hv[j][h]>=0){ if(hk[j][h]==r) return hv[j][h]; h=(h+1)%hsz[j]; } return -1; }
static int nu; static int cp[256], cm[256], ca[256]; static unsigned char cov[30240];
static int cnt[30240]; static int sizes[300]; static int sel[256];
static int covered_by_sel(int skip){ // is Z/L covered by selected primes except index skip
  static unsigned char t[30240]; memset(t,0,L);
  for(int i=0;i<nu;i++) if(sel[i] && i!=skip) for(int n=ca[i];n<L;n+=cm[i]) t[n]=1;
  for(int n=0;n<L;n++) if(!t[n]) return 0; return 1; }
static int pc[30240]; static int ncovers;  // how many selected primes cover each point
static void rec(int i, int nsel){
  if(i==nu){ if(nsel==0 || !covered_by_sel(-1)) return;
    for(int c=0;c<nu;c++) if(sel[c] && covered_by_sel(c)) return;
    sizes[nsel]=1; ncovers++; return; }
  if(nsel>20) return;
  // include i only if it adds a point that no selected prime covers (else every superset is non-minimal)
  int adds=0; for(int n=ca[i];n<L;n+=cm[i]) if(!pc[n]){adds=1;break;}
  if(adds){ sel[i]=1; for(int n=ca[i];n<L;n+=cm[i]) pc[n]++; rec(i+1,nsel+1);
            for(int n=ca[i];n<L;n+=cm[i]) pc[n]--; sel[i]=0; }
  // exclude i only if the selected primes and the primes after i can still cover Z/L
  unsigned char *t=malloc(L); memset(t,0,L);
  for(int c=0;c<nu;c++) if((c<i && sel[c]) || c>i) for(int n=ca[c];n<L;n+=cm[c]) t[n]=1;
  int ok=1; for(int n=0;n<L;n++) if(!t[n]){ok=0;break;} free(t);
  if(ok) rec(i+1,nsel);
}
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); L=atoi(argv[2]); S=atoi(argv[3]); uint64_t K0=strtoull(argv[4],0,10), K1=strtoull(argv[5],0,10);
  unsigned long long p; int m; NP=0; while(fscanf(f,"%llu %d",&p,&m)==2){ if(L%m) continue; P[NP]=p; Mo[NP]=m; NP++; }
  for(int j=0;j<NP;j++){ hsz[j]=4*Mo[j]+7; hk[j]=malloc(8*hsz[j]); hv[j]=malloc(4*hsz[j]); for(int h=0;h<hsz[j];h++) hv[j][h]=-1;
    // useful residue for class a: k*2^a + S == 0  ->  k == -S * (2^a)^{-1}; compute 2^a then invert by Fermat
    for(int a=0;a<Mo[j];a++){ uint64_t t=1,b=2,e=a; while(e){ if(e&1) t=(u128)t*b%P[j]; b=(u128)b*b%P[j]; e>>=1; }
      uint64_t inv=1; b=t; e=P[j]-2; while(e){ if(e&1) inv=(u128)inv*b%P[j]; b=(u128)b*b%P[j]; e>>=1; }
      uint64_t r = S>0 ? (P[j]-inv)%P[j] : inv;
      int h=(int)((r*0x9E3779B97F4A7C15ULL)>>40)%hsz[j]; while(hv[j][h]>=0) h=(h+1)%hsz[j]; hk[j][h]=r; hv[j][h]=a; } }
  long long ncov=0;
  for(uint64_t k=K0|1; k<K1; k+=2){
    nu=0; double w=0;
    for(int j=0;j<NP;j++){ int a=lookup(j,k%P[j]); if(a<0) continue; cp[nu]=j; cm[nu]=Mo[j]; ca[nu]=a; w+=1.0/Mo[j]; nu++; }
    if(w<1.0-1e-12) continue;                       // density bound: a cover needs sum 1/m >= 1
    memset(cov,0,L); for(int i=0;i<nu;i++) for(int n=ca[i];n<L;n+=cm[i]) cov[n]=1;
    int ok=1; for(int n=0;n<L;n++) if(!cov[n]){ok=0;break;} if(!ok) continue;
    ncov++; memset(sizes,0,sizeof sizes); memset(sel,0,sizeof sel); memset(pc,0,sizeof pc); ncovers=0; rec(0,0);
    printf("C %llu %d sizes",(unsigned long long)k,ncovers); for(int c=0;c<64;c++) if(sizes[c]) printf(" %d",c); printf("\n");
  }
  fprintf(stderr,"scan2 [%llu,%llu) sign %d: covered %lld\n",(unsigned long long)K0,(unsigned long long)K1,S,ncov);
}
