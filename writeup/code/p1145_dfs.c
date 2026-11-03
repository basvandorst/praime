#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
typedef unsigned __int128 u128;
typedef unsigned long long u64;
#define MAXP 64
#define MAXW 40
int N, W, np;
u64 P[MAXP]; int O[MAXP]; u64 *S[MAXP];
u64 *MASK[MAXP]; /* MASK[i] + r*W */
double LP[MAXP], LB;
int cap[MAXP];
u64 full[MAXW];
int used[MAXP], res[MAXP], order[MAXP];
int *forb[MAXP]; int stamp=1;
long long ncov=0, nnode=0; int stopflag=0; int exact_prune=0; int mrv=0; int nshift=1; int knap=0; int splitK=1, splitPart=0, splitDepth=2; long long splitCtr=0; int fix5=-1, fix7=-1; int capdyn[MAXP][MAXP]; double tlimit; clock_t t0; double w0; double wall(){struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts); return ts.tv_sec+ts.tv_nsec*1e-9-w0;}
u128 bestall=0; double report_below;
u128 bestset_a; 
static inline u64 mulmod(u64 a,u64 b,u64 m){return (u64)((u128)a*b%m);}
u64 inv64(u64 a, u64 m){ /* a,m coprime */
  __int128 t=0,nt=1; __int128 r=m,nr=a%m;
  while(nr){__int128 q=r/nr; __int128 tmp=t-q*nt;t=nt;nt=tmp; tmp=r-q*nr;r=nr;nr=tmp;}
  if(t<0)t+=m; return (u64)t;}
void pr128(u128 x,char*buf){char tmp[64];int n=0;if(!x){strcpy(buf,"0");return;}while(x){tmp[n++]='0'+(int)(x%10);x/=10;}for(int i=0;i<n;i++)buf[i]=tmp[n-1-i];buf[n]=0;}
int sel[MAXP], nsel;
#define HT (1<<16)
u64 hkey[HT]; long long hcnt[HT]; u128 hbest[HT]; int hused[HT];
void hrec(u64 key,u128 a){ u64 h=(key*0x9E3779B97F4A7C15ULL)>>48; while(hused[h]&&hkey[h]!=key) h=(h+1)&(HT-1); if(!hused[h]){hused[h]=1;hkey[h]=key;hcnt[h]=0;hbest[h]=a;} hcnt[h]++; if(a<hbest[h])hbest[h]=a; }
void crtmin1(int t){
  u128 M=6; for(int i=0;i<nsel;i++) M*=P[sel[i]];
  u128 D[MAXP+1]; u128 a=0;
  for(int i=0;i<nsel;i++){int j=sel[i]; u128 Mi=M/P[j]; u64 mi=(u64)(Mi%P[j]); u64 e=inv64(mi,P[j]); u64 s=S[j][(res[j]+t)%O[j]];
    u128 Ei=Mi*e%M; /* Mi*e may overflow? Mi<M<2^128, e<p ... use mulmod128 */
    /* safe mul: */ { u128 r=0,b=Mi%M; u64 ee=e; while(ee){ if(ee&1){ r+=b; if(r>=M) r-=M;} b<<=1; if(b>=M) b-=M; ee>>=1;} Ei=r; }
    u128 t=0,b=Ei; u64 ss=s; while(ss){ if(ss&1){t+=b; if(t>=M)t-=M;} b<<=1; if(b>=M)b-=M; ss>>=1;}
    a+=t; if(a>=M)a-=M; D[i]=t+t; if(D[i]>=M)D[i]-=M; }
  { u128 M6=M/6; u64 e=inv64((u64)(M6%6),6); u128 t=M6*e%M; a+=t; if(a>=M)a-=M; D[nsel]=t+t; if(D[nsel]>=M)D[nsel]-=M; }
  int k=nsel+1; u128 best=a; u64 g=0;
  for(u64 c=1;c<(1ULL<<k);c++){ int b=__builtin_ctzll(c); u64 ng=g^(1ULL<<b); 
    if(ng&(1ULL<<b)){ /* flip + to - : subtract D */ if(a>=D[b]) a-=D[b]; else a=a+M-D[b]; }
    else { a+=D[b]; if(a>=M)a-=M; }
    g=ng; if(a<best)best=a; }
  { u64 key=0; for(int i=0;i<nsel;i++) key|=1ULL<<sel[i]; hrec(key,best); }
  if(best<bestall){bestall=best;}
  if((double)best<report_below){char buf[64];pr128(best,buf); printf("COVER a=%s set",buf); for(int i=0;i<nsel;i++)printf(" %llu:%d",P[sel[i]],(res[sel[i]]+t)%O[sel[i]]); printf("\n"); fflush(stdout);}
}
void crtmin(){ for(int t=0;t<nshift;t++) crtmin1(t); }
void dfs(u64*cov,double lp,int unc){
  nnode++; if(splitK>1 && nsel==splitDepth){ long long c=splitCtr++; if(c%splitK!=splitPart) return; } if((nnode&0x3ffffff)==0){fprintf(stderr,"nodes=%lld covers=%lld depth=%d t=%.1f\n",nnode,ncov,nsel,wall());}
  if(stopflag) return; if((nnode&0xffff)==0 && wall()>tlimit){stopflag=1; return;}
  if(unc==0){ ncov++; crtmin(); return; }
  /* capacity prune */
  int capsum=0; int lvl = nsel<exact_prune? nsel : exact_prune-1; int *cd = capdyn[lvl<0?0:lvl];
  if(nsel<exact_prune){ for(int i=0;i<np;i++){ cd[i]=0; if(used[i]) continue; int b=0; u64*mk=MASK[i];
      for(int r=0;r<O[i];r++,mk+=W){ int c=0; for(int w=0;w<W;w++) c+=__builtin_popcountll(mk[w]&~cov[w]); if(c>b)b=c; }
      cd[i]=b; } }
  if(knap){ int idx[MAXP]; double ratio[MAXP]; int n=0; double rem=LB-lp;
    for(int i=0;i<np;i++) if(!used[i]&&LP[i]<=rem){ int b=exact_prune>0?cd[i]:cap[i]; if(b<=0) continue; idx[n]=i; ratio[n]=b/LP[i]; n++; }
    for(int x=1;x<n;x++){ int ix=idx[x]; double rx=ratio[x]; int y=x-1; while(y>=0&&ratio[y]<rx){idx[y+1]=idx[y];ratio[y+1]=ratio[y];y--;} idx[y+1]=ix; ratio[y+1]=rx; }
    double tot=0;
    for(int x=0;x<n;x++){ int i=idx[x]; int b=exact_prune>0?cd[i]:cap[i]; if(LP[i]<=rem){ tot+=b; rem-=LP[i]; } else { tot+=b*rem/LP[i]; break; } }
    if(tot<unc) return; capsum=unc; }
  else for(int i=0;i<np;i++) if(!used[i]&&lp+LP[i]<=LB) capsum+= exact_prune>0? cd[i] : cap[i];
  if(capsum<unc) return;
  int m=-1;
  if(mrv){ int bestc=1<<30; int avail[MAXP], na=0; for(int i=0;i<np;i++) if(!used[i]&&lp+LP[i]<=LB) avail[na++]=i;
    for(int w=0;w<W;w++){ u64 x=~cov[w]&full[w]; while(x){ int mm=w*64+__builtin_ctzll(x); x&=x-1; int c=0;
        for(int k=0;k<na;k++){ int i=avail[k]; if(!forb[i][mm%O[i]]) c++; }
        if(c<bestc){bestc=c;m=mm; if(c==0) return;} } } }
  else for(int w=0;w<W;w++){ u64 x=~cov[w]&full[w]; if(x){ m=w*64+__builtin_ctzll(x); break;} }
  u64 nc[MAXW]; int mystamp=++stamp;
  for(int i=0;i<np;i++){ if(used[i]||lp+LP[i]>LB) continue; int r=m%O[i]; if(forb[i][r]) continue; if((i==fix5||i==fix7)&&r!=0) continue;
    int nu=unc; u64*mk=MASK[i]+(size_t)r*W; for(int w=0;w<W;w++){ nc[w]=cov[w]|mk[w]; nu-=__builtin_popcountll(nc[w]^cov[w]); }
    used[i]=1; res[i]=r; sel[nsel++]=i; dfs(nc,lp+LP[i],nu); nsel--; used[i]=0;
    forb[i][r]=mystamp; }
  for(int i=0;i<np;i++){ int r=m%O[i]; if(forb[i][r]==mystamp) forb[i][r]=0; }
}
int main(int argc,char**argv){
  tlimit=atof(argv[1]); report_below=atof(argv[2]); exact_prune=argc>3?atoi(argv[3]):0; mrv=argc>4?atoi(argv[4]):0; int symm=argc>5?atoi(argv[5]):0; if(argc>7){splitK=atoi(argv[6]); splitPart=atoi(argv[7]);} if(argc>8) splitDepth=atoi(argv[8]); if(argc>9) knap=atoi(argv[9]);
  if(scanf("%d %lf %d",&N,&LB,&np)!=3) return 1; W=(N+63)/64;
  for(int i=0;i<N;i++) full[i/64]|=1ULL<<(i%64);
  for(int i=0;i<np;i++){ scanf("%llu %d",&P[i],&O[i]); LP[i]=log((double)P[i]); cap[i]=N/O[i];
    S[i]=malloc(sizeof(u64)*O[i]); for(int r=0;r<O[i];r++) scanf("%llu",&S[i][r]);
    MASK[i]=calloc((size_t)O[i]*W,sizeof(u64)); for(int r=0;r<O[i];r++) for(int m=r;m<N;m+=O[i]) MASK[i][(size_t)r*W+m/64]|=1ULL<<(m%64);
    forb[i]=calloc(O[i],sizeof(int)); }
  if(symm){ for(int i=0;i<np;i++){ if(P[i]==5&&O[i]==2) fix5=i; if(P[i]==7&&O[i]==3) fix7=i; } nshift=(fix5>=0?2:1)*(fix7>=0?3:1); }
  bestall=(u128)1<<126; bestset_a=bestall; t0=clock(); w0=0; w0=wall();
  u64 cov[MAXW]; memset(cov,0,sizeof cov); nsel=0;
  dfs(cov,log(6.0)*0,N);
  char buf[64]; pr128(bestall,buf);
  double el=wall();
  for(int h=0;h<HT;h++) if(hused[h]){ char b2[64]; pr128(hbest[h],b2); double lp=0; printf("SET covers=%lld best_a=%s logP=",hcnt[h],b2); for(int i=0;i<np;i++) if(hkey[h]>>i&1) lp+=LP[i]; printf("%.2f primes",lp/log(10)); for(int i=0;i<np;i++) if(hkey[h]>>i&1) printf(" %llu",P[i]); printf("\n"); }
  printf("DONE N=%d covers=%lld nodes=%lld best_a=%s time=%.1f %s\n",N,ncov,nnode,buf,el,el>tlimit?"TIMEOUT":"complete");
  return 0; }
