/* Sextics F = N/K with N(y) = y^6 + b4 y^4 + b3 y^3 + b2 y^2 + b1 y + b0 (b5 = 0 by shift), K = 72 (or 36).
   K = 72: N(y) mod 432 in {72, 360} for all y (F odd, 3 does not divide F). K = 36: N(y) mod 216 in {36, 180}.
   Shapes (b4, b3, b2) from random root sets; for each shape every b1 with |N| <= K*V on |y| <= H possible, and every
   b0 in the size interval by residue classes; filter: N has no root mod every prime 5..PS. usage:
   sx K NSHAPES SEED V H PS Y MINRUN RSPREAD */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <complex.h>
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
static inline int64_t md(int64_t x, int64_t p){ x%=p; return x<0?x+p:x; }
static uint64_t rs=88172645463325252ULL; static double urand(){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return (rs>>11)*(1.0/9007199254740992.0); }
int main(int argc, char **argv){
  int64_t K=atoll(argv[1]); long NS=atol(argv[2]); rs+=atoll(argv[3]); int64_t V=atoll(argv[4]); int H=atoi(argv[5]); int PS=atoi(argv[6]);
  int Y=atoi(argv[7]), MINRUN=atoi(argv[8]); double RS=atof(argv[9]);
  int64_t MOD = K==72 ? 432 : 216, T1 = K, T2 = K==72 ? 360 : 180;   /* N mod MOD must be T1 or T2 */
  int64_t B5 = K==72 ? 3 : 0;
  /* tables: for (b4,b3,b2,b1) mod m (m = 16 or 8 for the 2-part, 27 for the 3-part): is there b0 with N(y) mod m in the target set for all y */
  int m2 = K==72 ? 16 : 8; int t2 = (int)(T1 % m2); static uint8_t ok2[16][16][16][16], ok3[27][27][27][27]; static uint8_t any2[16][16][16], any3[27][27][27];
  for(int a4=0;a4<m2;a4++) for(int a3=0;a3<m2;a3++) for(int a2=0;a2<m2;a2++) for(int a1=0;a1<m2;a1++){ int ok=1; int64_t q0=0;
    for(int64_t y=1;y<m2;y++){ int64_t q=md((((((y+B5)*y+a4)*y+a3)*y+a2)*y+a1)*y,m2); if(q!=q0){ok=0;break;} }
    ok2[a4][a3][a2][a1]=ok; if(ok) any2[a4][a3][a2]=1; }
  int u1=(int)(T1%27), u2=(int)(T2%27);
  for(int a4=0;a4<27;a4++) for(int a3=0;a3<27;a3++) for(int a2=0;a2<27;a2++) for(int a1=0;a1<27;a1++){ int good=0;
    for(int bb=0;bb<27 && !good;bb++){ int ok=1; for(int64_t y=0;y<27;y++){ int64_t q=md((((((y+B5)*y+a4)*y+a3)*y+a2)*y+a1)*y+bb,27); if(q!=u1&&q!=u2){ok=0;break;} } if(ok) good=1; }
    ok3[a4][a3][a2][a1]=good; if(good) any3[a4][a3][a2]=1; }

  LIM=1ULL<<29; bits=malloc(LIM>>4); memset(bits,0xff,LIM>>4); bits[0]&=~2;
  for(uint64_t i=3;i*i<LIM;i+=2) if((bits[i>>4]>>((i>>1)&7))&1) for(uint64_t j=i*i;j<LIM;j+=2*i) bits[j>>4]&=~(1<<((j>>1)&7));
  int pr[32],np=0; for(int p=5;p<=PS;p+=2){ int q=1; for(int d=3;d*d<=p;d+=2) if(p%d==0) q=0; if(q) pr[np++]=p; }
  int nsm=0; while(nsm<np && pr[nsm]<=13) nsm++;
  int64_t M=MOD*5005; int64_t mods[6]={MOD,5,7,11,13}; int64_t crt[6];
  for(int i=0;i<5;i++){ int64_t p=mods[i], Mi=M/p, inv=1; while((Mi%p)*inv%p!=1) inv++; crt[i]=Mi*inv%M; }
  int64_t *cls=malloc(sizeof(int64_t)*M/8), *buf=malloc(sizeof(int64_t)*M/8);
  long long cand=0, tested=0, b1n=0; long long hist[128]={0}; int bestall=0;
  static long long SH[200000][5]; long nsh=0;
  if(getenv("SHAPEFILE")){ FILE *fp=fopen(getenv("SHAPEFILE"),"r"); while(nsh<200000 && fscanf(fp,"%lld %lld %lld %lld %lld",&SH[nsh][0],&SH[nsh][1],&SH[nsh][2],&SH[nsh][3],&SH[nsh][4])==5) nsh++; fclose(fp); NS=nsh; }
  for(long s=0;s<NS;s++){
    /* random root set: nr real roots (2, 4 or 6), the rest complex pairs */
    double complex r[6]; int nr = 2*(1+(int)(urand()*3)); int i=0;
    for(;i<nr;i++) r[i]=(urand()*2-1)*RS;
    for(;i<6;i+=2){ double re=(urand()*2-1)*RS, im=0.5+urand()*6; r[i]=re+im*I; r[i+1]=re-im*I; }
    double complex mean=0; for(i=0;i<6;i++) mean+=r[i]; mean/=6; for(i=0;i<6;i++) r[i]-=mean+B5/6.0;
    double complex c[7]={1,0,0,0,0,0,0};   /* prod (y - r_i), c[k] = coeff of y^k */
    for(i=0;i<6;i++){ for(int k=6;k>=1;k--) c[k]=c[k-1]-r[i]*c[k]; c[0]=-r[i]*c[0]; }
    /* after the loop c[] holds coefficients in reverse: compute properly */
    double complex q[7]={0}; q[0]=1; int deg=0;
    for(i=0;i<6;i++){ double complex nq[7]={0}; for(int k=0;k<=deg;k++){ nq[k+1]+=q[k]; nq[k]-=r[i]*q[k]; } deg++; memcpy(q,nq,sizeof(nq)); }
    int64_t b4=llround(creal(q[4])), b3=llround(creal(q[3])), b2=llround(creal(q[2])), b1c=llround(creal(q[1]));
    if(nsh){ b4=SH[s][0]; b3=SH[s][1]; b2=SH[s][2]; b1c=SH[s][3]; V=SH[s][4]; }
    if(getenv("SHAPE")){ long long x4,x3,x2,x1; sscanf(getenv("SHAPE"),"%lld %lld %lld %lld",&x4,&x3,&x2,&x1); b4=x4;b3=x3;b2=x2;b1c=x1; }
    { int found=0; for(int rad=0;rad<=20 && !found;rad++) for(int d4=-rad;d4<=rad && !found;d4++) for(int d3=-rad;d3<=rad && !found;d3++) for(int d2=-rad;d2<=rad && !found;d2++){
        if(abs(d4)!=rad && abs(d3)!=rad && abs(d2)!=rad) continue;
        int64_t x4=b4+d4,x3=b3+d3,x2=b2+d2; if(any2[md(x4,m2)][md(x3,m2)][md(x2,m2)] && any3[md(x4,27)][md(x3,27)][md(x2,27)]){ b4=x4;b3=x3;b2=x2;found=1; } }
      if(!found) continue; }
    /* h(y) = y^6 + b4 y^4 + b3 y^3 + b2 y^2 on the core window */
    int64_t hy[200]; for(int y=-H;y<=H;y++){ int64_t yy=y; hy[y+H]=(((((yy+B5)*yy+b4)*yy+b3)*yy+b2)*yy)*yy; }
    int64_t R1=getenv("R1")?atoll(getenv("R1")):(K*V)/(H>0?H:1)+1;
    for(int64_t b1=b1c-R1;b1<=b1c+R1;b1++){
      { int sym=1; for(int64_t kk=1;kk<=3;kk++){ int64_t u=kk, w=-1-kk;
          int64_t Nu=(((((u+B5)*u+b4)*u+b3)*u+b2)*u+b1)*u, Nw=(((((w+B5)*w+b4)*w+b3)*w+b2)*w+b1)*w; if(Nu!=Nw){sym=0;break;} }
        if(sym) continue; }
      if(!ok2[md(b4,m2)][md(b3,m2)][md(b2,m2)][md(b1,m2)] || !ok3[md(b4,27)][md(b3,27)][md(b2,27)][md(b1,27)]) continue;
      int64_t gmax=INT64_MIN,gmin=INT64_MAX; for(int y=-H;y<=H;y++){ int64_t g=hy[y+H]+b1*y; if(g>gmax)gmax=g; if(g<gmin)gmin=g; }
      int64_t elo=-K*V-gmin, ehi=K*V-gmax; if(elo>ehi) continue;
      /* mod MOD: values N(y) - b0 mod MOD for y in one period; need a b0 putting all of them in {T1, T2} */
      int64_t v0=-1; int dead=0; int64_t allowedb0[2]; int na=0;
      /* candidate b0 residues: from y = 0, b0 = T1 - q0 or T2 - q0 */
      int64_t qv[432]; for(int64_t y=0;y<MOD;y++){ int64_t yy=y; qv[y]=md((((((yy+B5)*yy+b4)*yy+b3)*yy+b2)*yy+b1)*yy,MOD); }
      for(int t=0;t<2;t++){ int64_t b0r=md((t?T2:T1)-qv[0],MOD); int ok=1;
        for(int64_t y=1;y<MOD;y++){ int64_t w=(qv[y]+b0r)%MOD; if(w!=T1 && w!=T2){ok=0;break;} }
        if(ok) allowedb0[na++]=b0r; }
      (void)v0; (void)dead;
      if(!na) continue; b1n++;
      uint64_t mm[32]; int bad=0;
      for(int k=0;k<np;k++){ int64_t p=pr[k]; uint64_t m=(1ULL<<p)-1;
        for(int64_t y=0;y<p;y++){ int64_t yy=y; int64_t v=md((((((yy+B5)*yy+b4)*yy+b3)*yy+b2)*yy+b1)*yy,p); m&=~(1ULL<<((p-v)%p)); }
        mm[k]=m; if(!m){bad=1;break;} }
      if(bad) continue;
      int nc=0; for(int t=0;t<na;t++) cls[nc++]=crt[0]*allowedb0[t]%M;
      for(int k=0;k<nsm;k++){ int p=pr[k]; int nb=0;
        for(int j=0;j<nc;j++) for(int rr=0;rr<p;rr++) if((mm[k]>>rr)&1) buf[nb++]=(cls[j]+crt[k+1]*rr)%M;
        memcpy(cls,buf,sizeof(int64_t)*nb); nc=nb; }
      for(int j=0;j<nc;j++) for(int64_t b0=cls[j]+((elo-cls[j])/M-1)*M; b0<=ehi; b0+=M){
        if(b0<elo) continue; cand++;
        int pass=1; for(int k=nsm;k<np;k++) if(!((mm[k]>>md(b0,pr[k]))&1)){pass=0;break;}
        if(!pass) continue; tested++;
        int best=0,bs=0,cur=0;
        for(int y=-Y;y<=Y;y++){ int64_t yy=y; int64_t N=(((((((yy+B5)*yy+b4)*yy+b3)*yy+b2)*yy+b1)*yy)+b0); int t=isp(N/K); cur=t?cur+1:0; if(cur>best){best=cur;bs=y-cur+1;} }
        hist[best<127?best:127]++; if(best>bestall) bestall=best;
        if(best>=MINRUN || getenv("DUMP")){ printf("%lld %lld %lld %lld %lld %lld %lld %d %d\n",(long long)K,(long long)B5,(long long)b4,(long long)b3,(long long)b2,(long long)b1,(long long)b0,bs,best); fflush(stdout); }
      }
    }
  }
  fprintf(stderr,"K=%lld shapes %ld b1ok %lld cand %lld tested %lld best %d | runs: ",(long long)K,NS,b1n,cand,tested,bestall);
  for(int L=30;L<128;L++) if(hist[L]) fprintf(stderr,"%d:%lld ",L,hist[L]); fprintf(stderr,"\n"); return 0; }
