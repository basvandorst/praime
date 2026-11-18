#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef unsigned __int128 u128; typedef uint64_t u64; typedef int64_t i64;
#define MAXD 64
static u64 N,NI,ONE,R2;
static inline u64 mont_inv(u64 n){u64 x=n;for(int i=0;i<6;i++)x*=2-n*x;return x;}
static inline u64 redc(u128 t){u64 m=(u64)t*NI;u128 mn=(u128)m*N;u64 r=(u64)((t-mn)>>64);return (t<mn)?r+N:r;}
static inline u64 mm(u64 a,u64 b){return redc((u128)a*b);}
static inline u64 ad(u64 a,u64 b){u64 s=a+b;return (s>=N||s<a)?s-N:s;}
static inline u64 tom(i64 v){u64 r;if(v>=0)r=(u64)v%N;else{r=(u64)(-v)%N;if(r)r=N-r;}return mm(r,R2);}
static void setn(u64 n){N=n;NI=mont_inv(n);ONE=(u64)(((u128)1<<64)%n);R2=(u64)(((u128)ONE*ONE)%n);}
typedef struct{int type;int d;i64 c[MAXD];i64 a;i64 ps[MAXD];int a2;}Mach;
static int test_exp(i64 b){u64 x=tom(b),bm=x;u64 e=N;int top=63-__builtin_clzll(e);for(int i=top-1;i>=0;i--){x=mm(x,x);if((e>>i)&1)x=mm(x,bm);}return x==bm;}
static int test_ring(int r,i64 a){u64 am=tom(a);u64 cur[MAXD],t[MAXD];memset(cur,0,sizeof(u64)*r);
  if(r==1){cur[0]=ad(ONE,am);}else{cur[0]=am;cur[1]=ONE;}
  u64 e=N;int top=63-__builtin_clzll(e);
  for(int i=top-1;i>=0;i--){
    for(int k=0;k<r;k++)t[k]=0;
    for(int x=0;x<r;x++){if(!cur[x])continue;int k=x+x;if(k>=r)k-=r;t[k]=ad(t[k],mm(cur[x],cur[x]));for(int y=x+1;y<r;y++){k=x+y;if(k>=r)k-=r;u64 v=mm(cur[x],cur[y]);t[k]=ad(t[k],ad(v,v));}}
    if((e>>i)&1){for(int k=0;k<r;k++){u64 v=mm(t[k],am);v=ad(v,t[k==0?r-1:k-1]);cur[k]=v;}}
    else for(int k=0;k<r;k++)cur[k]=t[k];}
  u64 want[MAXD];for(int k=0;k<r;k++)want[k]=0;want[0]=am;int j=(int)(N%r);want[j]=ad(want[j],ONE);
  for(int k=0;k<r;k++)if(cur[k]!=want[k])return 0;return 1;}
static int test_trace(Mach*m){int d=m->d;u64 c[MAXD],cur[MAXD],t[2*MAXD];for(int i=0;i<d;i++)c[i]=tom(m->c[i]);
  for(int i=0;i<d;i++)cur[i]=0;if(d==1)cur[0]=c[0];else cur[1]=ONE;
  u64 e=N;int top=63-__builtin_clzll(e);
  for(int i=top-1;i>=0;i--){
    for(int k=0;k<2*d-1;k++)t[k]=0;
    for(int x=0;x<d;x++){if(!cur[x])continue;for(int y=0;y<d;y++)t[x+y]=ad(t[x+y],mm(cur[x],cur[y]));}
    for(int k=2*d-2;k>=d;k--){u64 v=t[k];if(!v)continue;for(int j=0;j<d;j++)t[k-d+j]=ad(t[k-d+j],mm(v,c[j]));}
    if((e>>i)&1){u64 topc=t[d-1];for(int k=d-1;k>=1;k--)t[k]=ad(t[k-1],mm(topc,c[k]));t[0]=mm(topc,c[0]);}
    for(int k=0;k<d;k++)cur[k]=t[k];}
  u64 s=0;for(int k=0;k<d;k++)s=ad(s,mm(cur[k],tom(m->ps[k])));
  return s==tom(m->ps[1]);}
static void polypow(int d,u64*c,u64*base,u64*out){u64 cur[MAXD],t[2*MAXD];for(int i=0;i<d;i++)cur[i]=base[i];u64 e=N;int top=63-__builtin_clzll(e);
  for(int i=top-1;i>=0;i--){for(int k=0;k<2*d-1;k++)t[k]=0;
    for(int x=0;x<d;x++){if(!cur[x])continue;for(int y=0;y<d;y++)t[x+y]=ad(t[x+y],mm(cur[x],cur[y]));}
    for(int k=2*d-2;k>=d;k--){u64 v=t[k];if(!v)continue;for(int j=0;j<d;j++)t[k-d+j]=ad(t[k-d+j],mm(v,c[j]));}
    if((e>>i)&1){for(int k=0;k<2*d-1;k++)cur[k<d?k:0]=cur[k<d?k:0];u64 u[2*MAXD];for(int k=0;k<2*d-1;k++)u[k]=0;
      for(int x=0;x<d;x++){if(!t[x])continue;for(int y=0;y<d;y++){if(!base[y])continue;u[x+y]=ad(u[x+y],mm(t[x],base[y]));}}
      for(int k=2*d-2;k>=d;k--){u64 v=u[k];if(!v)continue;for(int j=0;j<d;j++)u[k-d+j]=ad(u[k-d+j],mm(v,c[j]));}
      for(int k=0;k<d;k++)cur[k]=u[k];}
    else for(int k=0;k<d;k++)cur[k]=t[k];}
  for(int k=0;k<d;k++)out[k]=cur[k];}
static int test_poly(Mach*m){int d=m->d;u64 c[MAXD],b1[MAXD],b2[MAXD],p1[MAXD],p2[MAXD];for(int i=0;i<d;i++){c[i]=tom(m->c[i]);b1[i]=0;b2[i]=0;}u64 am=tom(m->a);
  if(d==1){b2[0]=c[0];b1[0]=ad(c[0],am);}else{b2[1]=ONE;b1[1]=ONE;b1[0]=am;}
  polypow(d,c,b1,p1);polypow(d,c,b2,p2);if(m->type==5){int same=1;for(int k=0;k<d;k++)if(p2[k]!=b2[k])same=0;if(same)return 0;}p2[0]=ad(p2[0],am);for(int k=0;k<d;k++)if(p1[k]!=p2[k])return 0;return 1;}
static int test_gold(int strong){u64 c[2],b[2],p[2];c[0]=ONE;c[1]=ONE;b[0]=0;b[1]=ONE;polypow(2,c,b,p);
  int m5=(int)(N%5);if(m5==2||m5==3){if(strong)return p[0]==ONE&&p[1]==N-ONE;return ad(p[0],p[1])==0;}
  if(m5==1||m5==4){if(strong)return p[0]==0&&p[1]==ONE;return p[1]==ONE||1?(p[1]==ONE&&p[0]==0)||0:0;}return 0;}
static u64 RH[MAXD+1];
static u64 gcd64(u64 a,u64 b){while(b){u64 t=a%b;a=b;b=t;}return a;}
static int misfit(void){for(int r=2;r<MAXD;r++){if(gcd64(N,r)!=1)continue;u64 m=N%r;if((m*m)%r!=1)return r;}return 0;}
static int test_aring(i64 a){int r=misfit();if(!r)return 1;RH[r]++;return test_ring(r,a);}
static u64 CT[MAXD+1],CL[MAXD+1];
static int test_fib(void){u64 c[2],b[2],p[2];c[0]=ONE;c[1]=ONE;b[0]=0;b[1]=ONE;polypow(2,c,b,p);int m5=(int)(N%5);if(m5==0)return 0;if(m5==1||m5==4)return p[0]==0;return ad(p[0],p[1])==0;}
static inline u64 sb(u64 a,u64 b){return a>=b?a-b:a+(N-b);}
static int test_nega(int h){u64 cur[MAXD],t[MAXD];memset(cur,0,sizeof(u64)*h);cur[0]=ONE;cur[1]=ONE;u64 e=N;int top=63-__builtin_clzll(e);
  for(int i=top-1;i>=0;i--){for(int k=0;k<h;k++)t[k]=0;
    for(int x=0;x<h;x++){if(!cur[x])continue;int k=x+x;u64 v=mm(cur[x],cur[x]);if(k>=h)t[k-h]=sb(t[k-h],v);else t[k]=ad(t[k],v);for(int y=x+1;y<h;y++){k=x+y;v=mm(cur[x],cur[y]);v=ad(v,v);if(k>=h)t[k-h]=sb(t[k-h],v);else t[k]=ad(t[k],v);}}
    if((e>>i)&1){cur[0]=sb(t[0],t[h-1]);for(int k=1;k<h;k++)cur[k]=ad(t[k],t[k-1]);}else for(int k=0;k<h;k++)cur[k]=t[k];}
  u64 want[MAXD];for(int k=0;k<h;k++)want[k]=0;want[0]=ONE;int j=(int)(N%(2*h));if(j<h)want[j]=ad(want[j],ONE);else want[j-h]=sb(want[j-h],ONE);
  for(int k=0;k<h;k++)if(cur[k]!=want[k])return 0;return 1;}
static u64 NP[MAXD+1];
static const int SP[]={3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,0};
static u64 FC[256],FNO;
static int jac_small(int c){int m=(int)(N%c);if(m==0)return 0;int e=(c-1)/2,r=1,b=m;while(e){if(e&1)r=r*b%c;b=b*b%c;e>>=1;}int leg=(r==1)?1:-1;if((c&3)==3&&(N&3)==3)leg=-leg;return leg;}
static int test_frob(void){int c=0;for(int i=0;SP[i];i++){int j=jac_small(SP[i]);if(j==0){if(N==(u64)SP[i])continue;return 0;}if(j==-1){c=SP[i];break;}}
  if(!c){FNO++;return 0;}FC[c]++;u64 cm=tom(c),a=ONE,b=ONE;u64 e=N;int top=63-__builtin_clzll(e);
  for(int i=top-1;i>=0;i--){u64 a2=mm(a,a),b2=mm(b,b),ab=mm(a,b);a=ad(a2,mm(cm,b2));b=ad(ab,ab);if((e>>i)&1){u64 na=ad(a,mm(cm,b));b=ad(a,b);a=na;}}
  return a==ONE&&b==N-ONE;}
static u64 FT,FL;
static u64 PT,PL,PS;
static int test_census(Mach*m){if(m->a2){FT++;if(test_frob()){FL++;printf("two-cell liar %llu\n",(unsigned long long)N);fflush(stdout);}}{int m5=(int)(N%5);if(m5==2||m5==3){PT++;if(test_gold(0)){PL++;printf("pair liar %llu\n",(unsigned long long)N);fflush(stdout);}if(test_gold(1))PS++;}}for(int i=0;i<m->d;i++){int r=(int)m->c[i];if(gcd64(N,r)!=1)continue;u64 q=N%r;if((q*q)%r==1)continue;CT[r]++;if(r>=16&&(r&(r-1))==0){if(!test_nega(r/2))continue;NP[r]++;}if(test_ring(r,m->a)){CL[r]++;printf("ring %d liar %llu\n",r,(unsigned long long)N);fflush(stdout);}}return 0;}
static int run(Mach*ms,int nm){for(int i=0;i<nm;i++){Mach*m=&ms[i];int ok=m->type==0?test_exp(m->a):m->type==1?test_ring(m->d,m->a):m->type==2?test_trace(m):m->type==9?test_frob():m->type==8?test_fib():m->type==7?test_census(m):m->type==4?test_aring(m->a):m->type==6?test_gold((int)m->a):test_poly(m);if(!ok)return 0;}return 1;}
int main(int argc,char**argv){
  int ai=1;u64 K=0;int quiet=0;if(!strcmp(argv[ai],"-k")){K=strtoull(argv[ai+1],0,10);ai+=2;}if(!strcmp(argv[ai],"-q")){quiet=1;ai++;}
  int mode=!strcmp(argv[ai],"list")?0:!strcmp(argv[ai],"stdin")?2:1;ai++;const char*fn=0;u64 lo=0,hi=0;if(mode==0)fn=argv[ai++];else if(mode==1){lo=strtoull(argv[ai++],0,10);hi=strtoull(argv[ai++],0,10);}
  Mach ms[8];int nm=0;
  while(ai<argc){Mach*m=&ms[nm++];memset(m,0,sizeof(Mach));
    if(!strcmp(argv[ai],"exp")){m->type=0;m->a=atoll(argv[ai+1]);ai+=2;}
    else if(!strcmp(argv[ai],"ring")){m->type=1;m->d=atoi(argv[ai+1]);m->a=atoll(argv[ai+2]);ai+=3;}
    else if(!strcmp(argv[ai],"trace")){m->type=2;m->d=atoi(argv[ai+1]);ai+=2;for(int i=m->d-1;i>=0;i--)m->c[i]=atoll(argv[ai++]);
      int d=m->d;m->ps[0]=d;for(int k=1;k<d;k++){i64 v=0;for(int j=1;j<k;j++)v+=m->c[d-j]*m->ps[k-j];v+=(i64)k*m->c[d-k];m->ps[k]=v;}
      if(d==1)m->ps[1]=m->c[0];}
    else if(!strcmp(argv[ai],"gold")){m->type=6;m->a=atoll(argv[ai+1]);ai+=2;}
    else if(!strcmp(argv[ai],"frob")){m->type=9;ai+=1;}
    else if(!strcmp(argv[ai],"census2")){m->type=7;m->a=1;m->a2=1;m->d=atoi(argv[ai+1]);for(int i=0;i<m->d;i++)m->c[i]=atoll(argv[ai+2+i]);ai+=2+m->d;}
    else if(!strcmp(argv[ai],"census")){m->type=7;m->a=1;m->d=atoi(argv[ai+1]);for(int i=0;i<m->d;i++)m->c[i]=atoll(argv[ai+2+i]);ai+=2+m->d;}
    else if(!strcmp(argv[ai],"fib")){m->type=8;ai+=1;}
    else if(!strcmp(argv[ai],"aring")){m->type=4;m->a=atoll(argv[ai+1]);ai+=2;}
    else if(!strcmp(argv[ai],"poly")||!strcmp(argv[ai],"polym")){m->type=!strcmp(argv[ai],"polym")?5:3;m->d=atoi(argv[ai+1]);m->a=atoll(argv[ai+2]);ai+=3;for(int i=m->d-1;i>=0;i--)m->c[i]=atoll(argv[ai++]);}
    else{fprintf(stderr,"bad machine %s\n",argv[ai]);return 1;}}
  u64 cnt=0,tested=0;
  if(mode==0){FILE*f=fopen(fn,"rb");u64 n;while(fread(&n,8,1,f)==1){if(n%2==0||n<3)continue;setn(n);tested++;if(run(ms,nm)){cnt++;if(!quiet)printf("%llu\n",(unsigned long long)n);if(K&&cnt>=K)break;}}fclose(f);}
  else if(mode==2){char buf[64];u64 last=0;while(fgets(buf,64,stdin)){u64 n=strtoull(buf,0,10);if(n%2==0||n<3)continue;setn(n);tested++;last=n;if(run(ms,nm)){cnt++;if(!quiet){printf("%llu\n",(unsigned long long)n);fflush(stdout);}if(K&&cnt>=K)break;}if((tested&0xFFFFFF)==0)fprintf(stderr,"  .. %llu tested, at %llu, liars %llu\n",(unsigned long long)tested,(unsigned long long)n,(unsigned long long)cnt);}fprintf(stderr,"last number read %llu\n",(unsigned long long)last);}
  else{u64 sq=(u64)sqrt((double)hi)+2;char*small=calloc(sq+1,1);u64 np=0;u64*pr=malloc(sizeof(u64)*(sq/2+10));
    for(u64 i=2;i<=sq;i++){if(!small[i]){pr[np++]=i;for(u64 j=i*i;j<=sq;j+=i)small[j]=1;}}
    const u64 BL=1<<20;char*seg=malloc(BL);if(lo%2==0)lo++;if(lo<3)lo=3;int stop=0;
    for(u64 base=lo;base<hi&&!stop;base+=2*BL){memset(seg,0,BL);u64 end=base+2*BL;if(end>hi)end=hi;
      for(u64 k=1;k<np;k++){u64 p=pr[k];if(p*p>=end)break;u64 st=p*p;if(st<base)st=((base+p-1)/p)*p;if(st%2==0)st+=p;for(u64 j=st;j<end;j+=2*p)seg[(j-base)>>1]=1;}
      for(u64 n=base;n<end;n+=2){if(!seg[(n-base)>>1])continue;setn(n);tested++;if(run(ms,nm)){cnt++;if(!quiet)printf("%llu\n",(unsigned long long)n);if(K&&cnt>=K){stop=1;break;}}}}}
  fprintf(stderr,"tested %llu liars %llu\n",(unsigned long long)tested,(unsigned long long)cnt);for(int r=2;r<MAXD;r++)if(RH[r])fprintf(stderr,"  ring %d used %llu times\n",r,(unsigned long long)RH[r]);if(FT)fprintf(stderr,"  two-cell wheel with the doubling wheel: tested %llu liars %llu no wheel found %llu\n",(unsigned long long)FT,(unsigned long long)FL,(unsigned long long)FNO);for(int c=0;c<256;c++)if(FC[c])fprintf(stderr,"  two-cell wheel c = %d used %llu times\n",c,(unsigned long long)FC[c]);if(PT)fprintf(stderr,"  pair (doubling and Fibonacci wheel, n = 2 or 3 mod 5): tested %llu liars %llu strong liars %llu\n",(unsigned long long)PT,(unsigned long long)PL,(unsigned long long)PS);for(int r=2;r<MAXD;r++)if(NP[r])fprintf(stderr,"  ring %d: %llu passed the half ring with a sign first\n",r,(unsigned long long)NP[r]);for(int r=2;r<MAXD;r++)if(CT[r])fprintf(stderr,"  census ring %d: tested %llu liars %llu\n",r,(unsigned long long)CT[r],(unsigned long long)CL[r]);return 0;}
