#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned __int128 u128; typedef uint64_t u64; typedef int64_t i64;
#define MR 256
#define ME 24
static u64 N,NI,ONE,R2;
static inline u64 mont_inv(u64 n){u64 x=n;for(int i=0;i<6;i++)x*=2-n*x;return x;}
static inline u64 redc(u128 t){u64 m=(u64)t*NI;u128 mn=(u128)m*N;u64 r=(u64)((t-mn)>>64);return (t<mn)?r+N:r;}
static inline u64 mm(u64 a,u64 b){return redc((u128)a*b);}
static inline u64 sb(u64 a,u64 b){return a>=b?a-b:a+(N-b);}
static inline u64 ad(u64 a,u64 b){u64 s=a+b;return (s>=N||s<a)?s-N:s;}
static inline u64 tom(i64 v){u64 r;if(v>=0)r=(u64)v%N;else{r=(u64)(-v)%N;if(r)r=N-r;}return mm(r,R2);}
static void setn(u64 n){N=n;NI=mont_inv(n);ONE=(u64)(((u128)1<<64)%n);R2=(u64)(((u128)ONE*ONE)%n);}
static inline u64 fin(u128 lo,u64 hi){u128 s=lo+(u128)hi*R2;if(s<lo)s+=R2;u64 sh=(u64)(s>>64);if(sh>=N)sh%=N;return redc(((u128)sh<<64)|(u64)s);}
typedef struct{int np[2],ns[2];unsigned char p[2][MR][2],s[2][MR];}Out;
typedef struct{int len,nout,kind,r;Out o[MR];}Ring;
static Ring*mk(int kind,int r){Ring*R=calloc(1,sizeof(Ring));R->kind=kind;R->r=r;int len=r;R->len=len;R->nout=kind==1?(r+1)/2:len;
  for(int k=0;k<R->nout;k++){Out*o=&R->o[k];for(int i=0;i<len;i++){int j,sg=0;
      if(kind==2){j=k-i;if(j<0){j+=len;sg=1;}}else{j=((k-i)%len+len)%len;}
      if(i<j){o->p[sg][o->np[sg]][0]=i;o->p[sg][o->np[sg]][1]=j;o->np[sg]++;}else if(i==j){o->s[sg][o->ns[sg]++]=i;}}}
  return R;}
static inline u64 acc(Out*o,int sg,const u64*a){u128 lo=0;u64 hi=0;for(int t=0;t<o->np[sg];t++){u128 p=(u128)a[o->p[sg][t][0]]*a[o->p[sg][t][1]];lo+=p;hi+=lo<p;}
  hi=(hi<<1)|(u64)(lo>>127);lo<<=1;for(int t=0;t<o->ns[sg];t++){u64 v=a[o->s[sg][t]];u128 p=(u128)v*v;lo+=p;hi+=lo<p;}return fin(lo,hi);}
static void sq(Ring*R,const u64*a,u64*c){int len=R->len;for(int k=0;k<R->nout;k++){Out*o=&R->o[k];u64 v=acc(o,0,a);if(R->kind==2)v=sb(v,acc(o,1,a));c[k]=v;}
  if(R->kind==1)for(int k=R->nout;k<len;k++)c[k]=c[len-k];}
static int ringtest(Ring*R){int len=R->len,kind=R->kind;u64 cur[MR],t[MR];memset(cur,0,sizeof(u64)*len);
  if(kind==1){cur[0]=ad(ONE,ONE);cur[1]=ONE;cur[len-1]=ONE;}else{cur[0]=ONE;cur[1]=ONE;}
  u64 e=N;int top=63-__builtin_clzll(e);
  for(int i=top-1;i>=0;i--){sq(R,cur,t);
    if((e>>i)&1){
      if(kind==0)for(int k=0;k<len;k++)cur[k]=ad(t[k],t[k==0?len-1:k-1]);
      else if(kind==1)for(int k=0;k<len;k++)cur[k]=ad(ad(t[k],t[k]),ad(t[k==0?len-1:k-1],t[k==len-1?0:k+1]));
      else{cur[0]=sb(t[0],t[len-1]);for(int k=1;k<len;k++)cur[k]=ad(t[k],t[k-1]);}}
    else memcpy(cur,t,sizeof(u64)*len);}
  if(0){printf("%llu %d %d",(unsigned long long)N,R->r,kind);for(int k=0;k<len;k++)printf(" %llu",(unsigned long long)redc(cur[k]));printf("\n");}
  u64 want[MR];memset(want,0,sizeof(u64)*len);
  if(kind==0){want[0]=ONE;int j=(int)(N%len);want[j]=ad(want[j],ONE);}
  else if(kind==1){want[0]=ad(ONE,ONE);int j=(int)(N%len);want[j]=ad(want[j],ONE);want[(len-j)%len]=ad(want[(len-j)%len],ONE);}
  else{want[0]=ONE;int j=(int)(N%(2*len));if(j<len)want[j]=ad(want[j],ONE);else want[j-len]=sb(want[j-len],ONE);}
  for(int k=0;k<len;k++)if(cur[k]!=want[k])return 0;return 1;}
typedef struct{int r,e,f,trivial;int cos[MR];int M[ME][ME][ME];i64 u[ME];int no[ME],nd[ME],ng[ME];unsigned short lo[ME][MR],ld[ME][MR],lg[ME][2*MR];}Fold;
typedef struct{int r,ne;Fold*F[8];int pick[MR],first;Ring*full,*half;u64 T,S[8],U[8],small,LO,LS,HOME;}RingF;
static int isprime(int r){if(r<2)return 0;for(int d=2;d*d<=r;d++)if(r%d==0)return 0;return 1;}
static Fold*mkfold(int r,int e,const int*ind){Fold*F=calloc(1,sizeof(Fold));F->r=r;F->e=e;F->f=(r-1)/e;
  for(int x=1;x<r;x++)F->cos[x]=ind[x]%e;int rep[ME];for(int x=r-1;x>=1;x--)rep[F->cos[x]]=x;
  for(int i=0;i<e;i++)for(int j=0;j<e;j++){int v[MR];memset(v,0,sizeof(v));for(int x=1;x<r;x++)if(F->cos[x]==i)for(int y=1;y<r;y++)if(F->cos[y]==j)v[(x+y)%r]++;
    for(int z=1;z<r;z++)if(v[z]!=v[rep[F->cos[z]]]){fprintf(stderr,"table not constant on a class\n");exit(1);}
    if(v[0]!=(i==j?F->f:0)){fprintf(stderr,"unexpected constant term\n");exit(1);}
    for(int k=0;k<e;k++){F->M[i][j][k]=v[rep[k]]-v[0];for(int c=0;c<v[rep[k]];c++){F->lg[k][F->ng[k]++]=i*ME+j;if(i<j)F->lo[k][F->no[k]++]=i*ME+j;else if(i==j)F->ld[k][F->nd[k]++]=i*ME+j;}}}
  i64 P[MR],Q[MR];memset(P,0,sizeof(P));P[0]=1;
  for(int h=1;h<r;h++)if(F->cos[h]==0){for(int z=0;z<r;z++)Q[z]=P[z]+P[((z-h)%r+r)%r];memcpy(P,Q,sizeof(P));}
  for(int z=1;z<r;z++)if(P[z]!=P[rep[F->cos[z]]]){fprintf(stderr,"base not constant on a class\n");exit(1);}
  for(int k=0;k<e;k++)F->u[k]=P[rep[k]]-P[0];
  F->trivial=1;for(int k=1;k<e;k++)if(F->u[k]!=F->u[0])F->trivial=0;
  return F;}
static void mul(Fold*F,const u64*a,const u64*b,u64*out,int square){int e=F->e;u64 bR[ME],T[ME*ME];u128 dg=0;
  for(int j=0;j<e;j++)bR[j]=mm(b[j],R2);
  if(square){for(int i=0;i<e;i++)for(int j=i;j<e;j++)T[i*ME+j]=mm(a[i],bR[j]);}else{for(int i=0;i<e;i++)for(int j=0;j<e;j++)T[i*ME+j]=mm(a[i],bR[j]);}
  for(int i=0;i<e;i++)dg+=T[i*ME+i];
  u64 neg=redc(dg*(u64)F->f);
  for(int k=0;k<e;k++){u128 pos=0;
    if(square){const unsigned short*l=F->lo[k];int c=F->no[k];for(int t=0;t<c;t++)pos+=T[l[t]];pos<<=1;l=F->ld[k];c=F->nd[k];for(int t=0;t<c;t++)pos+=T[l[t]];}
    else{const unsigned short*l=F->lg[k];int c=F->ng[k];for(int t=0;t<c;t++)pos+=T[l[t]];}
    out[k]=sb(redc(pos),neg);}}
static int DUMP=0;
static int foldtest(Fold*F){int e=F->e;u64 B[ME],a[ME],t[ME];for(int k=0;k<e;k++){B[k]=tom(F->u[k]);a[k]=B[k];}
  u64 ex=N;int top=63-__builtin_clzll(ex);
  for(int i=top-1;i>=0;i--){mul(F,a,a,t,1);if((ex>>i)&1)mul(F,t,B,a,0);else memcpy(a,t,sizeof(u64)*e);}
  if(DUMP){printf("D %llu %d %d",(unsigned long long)N,F->r,e);for(int k=0;k<e;k++)printf(" %llu",(unsigned long long)redc(a[k]));printf("\n");}
  int c=F->cos[N%F->r];for(int k=0;k<e;k++)if(a[(k+c)%e]!=B[k])return 0;return 1;}
static RingF*mkring(int r){RingF*R=calloc(1,sizeof(RingF));R->r=r;int g=2;for(;;g++){int x=1,ok=1;for(int k=1;k<r-1;k++){x=x*g%r;if(x==1){ok=0;break;}}if(ok)break;}
  int ind[MR];{int x=1;for(int k=0;k<r-1;k++){ind[x]=k;x=x*g%r;}}
  int m=(r-1)/2,es[16],ne=0;for(int q=2;q<=m;q++)if(isprime(q)&&m%q==0){int e=q;while(m%e==0){es[ne++]=e;e*=q;}}
  for(int i=0;i<ne;i++)for(int j=i+1;j<ne;j++)if(es[j]<es[i]){int t=es[i];es[i]=es[j];es[j]=t;}
  R->ne=ne;for(int i=0;i<ne;i++){if(es[i]>ME){fprintf(stderr,"ring %d needs a fold of %d cells, too many\n",r,es[i]);exit(1);}R->F[i]=mkfold(r,es[i],ind);}
  for(int x=1;x<r;x++){R->pick[x]=-1;for(int i=0;i<ne;i++)if(!R->F[i]->trivial&&ind[x]%es[i]!=0){R->pick[x]=i;break;}}
  R->first=-1;for(int i=0;i<ne;i++)if(!R->F[i]->trivial){R->first=i;break;}
  R->full=mk(0,r);R->half=mk(1,r);
  return R;}
static u64 gcd64(u64 a,u64 b){while(b){u64 t=a%b;a=b;b=t;}return a;}
int main(int argc,char**argv){int ai=1,K=1,I=0,all=0,fit=0,show=0,three=0;
  while(ai<argc&&argv[ai][0]=='-'){if(!strcmp(argv[ai],"-s")){K=atoi(argv[ai+1]);I=atoi(argv[ai+2]);ai+=3;}else if(!strcmp(argv[ai],"-all")){all=1;ai++;}else if(!strcmp(argv[ai],"-fit")){fit=1;ai++;}else if(!strcmp(argv[ai],"-show")){show=1;ai++;}else if(!strcmp(argv[ai],"-three")){three=1;ai++;}else if(!strcmp(argv[ai],"-dump")){DUMP=1;ai++;}else return 1;}
  if(three){int nw=0,wr[64];Fold*W[64];u64 use[64]={0},tested=0,away=0,nowheel=0,home=0,last=0,line=0;
    for(int r=7;r<MR&&nw<64;r++)if(r%3==1&&isprime(r)){int g=2;for(;;g++){int x=1,ok=1;for(int k=1;k<r-1;k++){x=x*g%r;if(x==1){ok=0;break;}}if(ok)break;}
        int ind[MR];{int x=1;for(int k=0;k<r-1;k++){ind[x]=k;x=x*g%r;}}Fold*F=mkfold(r,3,ind);if(!F->trivial){wr[nw]=r;W[nw++]=F;}}
    char buf[64];
    while(fgets(buf,64,stdin)){u64 n=strtoull(buf,0,10);line++;if((int)(line%K)!=I%K)continue;if(n%2==0||n<3)continue;last=n;tested++;int res=-1;
      for(int i=0;i<nw;i++){int r=wr[i];int x=(int)(n%r);if(x==0){if(n==(u64)r)continue;res=0;away++;break;}if(W[i]->cos[x]==0)continue;setn(n);use[i]++;res=foldtest(W[i]);if(res){home++;printf("three-cell wheel %d home %llu\n",r,(unsigned long long)n);fflush(stdout);}break;}
      if(res<0)nowheel++;
      if((tested&0x3FFFFF)==0)fprintf(stderr,"  .. %llu numbers, at %llu\n",(unsigned long long)tested,(unsigned long long)n);}
    fprintf(stderr,"three-cell wheel: numbers read %llu, last %llu, came home %llu, turned away (divisible by a wheel size) %llu, no wheel %llu\n",(unsigned long long)tested,(unsigned long long)last,(unsigned long long)home,(unsigned long long)away,(unsigned long long)nowheel);
    fprintf(stderr,"  wheel sizes:");for(int i=0;i<nw;i++)fprintf(stderr," %d",wr[i]);fprintf(stderr,"\n");
    for(int i=0;i<nw;i++)if(use[i])fprintf(stderr,"  wheel from the ring of %d cells used %llu times\n",wr[i],(unsigned long long)use[i]);
    return 0;}
  int nr=argc-ai;RingF*Rs[64];for(int i=0;i<nr;i++){int r=atoi(argv[ai+i]);if(r<5||r>=MR||!isprime(r))return 1;Rs[i]=mkring(r);}
  char buf[64];u64 line=0,tested=0,last=0;
  while(fgets(buf,64,stdin)){u64 n=strtoull(buf,0,10);line++;if((int)(line%K)!=I%K)continue;if(n%2==0||n<3)continue;last=n;tested++;int set=0;
    for(int i=0;i<nr;i++){RingF*R=Rs[i];int r=R->r;if(gcd64(n,r)!=1)continue;int x=(int)(n%r);
      if(all){if(!set){setn(n);set=1;}for(int j=0;j<R->ne;j++)printf("%c %d %d %llu\n",foldtest(R->F[j])?'P':'F',r,R->F[j]->e,(unsigned long long)n);continue;}
      int fits=(x==1||x==r-1);if(fits&&!fit)continue;R->T++;
      if(!set){setn(n);set=1;}int j=fits?R->first:R->pick[x],go;
      if(j>=0){R->U[j]++;go=foldtest(R->F[j]);if(go){R->S[j]++;if(show)printf("fold %d %d passed by %llu\n",r,R->F[j]->e,(unsigned long long)n);}}
      else{R->LO++;go=ringtest(R->half);if(go)R->LS++;}
      if(go&&ringtest(R->full)){R->HOME++;printf("ring %d liar %llu\n",r,(unsigned long long)n);fflush(stdout);}}
    if((tested&0x3FFFFF)==0)fprintf(stderr,"  .. %llu numbers, at %llu\n",(unsigned long long)tested,(unsigned long long)n);}
  fprintf(stderr,"numbers read %llu, last %llu\n",(unsigned long long)tested,(unsigned long long)last);
  for(int i=0;i<nr;i++){RingF*R=Rs[i];fprintf(stderr,"  ring %d: tested %llu;",R->r,(unsigned long long)R->T);
    for(int j=0;j<R->ne;j++){if(R->F[j]->trivial)fprintf(stderr," fold of %d cells: says nothing;",R->F[j]->e);else fprintf(stderr," fold of %d cells: %llu tested, %llu passed;",R->F[j]->e,(unsigned long long)R->U[j],(unsigned long long)R->S[j]);}
    fprintf(stderr," no fold: %llu tested, %llu passed the half ring; came home on the full ring: %llu\n",(unsigned long long)R->LO,(unsigned long long)R->LS,(unsigned long long)R->HOME);}
  return 0;}
