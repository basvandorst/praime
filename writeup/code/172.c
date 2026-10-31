#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int K,w,m,L; int pr[]={2,3,5,7,11,13,17,19,23,29,31,37,41};
int left[64],right[64]; int nl,nr; int PL[64],PR[64]; long long found=0,maxf; long long nodes=0;
int cnt[4096]; char cov[2048]; int sh1=-2,sh2=-2;
static int check_full(){int seq[64],k=0;for(int i=0;i<nl;i++)seq[k++]=left[i];for(int i=nr-1;i>=0;i--)seq[k++]=right[i];
 memset(cov,0,L+1);for(int a=0;a<m;a++){int s=0;for(int b=a;b<m;b++){s+=seq[b];cov[s]=1;}}
 for(int s=1;s<=L;s++)if(!cov[s])return 0;
 if(found<maxf){for(int i=0;i<m;i++)printf("%d%c",seq[i],i<m-1?'.':'\n');fflush(stdout);}
 return 1;}
static int prune(){int lim=PL[nl]<PR[nr]?PL[nl]:PR[nr];for(int s=0;s<=lim;s++)if(!cnt[s])return 1;return 0;}
void rec(int ones,int mask,int depth){
 nodes++;
 if(nl+nr==m){ if(check_full())found++; return;}
 if(prune())return;
 int toleft = PL[nl]<=PR[nr];
 for(int e=-1;e<K;e++){int v;
  if(depth==0&&sh1>-2&&e!=sh1)continue;
  if(depth==1&&sh2>-2&&e!=sh2)continue;
  if(e<0){if(!ones)continue;v=1;}else{if(!(mask>>e&1))continue;v=pr[e];}
  if(toleft){left[nl]=v;PL[nl+1]=PL[nl]+v;nl++;for(int j=0;j<=nr;j++)cnt[PL[nl]+PR[j]]++;}
  else{right[nr]=v;PR[nr+1]=PR[nr]+v;nr++;for(int i=0;i<=nl;i++)cnt[PL[i]+PR[nr]]++;}
  rec(ones-(e<0),e<0?mask:mask&~(1<<e),depth+1);
  if(toleft){for(int j=0;j<=nr;j++)cnt[PL[nl]+PR[j]]--;nl--;}else{for(int i=0;i<=nl;i++)cnt[PL[i]+PR[nr]]--;nr--;}
 }}
int main(int c,char**v){K=atoi(v[1]);w=atoi(v[2]);maxf=c>3?atoll(v[3]):5;if(c>5){sh1=atoi(v[4]);sh2=atoi(v[5]);}
 m=K+w;L=w;for(int i=0;i<K;i++)L+=pr[i];
 PL[0]=PR[0]=0;nl=nr=0;cnt[0]=1;rec(w,(1<<K)-1,0);
 printf("K=%d w=%d L=%d shard=%d,%d sequences=%lld nodes=%lld\n",K,w,L,sh1,sh2,found,nodes);return 0;}
