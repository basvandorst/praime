#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
int n; int adj[256][256], deg[256]; int vis[256], path[256], avail[256];
long long nodes, lim;
int isp(long x){if(x<2)return 0;for(long d=2;d*d<=x;d++)if(x%d==0)return 0;return 1;}
int dfs(int len){
 if(len==n)return 1;
 if(++nodes>lim)return -1;
 int cur=path[len-1];
 int forced=-1, ones=0;
 for(int v=1;v<=n;v++){if(vis[v])continue; int a=avail[v]; int adjcur=0;
  for(int k=0;k<deg[v];k++)if(adj[v][k]==cur){adjcur=1;break;}
  int tot=a+adjcur;
  if(tot==0)return 0;
  if(a==0){ if(forced!=-1&&forced!=v)return 0; forced=v; }
  else if(a==1&&!adjcur){ ones++; if(ones>1)return 0; }
 }
 int cand[256],nc=0;
 for(int k=0;k<deg[cur];k++){int v=adj[cur][k];if(!vis[v])cand[nc++]=v;}
 if(forced!=-1){ int ok=0; for(int i=0;i<nc;i++)if(cand[i]==forced)ok=1; if(!ok)return 0; if(len+1<n){ cand[0]=forced; nc=1; } }
 for(int i=1;i<nc;i++){int x=cand[i],j=i;while(j>0&&avail[cand[j-1]]>avail[x]){cand[j]=cand[j-1];j--;}cand[j]=x;}
 for(int i=0;i<nc;i++){int v=cand[i];
  vis[v]=1;path[len]=v;for(int k=0;k<deg[v];k++)avail[adj[v][k]]--;
  int r=dfs(len+1);
  for(int k=0;k<deg[v];k++)avail[adj[v][k]]++;vis[v]=0;
  if(r)return r;}
 return 0;}
int main(int argc,char**argv){n=atoi(argv[1]);lim=atoll(argv[2]);
 for(int i=1;i<=n;i++){deg[i]=0;for(int j=1;j<=n;j++)if(i!=j&&isp(i+j)&&isp((long)i*i+(long)j*j))adj[i][deg[i]++]=j;}
 for(int i=1;i<=n;i++)avail[i]=deg[i];
 int res=0;
 for(int s=1;s<=n&&res!=1;s++){ if(n%2==1&&s%2==0)continue;
  memset(vis,0,sizeof vis);for(int i=1;i<=n;i++)avail[i]=deg[i];
  vis[s]=1;path[0]=s;for(int k=0;k<deg[s];k++)avail[adj[s][k]]--;
  nodes=0;int r=dfs(1); if(r==1)res=1; else if(r==-1)res=-1>res?res:-1;
  if(r==-1)res=-1;
  if(r==1)break;}
 if(res==1){printf("n=%d FOUND:",n);for(int i=0;i<n;i++)printf(" %d",path[i]);printf("\n");}
 else printf("n=%d %s\n",n,res==-1?"TIMEOUT(some start)":"NONE (exhaustive)");
 return 0;}
