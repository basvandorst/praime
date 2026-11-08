// Puzzle 1132: longest chain of k-digit primes, neighbours differing in one digit.
// Build: gcc -O2 -o p1132_path p1132_path.c ; graph file from p1132_graph.py K (writes gK.txt)
// Run: ./p1132_path g6.txt SEED MAXROT -1 -1 out.txt   (seed 1 gives 68873 for k = 6, 8360 for k = 5)
// Posa rotation-extension from a degree-1 prime, all other degree-1 primes but one removed, plus
// insertion of unvisited primes between consecutive path primes (direct, or with one 2-opt reversal).
// Posa rotation-extension + insertion of unvisited vertices between consecutive adjacent path vertices.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int n, *deg, **adj; static long long *lab;
static unsigned long long rs; static unsigned rnd(){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return (unsigned)(rs>>11); }
static int *ok,*vis,*pos,*p,len,target,ea,eb;
static int isadj(int a,int b){ for(int j=0;j<deg[a];j++) if(adj[a][j]==b) return 1; return 0; }
static void insert_at(int i,int v){ // put v between p[i] and p[i+1]
  memmove(p+i+2,p+i+1,sizeof(int)*(len-i-1)); p[i+1]=v; len++; vis[v]=1; for(int t=i+1;t<len;t++) pos[p[t]]=t; }
static void rev(int a,int b){ for(;a<b;a++,b--){ int t=p[a]; p[a]=p[b]; p[b]=t; pos[p[a]]=a; pos[p[b]]=b; } }
static int insertion_pass(){
  int ins=0;
  for(int v=0;v<n;v++){ if(!ok[v]||vis[v]||v==eb) continue;
    int done=0;
    for(int j=0;j<deg[v]&&!done;j++){ int u=adj[v][j]; if(!vis[u]) continue; int i=pos[u];
      if(i+1<len && isadj(v,p[i+1])){ insert_at(i,v); ins++; done=1; break; }
      if(i>0 && isadj(v,p[i-1])){ insert_at(i-1,v); ins++; done=1; break; } }
    // 2-opt insertion: a=p[i], b=p[j] neighbours of v, i<j
    for(int x=0;x<deg[v]&&!done;x++) for(int y=0;y<deg[v]&&!done;y++){ if(x==y) continue;
      int a=adj[v][x], b=adj[v][y]; if(!vis[a]||!vis[b]) continue; int i=pos[a], j=pos[b]; if(i>=j-1) continue;
      // (i) ..a, v, b, p[j-1] .. p[i+1], p[j+1].. needs p[i+1]~p[j+1] or j==len-1
      if(j==len-1 || isadj(p[i+1],p[j+1])){ rev(i+1,j); insert_at(i,v); ins++; done=1; break; }
      // (ii) ..p[i-1], p[j-1] .. p[i], v, b.. needs p[i-1]~p[j-1] (i>0)
      if(i>0 && isadj(p[i-1],p[j-1])){ rev(i,j-1); insert_at(j-1,v); ins++; done=1; break; } }
  }
  return ins; }
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); rs=strtoull(argv[2],0,10)*2654435761ULL+1; long long maxrot=atoll(argv[3]);
  ea = argc>4? atoi(argv[4]) : -1; eb = argc>5? atoi(argv[5]) : -1;
  if(fscanf(f,"%d",&n)!=1) return 1; deg=calloc(n,4); adj=malloc(sizeof(int*)*n); lab=malloc(8*n);
  char line[4096]; if(!fgets(line,sizeof line,f)) return 1;
  for(int i=0;i<n;i++){ if(!fgets(line,sizeof line,f)) return 1; char*s=line; lab[i]=strtoll(s,&s,10); int t[64],c=0; while(1){ while(*s==' ')s++; if(*s<'0'||*s>'9')break; t[c++]=strtol(s,&s,10);} deg[i]=c; adj[i]=malloc(4*c+4); memcpy(adj[i],t,4*c); }
  ok=calloc(n,4); int nl=0, leaves[256];
  for(int i=0;i<n;i++){ ok[i]= deg[i]>=2; if(deg[i]==1) leaves[nl++]=i; }
  if(ea<0) ea=leaves[rnd()%nl]; if(eb<0){ do eb=leaves[rnd()%nl]; while(eb==ea); }
  ok[ea]=ok[eb]=1; target=0; for(int i=0;i<n;i++) target+=ok[i];
  vis=calloc(n,4); pos=malloc(4*n); p=malloc(4*(n+1));
  len=0; p[len]=ea; pos[ea]=len++; vis[ea]=1;
  long long rots=0, stall=0;
  while(len<target && rots<maxrot){
    int e=p[len-1], bu=-1, bs=1<<30;
    for(int j=0;j<deg[e];j++){ int u=adj[e][j]; if(!ok[u]||vis[u]) continue; if(u==eb && len<target-1) continue;
      int r=0; for(int t=0;t<deg[u];t++) if(ok[adj[u][t]]&&!vis[adj[u][t]]) r++;
      int sc=r*16+(rnd()&15); if(sc<bs){bs=sc;bu=u;} }
    if(bu>=0){ p[len]=bu; pos[bu]=len++; vis[bu]=1; stall=0; continue; }
    if(stall>200){ int ins=insertion_pass(); stall=0; if(ins) continue; }
    int cand[64], nc=0; for(int j=0;j<deg[e];j++){ int w=adj[e][j]; if(ok[w]&&vis[w]&&pos[w]<len-2) cand[nc++]=w; }
    if(!nc){ int ins=insertion_pass(); if(!ins) break; continue; }
    int w=-1;
    for(int t=0;t<nc && w<0;t++){ int c=cand[(t+rnd())%nc]; int ne=p[pos[c]+1];
      for(int j=0;j<deg[ne];j++){ int u=adj[ne][j]; if(ok[u]&&!vis[u]&&!(u==eb&&len<target-1)){ w=c; break; } } }
    if(w<0) w=cand[rnd()%nc];
    int i=pos[w];
    for(int a=i+1,b=len-1;a<b;a++,b--){ int t=p[a]; p[a]=p[b]; p[b]=t; pos[p[a]]=a; pos[p[b]]=b; }
    rots++; stall++;
  }
  if(len==target-1 && isadj(p[len-1],eb)){ p[len]=eb; pos[eb]=len++; }
  fprintf(stderr,"ends %lld %lld target %d reached %d rotations %lld\n",lab[ea],lab[eb],target,len,rots);
  FILE*o=fopen(argc>6?argv[6]:"path.txt","w"); for(int i=0;i<len;i++) fprintf(o,"%lld\n",lab[p[i]]); fclose(o);
  printf("%d\n",len); return 0; }
