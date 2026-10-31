/* Puzzle 120, Q4 (writeup/120.md): residue-mod-3 test for an n x n full-torus primal square.
   The cells form a graph (grid neighbours plus mirror pairs (r,c)-(r,n-1-c), (r,c)-(n-1-r,c); with rivera=1 the
   mirror pairs between main-diagonal cells are dropped). In a primal square of 1..n^2 neighbours have opposite
   parity, so one checkerboard class holds the even numbers (orient chooses which). Colour each cell by its value
   mod 3: neighbours can not be 0-0 or 1-2 (sum > 3 divisible by 3), except one edge joining the cells of 1 (odd,
   colour 1) and 2 (even, colour 2), sum 3. The number of cells of each (parity, colour) is fixed by 1..n^2.
   The program searches all such colourings; "impossible" means no primal square exists for that n and orient.
   Build: gcc -O2 -o p120_mod3 p120_mod3.c      Usage: ./p120_mod3 n orient allow12 [rivera]
   ./p120_mod3 8 0 1 0 and ./p120_mod3 8 1 1 0 take about a second each. */
#include <stdio.h>
#include <stdlib.h>
int rivera = 0, n, N, nbr[400][8], nn[400], order[400], pos[400], col[400], par[400], need[2][3], have[2][3], allow12, used12;
long long nodes = 0, sols = 0;
void addedge(int a,int b){ for(int i=0;i<nn[a];i++) if(nbr[a][i]==b) return; nbr[a][nn[a]++]=b; nbr[b][nn[b]++]=a; }
int bt(int k){
  nodes++;
  if(k==N){ sols++; return 1; }
  int c = order[k], p = par[c];
  for(int x=0;x<3;x++){
    if(have[p][x] >= need[p][x]) continue;
    int bad12 = 0, ok = 1;
    for(int i=0;i<nn[c] && ok;i++){ int d = nbr[c][i]; if(pos[d] >= k) continue;
      int y = col[d];
      if(x==0 && y==0) ok = 0;
      else if((x==1 && y==2 && p==1 && par[d]==0) || (x==2 && y==1 && p==0 && par[d]==1)) bad12++;
      else if((x==1 && y==2) || (x==2 && y==1)) ok = 0; }
    if(!ok) continue;
    if(bad12 > 0 && !(allow12 && used12 == 0 && bad12 == 1)) continue;
    col[c] = x; have[p][x]++; used12 += bad12;
    /* remaining capacity check */
    int r = bt(k+1);
    have[p][x]--; used12 -= bad12;
    if(r) return r;
  }
  return 0;
}
int main(int argc,char**argv){
  n = atoi(argv[1]); int orient = atoi(argv[2]); allow12 = atoi(argv[3]); if(argc > 4) rivera = atoi(argv[4]); N = n*n;
  for(int r=0;r<n;r++) for(int c=0;c<n;c++){ int a=r*n+c;
    if(c+1<n) addedge(a,r*n+c+1); if(r+1<n) addedge(a,(r+1)*n+c);
    int dg = (r==c || r+c==n-1);
    if(c!=n-1-c && !(rivera && dg)) addedge(a,r*n+(n-1-c)); if(r!=n-1-r && !(rivera && dg)) addedge(a,(n-1-r)*n+c);
    par[a] = ((r+c)%2 == 0) ? orient : 1-orient;   /* par = parity of the number in the cell: 0 even, 1 odd */ }
  for(int v=1; v<=N; v++) need[v%2][v%3]++;
  int k=0; for(int i=0;i<n/2;i++) for(int j=0;j<n/2;j++){ int rr[2]={i,n-1-i}, cc[2]={j,n-1-j};
    for(int a=0;a<2;a++) for(int b=0;b<2;b++) order[k++]=rr[a]*n+cc[b]; }
  for(int i=0;i<N;i++) pos[order[i]]=i;
  int r = bt(0);
  printf("n=%d rivera=%d orient=%d allow12=%d: colouring %s (nodes %lld)\n", n, rivera, orient, allow12, r ? "EXISTS" : "impossible", nodes);
  if(r) for(int a=0;a<n;a++){ for(int b=0;b<n;b++) printf("%d", col[a*n+b]); printf("\n"); }
  return 0;
}
