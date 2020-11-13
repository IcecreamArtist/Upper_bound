# include <bits/stdc++.h>
# define pb push_back
# define mp make_pair
# define fi first
# define se second
# define INF 0x3f3f3f3f
# define cl(x) memset(x,0,sizeof x)
using namespace std;
inline int read(){
	int x=0; char c=getchar();
	for(;c<'0'||c>'0';c=getchar());
	for(;c>='0'&&c<='9';c=getchar())x=x*10+c-'0';
	return x;
}
/*
Use SPFA the shortest-path algo to solve this
For an edge (u,v):
dis[v][k+1] <- dis[u][k] if use the ticket once
dis[v][k] <- dis[u][k]+w of edge if not use the ticket
*/
int n,m;
vector<int> g1[maxn];
vector<int> w1[maxn];

vector<int> g2[maxn];
vector<int> w2[maxn];
int a,b,k;	//a->,b->,k tickets each
int dis1[maxn][maxk];
int dis2[maxn][maxk];
int dis3[maxn][maxk];
int dis3[maxn][maxk];
queue<pair<int,int> >Q;
bool inq[maxn][maxk];
void SPFA(vector<int> g[],vector<int> w[],int dis[maxn][maxk],int start){
	cl(inq);
	for(int i=1;i<=n;++i)for(int j=0;j<k;++j) dis[i][j] = INF;
	dis[start][0] = 0;
	Q.push(mp(start,0));
	while(!Q.empty()){
		int u = Q.front().fi,j = Q.front().se;Q.pop();
		inq[u][j] = 0;
		for(int i=0;i<g[u].size;++i){
			int v=g[u][i];
			if(dis[v][j] > dis[u][j]+w[u][i]){
				dis[v][j] = dis[u][j]+w[u][i];
				if(!inq[v][j]){
					Q.push(mp(v,j));
					inq[v][j] = 1;
				}
			}
			if(j<k&&dis[v][j+1]>dis[u][j]){
				dis[v][j+1]=dis[u][j];
				if(!inq[v][j+1]){
					Q.push(mp(v,j+1));
					inq[v][j+1]=1;
				}
			}
		}
	}
}
int main(){
	n=read(),m=read();
	a=read(),b=read();
	k=read();
	for(int i=1;i<=m;++i){
		int u,v,w;
		u=read(),v=read(),w=read();
		g1[u].pb(v);
		g1[u].pb(w);
		
		g2[v].pb(u);
		w2[v].pb(w);
	}
	
	return 0;
}