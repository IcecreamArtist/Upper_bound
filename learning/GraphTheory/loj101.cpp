# include <bits/stdc++.h>
using namespace std;
const int maxn = 5e3 + 10;
struct edge{
	int v,nxt;
	int flow,cap;
} e[maxn<<1]; int cnt=0;
int g[maxn];
inline void addedge(int u,int v,int flow){
	e[++cnt].v=v;
	e[cnt].flow=flow,e[cnt].cap=flow;
	e[cnt].nxt=g[u];
	g[u]=cnt;
}
int n,m,s,t;
int dfs(int x,int flow){
	
}
int Dinic(){
	
}
int main(){
	scanf("%d%d%d%d",&n,&m,&s,&t);
	for(int i=1;i<=m;++i){
		int u,v,c; scanf("%d%d%d",&u,&v,&c);
		addedge(u,v,c);
	}
	printf("%d\n",Dinic(s,t))
	return 0;
}