# include <bits/stdc++.h>
# define pb push_back
# define mp make_pair
# define fi first
# define se second
using namespace std;
const int maxn = 1e4+10;
const int maxm = 4e4+10;
const int maxk = 15;
const int INF = 1e9+10;
vector<int> g1[maxn];
vector<int> w1[maxn];
vector<int> g2[maxn];
vector<int> w2[maxn];
int n,m;
int a,b,k;
int dis1[maxn][maxk];
int dis2[maxn][maxk];
int dis3[maxn][maxk];
int dis4[maxn][maxk];
int disA[maxn],disB[maxn];

void Dij(vector<int> g[],vector<int> w[],int d[][maxk],int s){
    priority_queue<pair<int,pair<int,int> >,vector<pair<int,pair<int,int> > >,greater<pair<int,pair<int,int> > > > Q;
    for(int i=1;i<=n;++i) for(int j=0;j<=k;++j) d[i][j] = INF;
    d[s][0] = 0;
    Q.push(mp(d[s][0],mp(s,0)));
    while(!Q.empty()){
        int u = Q.top().se.fi,j=Q.top().se.se; Q.pop();
        for(int i=0;i<g[u].size();++i){
            int v=g[u][i];
            if(d[v][j] > d[u][j] + w[u][i]){
                d[v][j] = d[u][j] + w[u][i];
                Q.push(mp(d[v][j],mp(v,j)));
            }
            if(j<k&&d[v][j+1]>d[u][j]){
                d[v][j+1]=d[u][j];
                Q.push(mp(d[v][j+1],mp(v,j+1)));
            }
        }
    }
    for(int i=1;i<=n;++i)for(int j=1;j<=k;++j) d[i][j] = min(d[i][j],d[i][j-1]);
}
int main(){
//	freopen("F.in","r",stdin);
    scanf("%d%d",&n,&m);
    scanf("%d%d%d",&a,&b,&k);
	++a; ++b;
    for(int i=1;i<=m;++i){
        int u,v,w; scanf("%d%d%d",&u,&v,&w);
        ++u; ++v;
        g1[u].pb(v); w1[u].pb(w);
        g2[v].pb(u); w2[v].pb(w);
    }

    Dij(g1,w1,dis1,a);
    Dij(g1,w1,dis2,b);
    Dij(g2,w2,dis3,a);
    Dij(g2,w2,dis4,b);

    for(int i=1;i<=n;++i){
		disA[i] = dis1[i][0]+dis3[i][k];
		for(int j=1;j<=k;++j){
			disA[i] = min(disA[i],dis1[i][j]+dis3[i][k-j]);
		}
	}
	
	for(int i=1;i<=n;++i){
		disB[i] = dis2[i][0]+dis4[i][k];
		for(int j=1;j<=k;++j){
			disB[i] = min(disB[i],dis2[i][j]+dis4[i][k-j]);
		}
	}
	
	int pos = 1;
	while(((pos==a||pos==b)||(disA[pos]>INF||disB[pos]>INF)) && pos<=n)++pos;
	for(int i=pos;i<=n;++i){
		if(i==a||i==b||(disA[i]>INF||disB[pos]>INF))continue;
		//printf("For %d,res = %d=A:%d + B:%d\n",i,disA[i]+disB[i],disA[i],disB[i]);
		if(disA[pos]+disB[pos] > disA[i] + disB[i]){
			pos = i;
		}
	}
	if(disA[pos]>INF || disB[pos]>INF || pos > n) printf(">:(\n");
	else printf("%d %d\n",pos-1,disA[pos]+disB[pos]);
    return 0;
}

/*
4 5
0 1 2
0 2 2
1 2 2
2 3 2
3 0 2
3 1 2
*/
