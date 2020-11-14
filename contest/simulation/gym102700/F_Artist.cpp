#include<bits/stdc++.h>
using namespace std;
int n,m,a,b,k;
const int INF = 0x3f3f3f3f;
/*
 * 题意：有一张有向图。两个起点。要到同一个点之后返回。这条路径可以有k段免费，整个过程要最短路。
 * 思路：
 * dijkstra分别从两个起点到所有的点。dp记录到这个点，花费多少次免费票，最少的cost
 * 反向建图，再来一次dijkstra。
 * 暴力枚举所有的相遇点，得到结果。
 */
struct edge{
    int v,cost;
};

struct node{
    int u,cost,used;
    bool operator < (const node b)const{
        return cost>b.cost;
    }
};

const int N = 1e5+6;
int dp1[N][11],dp2[N][11],dp3[N][11],dp4[N][11];
vector<edge>G1[N];
vector<edge>G2[N];

void dijkstra(vector<edge>G[],int st,int dp[][11]){
    priority_queue<node> pq;
    for(int i=1;i<=n;++i) for(int j=0;j<=k;++j) dp[i][j] = INF;
    dp[st][0]=0;
    pq.push(node{st,0,0});
    while(!pq.empty()){
        node cur=pq.top();pq.pop();
        for(auto i:G[cur.u]){
            if(dp[i.v][cur.used]>cur.cost+i.cost){
                dp[i.v][cur.used] = cur.cost+i.cost;
                pq.push(node{i.v,dp[i.v][cur.used],cur.used});
            }
            if(cur.used<k&&dp[i.v][cur.used+1]>cur.cost){
                dp[i.v][cur.used+1] = cur.cost;
                pq.push(node{i.v,dp[i.v][cur.used+1],cur.used+1});
            }
        }
    }
    // 能少用票固然少用票。票不必用完
    for(int i=1;i<=n;++i) for(int j=1;j<=k;++j) dp[i][j] = min(dp[i][j],dp[i][j-1]);
}

int main(){
    scanf("%d%d",&n,&m);
    scanf("%d%d%d",&a,&b,&k);
    a++,b++;
    for(int i=0;i<m;++i){
        int u,v,c;scanf("%d%d%d",&u,&v,&c);
        u++,v++;
        G1[u].push_back(edge{v,c});
        G2[v].push_back(edge{u,c});
    }
    dijkstra(G1,a,dp1);
    dijkstra(G1,b,dp2);
    dijkstra(G2,a,dp3);
    dijkstra(G2,b,dp4);
    int ans = INF,pos = -1;
    for(int i=1;i<=n;++i){
        // 枚举相遇点
        if(i==a||i==b) continue;
        int ans1 = INF,ans2 = INF;
        // 从a出发到当前点的最短距离
        // 以及从b出发到当前点的最短距离
        for(int j=0;j<=k;++j){
            // 枚举使用的票数
            ans1=min(ans1,dp1[i][j]+dp3[i][k-j]);
            ans2=min(ans2,dp2[i][j]+dp4[i][k-j]);
        }
        if(ans1+ans2<ans){
            ans = ans1 + ans2;
            pos = i;
        }
    }
    if(pos == -1) puts(">:(");
    else printf("%d %d\n",pos-1,ans);
}