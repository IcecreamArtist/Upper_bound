// 集卡有奖

#include<bits/stdc++.h>
#define lowbit(x) ((x)&-(x))
using namespace std;
typedef long long ll;
const int maxn = 1e5+5;
int a[maxn],c[maxn],L[maxn],ori[maxn];
ll dp[maxn];
int n,m;

int discretize(int n){
    for(int i=0;i<n;++i) c[i] = a[i];
    sort(c,c+n);
    int maxx = 0;
    for(int i=0;i<n;++i){
        L[i] = lower_bound(c,c+n,a[i])-c+1;
        maxx = max(maxx,L[i]);
        ori[L[i]] = a[i];
    }
    return maxx;
}

ll tree[maxn];

void add(int x,ll d){
    while(x<=n){
        tree[x] = max(tree[x],d);
        x += lowbit(x);
    }
}

ll query(int x){
    ll maxx = 0;
    while(x>0){
        maxx = max(maxx,tree[x]);
        x -= lowbit(x);
    }
    return maxx;
}

int main(){
    scanf("%d%d",&n,&m);
    for(int i=0;i<n;++i){
        scanf("%d",&a[i]);
    }
    int maxx = discretize(n);
    // 原始版本
    /*
    // 收集了i张卡片
    for(int i=1;i<=m;++i){
        // 到第j张卡片
        for(int j=1;j<=n;++j){
            for(int k=1;k<=j;++k){
                if(a[k]>=a[j]) continue;
                // 以第j张卡片结尾
                // 上一轮中，序号比j小，值比a[j]小，中的最大值。
                dp[i][j] = max(dp[i][j], dp[i-1][k]+a[j]);
            }
        }
    }
    */
    // BIT
    // dp[j]以第j张卡片结尾。
    // BIT 小于x(离散化后的a[j])，（存的是dp值）dp值最大

    for(int i=1;i<=m;++i){
        // 卡堆的卡片数量
        for(int j=0;j<n;++j){
            // 轮到第几张卡
            // 在上一个卡堆卡片数量的基础上，到这张卡之前的树状数组中查询
            ll lastdp = dp[j];
            ll tmp = query(L[j]-1);
            if(i!=1 && tmp==0) dp[j] = 0;
            // 有，续；或当前是第一轮，放入
            else dp[j] = tmp + a[j];
            add(L[j],lastdp); // 加入上一轮，这张卡入树状数组
        }
        for(int j=0;j<=n;++j) tree[j] = 0; //  清空一代的BIT
    }


    ll maxxx = 0;
    for(int i=0;i<n;++i) maxxx = max(dp[i],maxxx);
    if(maxxx) printf("%lld\n",maxxx);
    else puts("-1");
}
