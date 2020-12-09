# 区间DP
## 应用：石子合并
```c++
#include<bits/stdc++.h>
using namespace std;
int a[104];
int dp[104][104],dp2[103][104];
int sum[104];
const int inf=1e9+7;
int main(){
	int n;
	while(~scanf("%d",&n)){
		for(int i=1;i<=n;i++){
			scanf("%d",&a[i]);
			sum[i]=sum[i-1]+a[i];//前缀和
			dp[i][i]=0;dp2[i][i]=0;
		} 
		for(int k=1;k<=n-1;k++){//枚举长度
			for(int i=1;i<n;i++){
				int j=i+k;
				if(j>n) break;
				dp[i][j]=inf;
				dp2[i][j]=0;
				int tmp=sum[j]-sum[i-1];
				for(int m=i;m<j;m++){
					dp[i][j]=min(dp[i][m]+dp[m+1][j]+tmp,dp[i][j]);
					dp2[i][j]=max(dp2[i][m]+dp2[m+1][j]+tmp,dp2[i][j]);
				}
			}
		}
		printf("%d %d\n",dp[1][n],dp2[1][n]);
	}
}
```
## 应用：环形石子合并
```c++
#include<bits/stdc++.h>
using namespace std;
int stone[102];
int sum[202];
int dp[202][202];
int dp2[202][202];
const int inf=1e9;
int main(){
    int n;
    scanf("%d",&n);
    for(int i=1;i<=n;i++){
        scanf("%d",&stone[i]);
        sum[i]=sum[i-1]+stone[i];
    }
    for(int i=n+1;i<=2*n;i++)sum[i]=sum[i-1]+stone[i-n];
    for(int len=1;len<=n-1;len++){
        for(int i=1;i<=2*n-1;i++){
            int j=i+len;
            if(j>2*n)break;
            dp2[i][j]=inf;
            for(int k=i;k<j;k++){
        dp[i][j]=max(dp[i][j],dp[i][k]+dp[k+1][j]+sum[j]-sum[i-1]);
        dp2[i][j]=min(dp2[i][j],dp2[i][k]+dp2[k+1][j]+sum[j]-sum[i-1]);
            }
        }
    }
    int maxn=0,minn=inf;
    //最后寻找最好的段落！
    for(int i=1;i<=n;i++){
        maxn=max(dp[i][i+n-1],maxn);
        minn=min(dp2[i][i+n-1],minn);
    }
    printf("%d\n%d\n",minn,maxn);
} 
```
# 树形DP
```c++
ll dp[maxn][3];
void dfs(int u,int fa){
	dp[u][0]=dp[u][1]=1;
	for(int i=head[u];i;i=e[i].next){
		int v=e[i].to;
		if(v==fa) continue;
		dfs(v,u);
		dp[u][1]=(dp[u][1]*(dp[v][0]+dp[v][1]))%mod;
		dp[u][0]=(dp[u][0]+dp[v][0]+dp[v][1]-1)%mod;
	}
}
```
## 应用：换根DP
dp[i]=dp[u]-pt[i]+(n-pt[i])
```c++
void dfs(int u,int fa){
    for(auto i:G[u]) if(i-fa) dfs(i,u),dis[u]+=dis[i]+pt[i],pt[u]+=pt[i];
    pt[u]+=1;
}
void DP(int u,int fa){
    for(auto i:G[u]) if(i-fa) dp[i]=dp[u]-pt[i]+(n-pt[i]),DP(i,u);
}
int main(){
    dfs(1,-1);
    dp[1]=dis[1];
    DP(1,-1);
    int ans=inf;
    FOR(i,1,n) ans=min(ans,dp[i]);
    printf("%d\n",ans);
}
```
# 数位DP
```c++
#include<bits/stdc++.h>
using namespace std;
const int N = 30;
typedef long long ll;
ll dp[N][3]; // 0：前面是其他,1：前面是4,2：前面是48
int digit[N];

int dfs(int pos,int pre,bool ismax){
    if(pos==0) return 1;
    int res = 0;
    if(!ismax&&dp[pos][pre]!=-1) return dp[pos][pre];
    int mx = ismax?digit[pos]:9;
    for(int i=0;i<=mx;++i){
        if(i==6&&pre==2) continue;
        if(i==8&&pre==1) res += dfs(pos-1,2,ismax&&(i==mx));
        else if(i==4) res += dfs(pos-1,1,ismax&&(i==mx));
        else res += dfs(pos-1,0,ismax&&(i==mx));
    }
    if(!ismax) dp[pos][pre] = res;
    return res;
}

int solve(int n){
    memset(digit,0,sizeof(digit));
    int len = 0;
    while (n){
        digit[++len] = n % 10;
        n/=10;
    }
    return dfs(len,0,1);
}

int main(){
    ll n;scanf("%lld",&n);
    memset(dp,-1,sizeof(dp));
    ll ans = solve(n);   // 1~n中没有486的数字的个数
    printf("%lld\n",n-ans+1);    // n-ans，即有的个数
}
```
# DP优化
## 单调队列/单调栈优化