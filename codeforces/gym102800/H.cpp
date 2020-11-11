# include <bits/stdc++.h>
# define ll long long
using namespace std;
int n,m,k;
int mu[maxs];
ll ans[maxs];
bool np[maxs];
int p[maxs],cnt;
//mul[x] = the nums of multiples of x.
int a[maxs],mul[maxs],c[maxs];
void solve(){
	scanf("%d%d%d",&n,&m,&k);
	int mx = 0;
	for(int i=1;i<=n;++i) scanf("%d",a+i),mx=max(mx,a[i]),c[a[i]]++;
	for(int i=1;i<=mx;++i){
		for(int j=i;j<=mx;j+=i)if(c[j])mul[i] += c[j];
	}
	
	for(int x=1;x<=n && x <= mx;++x){
		for(int d=1;(1ll*d*x)<=mx;++d){
			ans[x] += 1ll*mul[d*x] * mul[d*x] * mu[d];
		}
	}
	
	for(int i=1;i<=k;++i){
		int x; scanf("%d",&x);
		printf("%d\n",ans[x]);
	}
}
void getmu(){
	np[1]=np[2]=1;
	for(int i=2;i<=100000;++i){
		if(!np[i])p[++cnt]=i,mu[i]=-1;
		for(int j=1;j<=cnt&&p[j]*i<=mxas;++j){
			np[i*p[j]]=1;
			if(i % p[j] == 0){mu[i*p[j]]=0;break;}
			else mu[i*p[j]]=-mu[i];
		}
	}
}
int main(){
	getmu();
	int T; scnaf("%d",&T);
	while(T--) solve();
	return 0;
}