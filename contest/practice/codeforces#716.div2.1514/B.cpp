# include <bits/stdc++.h>
# define fo(i,a,b) for(int i=(a);i<=(b);++i)
# define ll long long
# define pb push_back
# define fi first
# define se second
using namespace std;
const ll MOD = 1e9+7;
ll qpow(ll x,ll n){
    if (n==0)return 1;
    if (n==1)return x;
    ll d = qpow(x,n>>1);
    if(n&1)return d*d%MOD*x%MOD;
    else return d*d%MOD;
}
void solve(){
    int n,k; scanf("%d%d",&n,&k);
    printf("%lld\n",qpow(n,k));
}
int main(){
    int T; scanf("%d",&T);
    while(T--)solve();
    return 0;
}