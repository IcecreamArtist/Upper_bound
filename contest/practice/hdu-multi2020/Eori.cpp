# include <bits/stdc++.h>
# define ll long long
using namespace std;
const ll MOD = 1e9+9;
const ll sq5 = 383008016;
ll qpow(ll x,ll n){
    if(n>MOD) n%=(MOD-1);
    ll ans = 1;
    while(n){
        if(n&1) ans=ans*x%MOD;
        x = x*x%MOD;
        n>>=1;
    }
    return ans;
}
const int maxn = 1e5+10;
ll jc[maxn];
ll ijc[maxn];
ll N,C;
int K;
const int alpha = 691504013;
const int beta = 308495997;

void solve(){
    scanf("%lld%lld%d",&N,&C,&K);
    ll an = 1;
    ll bn = qpow(beta,K*C);
    ll ac = qpow(alpha,C);
    ll bc = qpow(qpow(beta,C),MOD-2);
    ll abn = qpow(an*bn%MOD,N+1);
    ll acN = qpow(ac,N+1);
    ll bcN = qpow(bc,N+1);
    int ans = 0;
    //489317468
    int res;
    for(int k=0;k<=K;++k){
        ll cnk = jc[K]*ijc[K-k]%MOD*ijc[k]%MOD;
        ll t=an*bn%MOD;
            res = (abn-1>=0?abn-1:abn-1+MOD)*qpow(t-1,MOD-2)%MOD*cnk%MOD;
            if((K-k)&1) (ans-res+MOD>=MOD?ans-=res:ans=ans-res+MOD);
            else (ans+res>=MOD?ans=ans+res-MOD:ans=ans+res);
            an = an*ac%MOD;
            bn = bn*bc%MOD;
            abn = abn * acN%MOD * bcN%MOD;
    }
    printf("%lld\n",qpow(qpow(sq5,K),MOD-2)*ans%MOD);
}
int main(){
    //freopen("test.in","r",stdin);
    //freopen("test.out","w",stdout);
    jc[0] = 1;
    for(int i=1;i<maxn;++i)jc[i]=1ll*jc[i-1]*i%MOD;
    ijc[maxn-1] = qpow(jc[maxn-1],MOD-2);
    for(int i=maxn-2;i>=0;--i)ijc[i] = 1ll*ijc[i+1]*(i+1)%MOD;
    
    
    int T; scanf("%d",&T);
    while(T--) solve();
    return 0;
}
