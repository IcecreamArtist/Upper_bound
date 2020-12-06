# include <bits/stdc++.h>
# define ll long long
using namespace std;
const ll MOD = 1e9+9;
const ll sq5 = 383008016;
inline ll ps(ll x){return x>=MOD?x-MOD:x;}
ll qpow(ll x,ll n,ll d=1){
    if(n>MOD)n%=(MOD-1);
//    if(n==0)return 1;
//    if(n==1)return x;
//    d=qpow(x,n>>1);
//    if(n&1){
//        return x*d%MOD*d%MOD;
//    }
//    else{
//        return d*d%MOD;
//    }
    ll ans = 1;
    x%=MOD;
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
ll alpha,beta;
void solve(){
    scanf("%lld%lld%d",&N,&C,&K);
    ll ans = 0;
    ll an = 1;
    ll bn = qpow(beta,K*C);
    ll ac = qpow(alpha,C);
  //  ll abn = qpow(an*bn%MOD,N+1);
    ll inv_bc = qpow(qpow(beta,C),MOD-2);
    for(int k=0;k<=K;++k){
        ll cnk = jc[K]*ijc[K-k]%MOD*ijc[k]%MOD;
        ll res = (abn-1+MOD)%MOD*qpow(ps(an*bn%MOD-1+MOD),MOD-2)%MOD*cnk%MOD;
        if((K-k)&1) ans=ps(ans-res+MOD);
        else ans=ps(ans + res);
        an = an*ac%MOD;
        bn = bn*inv_bc%MOD;
     //   abn = abn*ac%MOD*inv_bc%MOD;
    }
    cout << ans*qpow(qpow(sq5,K),MOD-2)%MOD << '\n';
}
int main(){
    int T; scanf("%d",&T);
    alpha = (1+sq5) %MOD * qpow(2,MOD-2)%MOD;
    beta = (1-sq5+MOD)%MOD * qpow(2,MOD-2)%MOD;
    //for()
    jc[0] = 1;
    for(int i=1;i<maxn;++i)jc[i]=1ll*jc[i-1]*i%MOD;
    ijc[maxn-1] = qpow(jc[maxn-1],MOD-2);
    for(int i=maxn-2;i>=0;--i)ijc[i] = 1ll*ijc[i+1]*(i+1)%MOD;
    
    while(T--) solve();
    return 0;
}
