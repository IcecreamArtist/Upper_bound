# include <bits/stdc++.h>
# define ll long long
using namespace std;
const int maxn = 5000010;
const int maxs = 5000010;
int mu[maxn]={0,1};
bool nprime[maxn];
int prime[maxn],num;
int pre[maxn];
void getprime(){
    for (int i=2;i<=maxs;++i){
        if (!nprime[i]){ prime[num++]=i; mu[i]=-1; }
        for (int j=0;j<num&&prime[j]*i<=maxs;++j){
            nprime[i*prime[j]] = true;
            if (i%prime[j] == 0){
                mu[i*prime[j]] = 0;
                break;
            }
            else mu[i*prime[j]] = -mu[i];
        }
    }
}

const int MOD = 1e9+7;
inline void adm(int& x,int y){x+=y-(x+y>=MOD?MOD:0);}
int a,k;
int main(){
    getprime();   
    for (int i=1;i<=5000000;++i) adm(pre[i],pre[i-1] + mu[i]);
    cin >> a >> k;
    ll ans = 0;
    ll powSigma = 1;
    for(int d=1;d<=k;++d){
        (powSigma=powSigma*a)%=MOD;
        (ans += powSigma*pre[k/d]) %= MOD;
    }
    cout << ans << endl;
    return 0;
}
