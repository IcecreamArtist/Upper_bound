#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
ll f[100000];

void init(){
    f[1]=1;
    for(int i=2;i<=100000;++i) f[i]=f[i-1]+f[i-2];
}

ll quickp(ll a,int k){
    ll ans = 1;
    while(k){
        if(k&1) ans*=a;
        a*=a;
        k>>=1;
    }
    return ans;
}

int main(){
    int t;scanf("%d",&t);
    init();
    while(t--){
        ll ans = 0,tmp;
        int n,c,k;scanf("%d%d%d",&n,&c,&k);
        for(int i=0;i<=n;++i) cout<<f[i*c]<<" ";
        cout<<endl;
        for(int i=0;i<=n;++i) tmp=quickp(f[i*c],k),ans+=tmp,cout<<tmp<<" ";
        cout<<endl;
        cout<<ans<<endl;
    }
}
