//
// Created by Artis on 2021/5/11.
//
#include<bits/stdc++.h>
using namespace std;
const int maxn = 1e5+5;
typedef long long ll;
int a[maxn];
ll b[maxn];
priority_queue<ll,vector<ll>,greater<ll> > pq; // 放的是设置为正号，到当前之前的数

int main(){
    int t;scanf("%d",&t);
    while(t--){
        int n;scanf("%d",&n);
        for(int i=1;i<=n;++i) {
            scanf("%lld",&a[i]);
            b[i+1]=b[i]+1ll*a[i];
        }
        n++;
        while(!pq.empty()) pq.pop();
        ll ans = 0;
        for(int i=1;i<n;i+=2){
            if(pq.size()==0){
                ans += 1ll*b[i+1]-b[i];
                pq.push(b[i+1]);
                continue;
            }
            ll ans1 = 1ll*b[i+1]-b[i]; // -+
            ll ans2 = 1ll*b[i+1]+b[i]-2ll*pq.top(); // ++
            if(ans1>ans2){
                pq.push(b[i+1]);
                ans += ans1;
            }else{
                pq.pop();
                pq.push(b[i+1]);
                pq.push(b[i]);
                ans += ans2;
            }
        }
        printf("%lld\n",ans);
    }
}
