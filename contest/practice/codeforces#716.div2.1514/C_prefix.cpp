//
// Created by Lenovo on 2021/4/19.
//
#include<bits/stdc++.h>
using namespace std;

int gcd(int a,int b){
    return b?gcd(b,a%b):a;
}
typedef long long ll;
int main(){
    int n;cin>>n;
    ll tmp=1;int cnt = 1,ans = 1;
    for(int i=2;i<=n-1;++i){
        if(gcd(i,n)==1) {
            tmp=tmp*i%n;
            cnt++;
            if(tmp%n==1) ans = cnt;
        }
    }
    cout<<ans<<endl;
    for(int i=1;i<=n-1;++i){
        if(gcd(i,n)==1) {
            cout<<i<<" ";
            ans--;
        }
        if(ans==0) break;
    }
}
