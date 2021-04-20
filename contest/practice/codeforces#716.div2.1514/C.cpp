//
# include <bits/stdc++.h>
# define fo(i,a,b) for(int i=(a);i<=(b);++i)
# define ll long long
# define pb push_back
# define fi first
# define se second
using namespace std;
ll gcd(ll x,ll y){
    return y?gcd(y,x%y):x;
}
int main(){
    int n; scanf("%d",&n);
    ll mul = 1;
    int cnt = 0;
    fo(i,1,n-1){
        if(gcd(i,n) == 1){
            mul = mul * i % n;
            ++cnt;
        }
    }
    if(mul == 1){
        printf("%d\n",cnt);
        fo(i,1,n-1){
            if(gcd(i,n)==1)printf("%d ",i);
        }
        printf("\n");
        return 0;
    }
    printf("%d\n",cnt-1);
    fo(i,1,n-1){
        if(mul != i && gcd(i,n) == 1){
            printf("%d ",i);
        }
    }
    printf("\n");
    return 0;
}