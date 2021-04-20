# include <bits/stdc++.h>
# define fo(i,a,b) for(int i=(a);i<=(b);++i)
# define ll long long
# define pb push_back
# define fi first
# define se second
using namespace std;
const int maxn = 1e4+10;
int n;
int a[maxn];
void solve(){
    scanf("%d",&n);
    bool tag = 0;
    fo(i,1,n){
        scanf("%d",a+i);
        int x = sqrt(a[i]);
        if(x*x != a[i]){
            tag = 1;
        }
    }
    if(tag) puts("YES");
    else puts("NO");
}
int main(){
    int T; scanf("%d",&T);
    while(T--)solve();
    return 0;
}