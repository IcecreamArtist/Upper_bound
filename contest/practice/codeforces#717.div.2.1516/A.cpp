# include <bits/stdc++.h>
# define pb push_back
# define fi first
# define se second
# define ll long long
# define fo(i,a,b) for(int i=(a);i<=(b);++i)
using namespace std;
const int maxn = 110;
int a[maxn];
int n,k;
void solve(){
    scanf("%d%d",&n,&k);
    fo(i,1,n) scanf("%d",a+i);
    int sum = 0;
    int rp = n;
    fo(i,1,n){
        if(k < a[i]){
            a[i] -= k;
            sum += k;
            break;
        }
        sum += a[i];
        k -= a[i];
        a[i] = 0;
    }
    a[n] += sum;
    fo(i,1,n) printf("%d ",a[i]);
    printf("\n");
}
int main(){
    int T; scanf("%d",&T);
    while(T--)solve();
    return 0;
}