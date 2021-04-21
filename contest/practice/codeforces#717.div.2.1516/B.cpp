# include <bits/stdc++.h>
# define pb push_back
# define fi first
# define se second
# define ll long long
# define fo(i,a,b) for(int i=(a);i<=(b);++i)
using namespace std;
const int maxn = 2010;
int a[maxn];
int xa[maxn];
int n;
void solve(){
    scanf("%d",&n);
    fo(i,1,n) scanf("%d",a+i);
    fo(i,1,n) xa[i] = xa[i-1]  ^ a[i];
    if(xa[n] == 0){
        printf("YES\n");
        return;
    }
    int p=-1;
    fo(i,1,n-1){
        if(xa[i] == xa[n]){
            p=i;
            break;
        }
    }
    if(p==-1){
        printf("NO\n");
        return;
    }
    fo(i,p+1,n-1){
        if(xa[i] == 0){
            printf("YES\n");
            return;
        }
    }
    printf("NO\n");
    return;
}
int main(){
    int T; scanf("%d",&T);
    while(T--)solve();
    return 0;
}