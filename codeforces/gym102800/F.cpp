#include<bits/stdc++.h>
#define For(i,l,r) for(int i=l;i<=r;i++)
#define ll long long
using namespace std;
const int maxn = 1e5+10;
inline int lowbit(int x){return x&(-x);}
int n;
ll T[maxn<<1];
void add(int x){
    for(;x<=100000;x+=lowbit(x))T[x]++;
}
ll query(int x){
    ll s = 0;
    for(;x>0;x-=lowbit(x))s+=T[x];
    return s;
}
int a[maxn];
ll pre[maxn];
ll temp[maxn];
void solve(){
    memset(T,0,sizeof T);
    scanf("%d",&n);
    for(int i=1;i<=n;++i){ scanf("%d",a+i); ++ a[i]; }
    ll sum = 0;
    for(int i=n;i>0;--i){   //Back, < a[i]
        pre[i] = query(a[i]-1);
        sum += pre[i];
        add(a[i]);
    }
    int m; scanf("%d",&m);
    ll ans = sum;
    for(int i=1;i<=n;i++){
        temp[i] = pre[i];
    }
    while(m--){
        int p,q; scanf("%d%d",&p,&q);
        int num1=0, num2=0;
        for(int i=q-1;i>=p+1;i--) {
            if (a[q]>a[i]) {
                num1++;
            }
            else if(a[q]<a[i]) {
                pre[i]--;
            }
        }
        for(int i=p+1;i<=q-1;i++){
            if(a[p]>a[i]) {
                num2++;
            }
            else if(a[p]<a[i]){
                pre[i]++;
            }
        }
        if(a[p]>a[q]){
            pre[p]=pre[p]-num2-1;
            pre[q]=pre[q]+num1;
        }
        else if(a[p]<a[q]){
            pre[p] = pre[p] - num2;
            pre[q] = pre[q] + num1 + 1;
        }
        else {
            pre[p] -= num2;
            pre[q] +=num1;
        }
        swap(a[p],a[q]);
        swap(pre[p],pre[q]);

        for(int i=p;i<=q;i++){
            ans = ans - temp[i] + pre[i];
        }
        for(int i=p;i<=q;i++){
            temp[i]=pre[i];
        }
        sum = min(sum,ans);
    }
    printf("%lld\n",sum);
}
int main(){
	//freopen("F.in","r",stdin);
	//freopen("a.out","w",stdout);
    int t;
    scanf("%d",&t);
    while(t--){
        solve();
    }
}
