# include <bits/stdc++.h>
# define pb push_back
# define fi first
# define se second
# define ll long long
# define fo(i,a,b) for(int i=(a);i<=(b);++i)
using namespace std;
const int maxn = 110;
int n;
int a[maxn];
bitset<maxn*2000> pa;
int main(){
    scanf("%d",&n);
    pa[maxn*1000] = 1;
    fo(i,1,n){
        scanf("%d",a+i);
        pa = (pa << a[i]) | (pa >> a[i]);
    }
    
    if(!pa[maxn*1000]){
        printf("0\n");
        return 0;
    }
    while(1){
        fo(i,1,n){
            if(a[i] &1){
                printf("1\n%d\n",i);
                return 0;
            }
            a[i] >>= 1;
        }
    }
    return 0;
}