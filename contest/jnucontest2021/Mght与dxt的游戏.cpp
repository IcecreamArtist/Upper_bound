//
// Created by Artis on 2021/5/11.
//

#include<bits/stdc++.h>
using namespace std;
const int maxn = 1e6+5;
char str[maxn];
int l[maxn],r[maxn],n,pos[maxn];

int dfs(int u){
    if(l[u]) l[u] += dfs((u-1+n)%n);
    return l[u];
}

int dfs2(int u){
    if(r[u]) r[u] += dfs2((u+1)%n);
    return r[u];
}

int main(){
    int t;scanf("%d",&t);
    while(t--){
        scanf("%d",&n);
        for(int i=0;i<n;++i) l[i]=r[i]=0;
        scanf("%s",str);
        int flg = 0;
        for(int i=1;i<n;++i) if(str[i]!=str[i-1]) flg=1;
        if(!flg) {
            cout<<n<<endl;
            continue;
        }
        // 可以简化：记录每个位置右边第一个不一样的字母的位置
        // 先弄出最后一个右边第一个不一样
        for(int i=0;i<n-1;++i){
            if(str[i]!=str[n-1]){
                pos[n-1]=i;
                if(str[pos[n-1]]>str[n-1]) r[n-1]=1;
                break;
            }
        }
        for(int i=n-2;i>=0;--i) {
            if(str[i]==str[i+1]) pos[i]=pos[i+1];
            else pos[i]=i+1;
            if(str[pos[i]]>str[i]) r[i]=1;
        }

        for(int i=0;i<n;++i){
            if(r[i]==0){
                // 当前不能往右。那么当前的右能往左。
                l[(i+1)%n] = 1;
            }
        }

        for(int i=n-1;i>=0;--i){
            // 如果上一个为0
            if(l[(i+1)%n]==0) dfs(i);
        }
        for(int i=0;i<n;++i){
            // 如果上一个为0
            if(r[(i-1+n)%n]==0) dfs2(i);
        }
        int ans = 0;
        for(int i=0;i<n;++i){
            if(l[i]%2==0&&r[i]%2==0) ans++;
        }
      //  for(int i=0;i<n;++i) cout<<l[i]<<" ";
      //  cout<<endl;
      //  for(int i=0;i<n;++i) cout<<r[i]<<" ";
      //  cout<<endl;
        printf("%d\n",ans);
    }
}
