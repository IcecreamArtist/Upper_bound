# include <bits/stdc++.h>
# define pb push_back
# define fi first
# define se second
# define ll long long
# define fo(i,a,b) for(int i=(a);i<=(b);++i)
using namespace std;
const int maxn = 1e5+10;
int n,q;
int a[maxn];
int npri[maxn];
int p[maxn],cnt = 0;
int nxt[maxn];
int pnxt[maxn][400];
int fpin[400];
int spin[400];
int r[maxn][20];
int main(){
    fo(i,1,maxn-1) npri[i] = 1;
    npri[2] = 1;
    fo(i,2,maxn-1){
        if(npri[i]) p[cnt++] = i;
        for(int j=0;j<cnt && p[j] * i < maxn;++j){
            npri[p[j] * i] = 0;
        }
    }
    scanf("%d%d",&n,&q);
    fo(i,1,n){
        scanf("%d",a+i);
    }
    fo(j,0,cnt-1)nxt[j] = pnxt[n+1][j] = n+1;
    for(int i=n;i;--i){
        for(int j=0;j<cnt;++j){
            if(a[i] % p[j] == 0){
                printf("G %d %d\n",i,j);
                pnxt[i][j] = nxt[j] = i;
            }
            else pnxt[i][j] = nxt[j];
        }
    }
    
    fo(i,1,n){
        fo(j,0,3){
            printf("%d ",pnxt[i][j]);
        }
        printf("\n");
    }

    fo(j,0,cnt-1) fpin[j] = pnxt[1][j],spin[j]=pnxt[fpin[j]][j];
    fo(l,1,n){
        r[l][0] = n+1;
        fo(j,0,cnt-1) r[l][0] = min(r[l][0],spin[j]-1);
        //printf("Get r[%d][0] = %d\n",l,r[l][0]);
        fo(j,0,cnt-1){
            if(a[l] % p[j] == 0){
                fpin[j] = spin[j];
                spin[j] = pnxt[spin[j]][j];
            }
        }
    }
    for(int k=0;k<20;++k) r[n+1][k] = n+1;
    for(int i=n;i;--i){
        for(int l=1;l<20;++l){
            r[i][l] = r[r[i][l-1]][l-1];
        }
    }
    while(q--){
        int x,y; scanf("%d%d",&x,&y);
        int mint = 0;
        while(x != y){
            int k = 0;
            while(k<20 && r[x][k] < y)++k;
            mint += 1<<k;
            x = r[x][k] + 1;
        }
        printf("%d\n",mint);
    }
    return 0;
}