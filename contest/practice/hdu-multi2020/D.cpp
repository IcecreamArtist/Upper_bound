#include<bits/stdc++.h>
using namespace std;
#define FOR(i,l,r) for(int i=l;i<=r;i++)
typedef long long ll;

int main(){
    
    int t;scanf("%d",&t);
    while(t--){
        int n;scanf("%d",&n);
        if(n==1) puts("26");
        else if(n==2) puts("676");
        else if(n==3) printf("%d\n",26*26*26);
        else printf("%d\n",26*25*24);
    }
    return 0;
}
