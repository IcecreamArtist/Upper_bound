# include <bits/stdc++.h>
using namespace std;
const int MOD = 1e9+9;
int main(){
    for(int i=1;i<MOD;++i){
        if(1ll*i*i%MOD == 5){
            printf("%d\n",i);
        }
    }
    return 0;
}
