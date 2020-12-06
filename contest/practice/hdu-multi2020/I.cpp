# include <bits/stdc++.h>
# define fo(i,a,b) for(int i=(a);i<=(b);++i)
using namespace std;
const int maxn = 5e4+10;
int n;
struct robot{
    int p,a;
} r[maxn];
map<pair<int,int>,int> S;
bool cmp(robot x,robot y){
    return x.p == y.p ? x.a > y.a : x.p < y.p;
}
int sta[maxn],pos;
inline bool judge(int a,int b,int c){
    if(r[a].a == r[b].a){
        return 0;
    }
    else{
        if(r[b].a == r[c].a){
            //Actually impossible
            return 0;
        }
        else{
            if(r[b].a < r[a].a) return 0;
            return 1ll*(r[a].p - r[b].p) * (r[c].a - r[b].a) < 1ll*(r[b].p - r[c].p) * (r[b].a - r[a].a);
        }
    }
}
void solve(){
    scanf("%d",&n);
    fo(i,1,n) scanf("%d%d",&r[i].p,&r[i].a),S[make_pair(r[i].p,r[i].a)]++;
    sort(r+1,r+n+1,cmp);
    
    pos = 0;
    sta[pos++] = 1;
    fo(i,2,n){
        if(r[i].p == r[i-1].p) continue;
        while(pos > 1 && !judge(i,sta[pos-1],sta[pos-2])){--pos; }
        if(pos == 1 && r[i].a >= r[sta[pos-1]].a){
            --pos;
        }
        sta[pos++] = i;
    }
    int cnt = 0;
    for(int i=0;i<pos;++i){
        if(S[make_pair(r[sta[i]].p,r[sta[i]].a)]==1) ++cnt;
    }
    printf("%d\n",cnt);
    S.clear();
}
int main(){
    int T; scanf("%d",&T);
    while(T--) solve();
    return 0;
}
