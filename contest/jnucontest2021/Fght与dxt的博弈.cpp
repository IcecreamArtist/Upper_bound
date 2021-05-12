//
// Created by Artis on 2021/5/10.
//

/*
 *
 * 第一种石子：有a组。每组两堆。每次在一组的一堆中拿x个或两堆同时拿x个。
 * 第一种石子是a组威佐夫博弈。把结果异或起来即可。
 *
 * 第二种石子：有b堆。每次从一堆中拿走1~i个石子。（b<=100000）。
 * 巴什博弈
 *
 * 第三种石子：有c堆。每次从一堆中拿走1~(当前石子数-i)个石子。
 * 把石子总数转化为max(n-i,0)。nim
 *
 * 第四种石子：有d堆。数量都不同。每次能从第i堆里拿x个，且该堆排名不变。
 * 也就是说这些堆的排名自始至终不能改变。
 * 排序，问题转化为前后两两差的nim。
 */

#include <bits/stdc++.h>
using namespace std;
const int maxn = 104;
int sg[maxn][maxn];
int s[3*maxn]; // 后继的sg的集合
int s4[100006];
int main(){
   // freopen("21.in","r",stdin);
    int t;scanf("%d",&t);
    double gold = (1+sqrt(5))/2;
    for(int i=0;i<maxn;++i){
        for(int j=0;j<maxn;++j){
            // (i,j)
            int a = min(i,j),b=max(i,j);
            double kk=(double)(b-a);
            int test=(int)(kk*gold);
            if(test==a) {
                sg[i][j] = 0;
              //  cout<<i<<" "<<j<<endl;
                continue;
            }
            memset(s,0,sizeof(s));
            for(int k=0;k<i;++k) s[sg[k][j]]=1;
            for(int k=0;k<j;++k) s[sg[i][k]]=1;
            for(int k=1;k<=min(i,j);++k) s[sg[i-k][j-k]]=1;
            for(int k=0;k<3*maxn;++k)
                if(!s[k]) {sg[i][j]=k;break;}
        }
    }

    while(t--){
        int a,b,c,d;scanf("%d%d%d%d",&a,&b,&c,&d);
        int flg=0;
        for(int i=1;i<=a;++i){
            int x,y;scanf("%d%d",&x,&y);
            flg ^= sg[x][y];
        }
        for(int i=1;i<=b;++i){
            int x;scanf("%d",&x);
            flg ^= x%(i+1);
        }
        for(int i=1;i<=c;++i){
            int x;scanf("%d",&x);
            flg ^= max(0,x-i);
        }
        for(int i=1;i<=d;++i){
            scanf("%d",&s4[i]);
        }
        sort(s4+1,s4+d+1);
        s4[0]=-1;
        if(d&1) for(int i=0;i<d;i+=2){
            flg ^= (s4[i+1]-s4[i]-1);
        }
        else for(int i=1;i<d;i+=2)
            flg ^= (s4[i+1]-s4[i]-1);
        if(flg) cout<<"ght"<<endl;
        else cout<<"dxt"<<endl;
    }
}
