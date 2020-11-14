#include <bits/stdc++.h>
using namespace std;
/*
 * 题意：一个地图，S是起点，E是终点。
 * 从S走到E，
 */
const int N = 2e3+5;
int n,m;
int dir[4][2]={
        {1,0},
        {0,-1},
        {0,1},
        {-1,0}
};
bool check(int x,int y){
    return x<=n&&x>=1&&y<=m&&y>=1;
}
char mp[N][N];
int vis[N][N];  // 记录曾经是否走过这个点
int sx,sy;  // start
struct node{
    int x,y;
    int fa;   // 前驱
};
node path[N*N];  // 要记录路径，须手工模拟队列
char ans[N*N];

void print(int pos){
    node cur,last;
    int cnt = 0;
    cur = path[pos];
    pos = cur.fa;
    while(pos!=-1){
        last = path[pos];
        if(cur.x>last.x) ans[cnt++] = 'D';
        if(cur.x<last.x) ans[cnt++] = 'U';
        if(cur.y>last.y) ans[cnt++] = 'R';
        if(cur.y<last.y) ans[cnt++] = 'L';
        pos = last.fa;
        cur = last;
    }
    printf("%d\n",cnt);
    for(int i=cnt-1;i>=0;--i) printf("%c",ans[i]);
}

void BFS(){
    node start,next;
    start.x = sx,start.y = sy,start.fa=-1;
    int front = 0,rear = 0;  // 队列头以及队列尾
    path[rear++] = start;
    vis[sx][sy] = 1;
    while(front<rear){
        start = path[front++];
        for(int i=0;i<4;++i){
            next.x=start.x+dir[i][0],next.y=start.y+dir[i][1];
            if(check(next.x,next.y)&&!vis[next.x][next.y]&&mp[next.x][next.y]!='X'){
                next.fa = front-1;
                path[rear++]=next;
                vis[next.x][next.y]=1;
                if(mp[next.x][next.y]=='E') {
                    print(rear-1);return;
                }
            }
            while(check(next.x,next.y)&&mp[next.x][next.y]=='X')
                next.x=next.x+dir[i][0],next.y=next.y+dir[i][1];
            if(check(next.x,next.y)&&!vis[next.x][next.y]){
                next.fa = front-1;
                path[rear++]=next;
                vis[next.x][next.y]=1;
                if(mp[next.x][next.y]=='E'){
                    print(rear-1);return;
                }
            }
        }
    }
    puts("-1");
}

int main(){
    scanf("%d%d",&n,&m);
    for(int i=1;i<=n;++i) scanf("%s",mp[i]+1);
    for(int i=1;i<=n;++i) for(int j=1;j<=m;++j) if(mp[i][j]=='S') sx = i,sy = j;
    BFS();
}