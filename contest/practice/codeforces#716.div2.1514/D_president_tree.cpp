#include<bits/stdc++.h>
using namespace std;
const int maxn = 2e6+5;
int A[maxn],root[maxn],cnt=1;
int c[maxn],L[maxn],ori[maxn];

#define num(x) tree[x].num
#define ls(x) tree[x].ls
#define rs(x) tree[x].rs

struct{
    int num;
    int ls,rs;
}tree[maxn<<2];

void build(){
    root[0] = 1;
    num(1) = 0;
    ls(1)=1;
    rs(1)=1;
}

void update(int a,int l,int r,int rt,int nrt){
    num(nrt) = num(rt);
    if(l==r) num(nrt)++;
    else{
        ls(nrt) = ls(rt),rs(nrt) = rs(rt);
        int mid = (l+r)>>1;
        if(a<=mid) ls(nrt)=++cnt,update(a,l,mid,ls(rt),ls(nrt));
        else rs(nrt)=++cnt,update(a,mid+1,r,rs(rt),rs(nrt));
        num(nrt) = num(ls(nrt)) + num(rs(nrt));
    }
}

// 查询区间第a大
int query(int a,int l,int r,int p,int q){
    //  cout<<a<<" "<<l<<" "<<r<<" "<<p<<" "<<q<<endl;
    if(l==r) return num(q)-num(p); // 该数的出现次数
    else{
        int mid = (l+r)>>1;
        if(a<=num(ls(q))-num(ls(p))) return query(a,l,mid,ls(p),ls(q));
        else return query(a-(num(ls(q))-num(ls(p))),mid+1,r,rs(p),rs(q));
    }
}

int discretize(int n){
    for(int i=0;i<n;++i) c[i]=A[i];
    sort(c,c+n);
    int len = unique(c,c+n)-c;
    int maxx = 0;
    for(int i=0;i<n;++i){
        L[i] = lower_bound(c,c+len,A[i])-c+1; // 查找
        maxx = max(maxx,L[i]);
        ori[L[i]] = A[i]; // 离散化后的数对应的原数
    }
    return maxx;
}

int main(){
    int n,m;scanf("%d%d",&n,&m);
    for(int i=0;i<n;++i) scanf("%d",&A[i]);
    build();
    int maxx = discretize(n);

    for(int i=0;i<n;++i){
        root[i+1] = ++cnt;
        update(L[i],1,maxx,root[i],root[i+1]);
    }

    for(int i=1;i<=m;++i){
        int l,r;scanf("%d%d",&l,&r);
        int ans = query(((l+r)>>1)-l+1,1,maxx,root[l-1],root[r]); // 查询中位数
      //  cout<<ans<<endl;
        cout<<max(1,ans-(r-l+1-ans))<<endl;
    }
}
