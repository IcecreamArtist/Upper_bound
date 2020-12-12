# 1. 莫队算法
## 特点
- O2优化
- 离线询问，进行排序与分块，对每个块集中处理
- 已知一个区间的答案，推出所有区间的答案
- 时间复杂度$O(n\sqrt{n})$
## 适用问题
- 给一个区间，多次查询（O(n^2)），允许离线
- 答案从(i,j)转换到(i-1,j)、(i,j-1)、(i+1,j)、(i,j+1)的复杂度是O(1)
- 即：已知一个区间的答案，求下一个区间的答案时，可直接根据所在数据$a[i]$处理更新
## 普通莫队
- 将i分块，在每个块内，从1到n扫j。

```cpp
/*
 * 莫队模板题
 * 一个数组,N(5e4)个数字,M(5e4)次查询：区间中抽两个数字,两者相等的概率。
 * 先上暴力做法
 */
#include <bits/stdc++.h>
using namespace std;
int n,m,maxn;
typedef long long ll;
const int N = 5e4+5;
int a[N],cnt[N];
ll sum;
struct node{
    int l,r,id;
    bool operator < (const node &b) const{
        if(l/maxn!=b.l/maxn) return l<b.l;
        return r<b.r;
    }
}q[N];

void add(int i){
    sum += cnt[i];
    cnt[i]++;
}

void del(int i){
    cnt[i]--;
    sum -= cnt[i];
}

ll gcd(ll a,ll b){
    return b?gcd(b,a%b):a;
}
ll ans1[N],ans2[N];

int main(){
    scanf("%d%d",&n,&m);
    maxn = sqrt(n);
    for(int i=1;i<=n;++i) scanf("%d",&a[i]);
    for(int i=0;i<m;++i) scanf("%d%d",&q[i].l,&q[i].r),q[i].id = i;
//    下面是暴力做法
//    for(int i=1;i<=m;++i) {
//        memset(cnt,0,sizeof(cnt));
//        int fz = 0,fm = (q[i].r-q[i].l+1)*(q[i].r-q[i].l)/2;
//        if(fm==0) {puts("0/1");continue;}
//        for(int j = q[i].l;j<=q[i].r;++j){
//            fz += cnt[a[j]];
//            cnt[a[j]]++;
//        }
//        if(fz==0) {puts("0/1");continue;}
//        int g = gcd(fz,fm);
//        printf("%d/%d\n",fz/g,fm/g);
//    }
    // 下面是莫队算法
    sort(q,q+m);
    for(int i=0,l=1,r=0;i<m;++i){
        if(q[i].l==q[i].r){ans1[q[i].id]=0,ans2[q[i].id]=1;continue;}
        while(l>q[i].l) add(a[--l]);
        while(r<q[i].r) add(a[++r]);
        while(l<q[i].l) del(a[l++]);
        while(r>q[i].r) del(a[r--]);
        ans1[q[i].id] = sum;
        ans2[q[i].id] = (ll)(r-l+1)*(r-l)/2;
    }
    for(int i=0;i<m;++i){
        if(ans1[i]!=0){
            ll g = gcd(ans1[i],ans2[i]);
            ans1[i] /= g, ans2[i] /= g;
        }
        else ans2[i] = 1;
        printf("%lld/%lld\n",ans1[i],ans2[i]);
    }
}
```
## 带修改莫队
[JNUOJ1228 486的区间](https://jnuacmoj.jnu.edu.cn/JudgeOnline/problem.php?id=1228)
带修改莫队的特殊之处，就在于多了一个特征：时间。
将i和j分块，对每个i块，再对每个j块，从1~n扫t。
```cpp
/*
 * 题意：n个数，每个数为ai，m个询问，求：l到r区间ai*num(ai)的和
 * 带修莫队模板
 */
#include <bits/stdc++.h>
using namespace std;
const int N = 1e4+6;
int n,m,a[N],sqn,qcnt,ccnt,num[N*10];
long long ans[N],sm;
 
struct query{
    int l,r,id,t;
    bool operator <(const query &b)const{
        if(l/sqn!=b.l/sqn) return l/sqn<b.l/sqn;
        if(r/sqn!=b.r/sqn) return r/sqn<b.r/sqn;
        return t<b.t;
    }
}q[N];
 
struct modify{
    int pos,pre,nex;
}c[N];
 
template <typename _Tp>
inline void IN(_Tp& dig){
    char c;dig = 0;
    while(c=getchar(),!isdigit(c));
    while(isdigit(c)) dig = dig * 10 + c -'0',c = getchar();
}
 
void add(int p){
    sm-=(long long)num[a[p]]*num[a[p]]*a[p];
    sm+=(long long)(++num[a[p]])*num[a[p]]*a[p];
}
 
void del(int p){
    sm-=(long long)num[a[p]]*num[a[p]]*a[p];
    sm+=(long long)(--num[a[p]])*num[a[p]]*a[p];
}
 
int main(){
    //freopen("input.txt","r",stdin);
    IN(n);IN(m);sqn = ceil(pow(n,(double)2/(double)3));
    for(int i=1;i<=n;++i) IN(a[i]);
    for(int i=1,x,y;i<=m;++i) {
        char str[4];
        if(scanf("%s",str),IN(x),IN(y),str[0]=='Q')
            q[++qcnt].l = x,q[qcnt].r = y,q[qcnt].id = qcnt,q[qcnt].t = ccnt;
        else
            c[++ccnt].pos = x,c[ccnt].nex = y,c[ccnt].pre = a[x],a[x] = y;
    }
    sort(q+1,q+1+qcnt);
    for(int i=ccnt;i>=1;--i) a[c[i].pos] = c[i].pre;
    for(int i=1,l=1,r=0,t=0;i<=qcnt;++i){
        while(l>q[i].l) add(--l);
        while(r<q[i].r) add(++r);
        while(l<q[i].l) del(l++);
        while(r>q[i].r) del(r--);
        while(t<q[i].t){
            int p = c[++t].pos;
            if(l<=p&&r>=p) del(p);
            a[p] = c[t].nex;
            if(l<=p&&r>=p) add(p);
        }
        while(t>q[i].t){
            int p = c[t].pos;
            if(l<=p&&r>=p) del(p);
            a[p] = c[t].pre;
            if(l<=p&&r>=p) add(p);
            t--;
        }
        ans[q[i].id] = sm;
    }
    for(int i=1;i<=qcnt;++i) printf("%lld\n",ans[i]);
}
 
/**************************************************************
    Problem: 1228
    User: IcecreamArtist
    Language: C++
    Result: Accepted
    Time:36 ms
    Memory:2696 kb
****************************************************************/
```

# 2. CDQ分治
## 三维偏序
二维偏序可先以其中一个特征排序，再用树状数组对第二个特性进行统计，累计答案。
二维偏序也可以用分治来做，就不需要套BIT了，见用归并排序求逆序对。（逆序对也可以用树状数组来做，BIT和分治本质上都是降维）
三维偏序则在分治的基础上加多了一维，则是分治上用树状数组处理第二个特性。
用第一个特征将数组排序，然后对这个区间进行分治。每到一个区间（l,r），都对（l,mid）和（mid+1,r）内部再以第二个特征进行排序。然后利用左区间，进行右区间的贡献累计（也可能是利用右区间，对左区间进行贡献累计，取决于第一维特征）。这个累计过程使用BIT，道理类似于归并排序求逆序对。记得BIT用完过后要进行撤销。
## 
[洛谷3157 动态逆序对](https://www.luogu.com.cn/problem/P3157)

```cpp
/*
 * CDQ分治模板题：三维偏序
 * 题意：给一个数组(1e5)，m(5e4)次操作，每次删除指定数字，每次删除之前输出逆序对个数。
 * 思路：每次删除指定数字，等价于时间反向，每次加入指定数字。
 * 对于一个数字的加入，该时间点增加的逆序对数等于比他早加入、位置在他之前、比他大的数+比他早加入、位置在他之后、比他小的数。
 * 分治过程中用左区间对右区间进行贡献即可。
 */
#include <bits/stdc++.h>
using namespace std;
const int maxn = 1e5+6;
typedef long long ll;
int n,m,c[maxn],p[maxn];
struct data{int t,p,v;}a[maxn],b[maxn];
ll ans[maxn];
inline int read(){int x=0,f=1;char ch=getchar();while(ch<'0'||ch>'9') {if(ch=='-') f=-1;ch=getchar();}while(ch>='0'&&ch<='9') x=x*10+ch-'0',ch=getchar();return x*f;}
inline void insert(int x,int y){for(;x<=n;x+=(x&(-x))) c[x]+=y;}
inline int query(int x){int sm = 0;for(;x;x-=(x&(-x))) sm+=c[x];return sm;}

void cdq(int l,int r){
    if(l>=r) return;
    int mid = (l+r)>>1,tmp,ll=l-1,rr=mid;
    for(int i=l;i<=r;++i){
        if(a[i].t<=mid) b[++ll] = a[i];
        else b[++rr] = a[i];
    }
    // 将数组重新排序：mid左时间都早于右边；同时，左之中的数字按照位置从前往后排序，右之中的数字也是按位置从前往后排序
    for(int i=l;i<=r;++i) a[i]=b[i];
    tmp=l;
    for(int i=mid+1;i<=r;++i){
        // 每一个右边的，都要计算上左边中出现在他之前且比他大的
        for(;tmp<=mid&&a[tmp].p<a[i].p;tmp++) insert(a[tmp].v,1);
        ans[a[i].t]+=tmp-l-query(a[i].v);
    }
    for(int i=l;i<tmp;++i) insert(a[i].v,-1); // 撤销BIT
    tmp=mid;
    for(int i=r;i>=mid+1;--i){
        // 每一个右边的，还要计算上左边中出现在他之后且比他小的
        for(;tmp>=l&&a[tmp].p>a[i].p;tmp--) insert(a[tmp].v,1);
        ans[a[i].t]+=query(a[i].v-1);
    }
    for(int i=tmp+1;i<=mid;++i) insert(a[i].v,-1); // 撤销BIT
    cdq(l,mid);cdq(mid+1,r);
}

int main(){
    n=read(),m=read();
    int time = n,x;
    for(int i=1;i<=n;++i) a[i].v=read(),a[i].p=i,p[a[i].v]=i;
    for(int i=1;i<=m;++i) x=read(),a[p[x]].t=time--;
    for(int i=1;i<=n;++i) if(!a[i].t) a[i].t=time--;
    cdq(1,n);
    // 依次删除每个数，等价于反向加入每个数
    for(int i=1;i<=n;++i) ans[i]+=ans[i-1]; // 反向加入每个数，到第几次的结果。
    for(int i=n;i>=n-m+1;--i) printf("%lld\n",ans[i]);
}

```
# 快读快输
```c++
inline int read() {
  int x = 0, f = 1;
  char ch = getchar();
  while (ch < '0' || ch > '9') {
    if (ch == '-') f = -1;
    ch = getchar();
  }
  while (ch <= '9' && ch >= '0') {
    x = 10 * x + ch - '0';
    ch = getchar();
  }
  return x * f;
}
void print(int x) {
  if (x < 0) putchar('-'), x = -x;
  if (x >= 10) print(x / 10);
  putchar(x % 10 + '0');
}

```
# 尺取法
```c++
int main(){
	ull s;
	int t;
	scanf("%d",&t);
	while(t--){
		int n;
		scanf("%d%llu",&n,&s);
		for(int i=1;i<=n;i++){
			scanf("%d",&a[i]);
		}
		deque<int>q;
		q.clear();
		ull sum=0;
		int i=1;
		int len=0,ans=inf;
		while(1){
			while(i<=n&&sum<s){
				sum+=a[i];
				q.push_back(a[i++]);
				len++;
			}
			if(sum<s) break;
			ans=min(ans,len);
			len--;
			sum-=q.front();
			q.pop_front();
		}
		if(ans==inf) printf("0\n");
		else printf("%d\n",ans);
	}
}
```
```c++
int main(){
	int p;
	int num=0;
	scanf("%d",&p);
	map<int,int>mp; 
	for(int i=1;i<=p;i++){
		scanf("%d",&a[i]);
		if(!mp[a[i]]) num++;
		mp[a[i]]=1;
	}
	mp.clear();
	deque<int>q;
	int cnt=0;
	int i=1;
	int ans=inf;
	while(1){
		while(i<=p&&cnt<num){
			if(mp[a[i]]==0) cnt++;
			mp[a[i]]++;
			q.push_back(a[i++]);
		}
		if(cnt<num) break;
		int len=q.size();
		ans=min(ans,len);
		mp[q.front()]--;
		if(mp[q.front()]==0) cnt--;
		q.pop_front();
	}
	printf("%d\n",ans);
}

```
