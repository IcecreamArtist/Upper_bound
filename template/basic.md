# 排序
## 计数排序
```c++
const int N = 100010;
const int W = 100010;

int n, w, a[N], cnt[W], b[N];

void counting_sort() {
  memset(cnt, 0, sizeof(cnt));
  for (int i = 1; i <= n; ++i) ++cnt[a[i]];
  for (int i = 1; i <= w; ++i) cnt[i] += cnt[i - 1];
  for (int i = n; i >= 1; --i) b[cnt[a[i]]--] = a[i];
}
```
## 基数排序
```c++
const int N = 100010;
const int W = 100010;
const int K = 100;

int n, w[K], k, cnt[W];

struct Element {
  int key[K];
  bool operator<(const Element& y) const {
    // 两个元素的比较流程
    for (int i = 1; i <= k; ++i) {
      if (key[i] == y.key[i]) continue;
      return key[i] < y.key[i];
    }
    return false;
  }
} a[N], b[N];

void counting_sort(int p) {
  memset(cnt, 0, sizeof(cnt));
  for (int i = 1; i <= n; ++i) ++cnt[a[i].key[p]];
  for (int i = 1; i <= w[p]; ++i) cnt[i] += cnt[i - 1];
  // 为保证排序的稳定性，此处循环i应从n到1
  // 即当两元素关键字的值相同时，原先排在后面的元素在排序后仍应排在后面
  for (int i = n; i >= 1; --i) b[cnt[a[i].key[p]]--] = a[i];
  memcpy(a, b, sizeof(a));
}

void radix_sort() {
  for (int i = k; i >= 1; --i) {
    //借助计数排序完成对关键字的排序
    counting_sort(i);
  }
}
```
## 快速排序
```c++
void qsort(int l,int r)
{
    int mid=a[(l+r)/2];
    int i=l,j=r;
    do{
        while(a[i]<mid) i++;
        while(a[j]>mid) j--;
        if(i<=j)
        {
            swap(a[i],a[j]);
            i++;
            j--;
        }
    }while(i<=j);
    if(l<j) qsort(l,j);
    if(i<r) qsort(i,r);
}
```
### 应用：线性找第k大的数
```c++
// 模板的T参数表示元素的类型，此类型需要定义小于（<）运算
template <typename T>
// arr为查找范围数组，rk为需要查找的排名（从0开始），len为数组长度
T find_kth_element(T arr[], int rk, const int len) {
  if (len <= 1) return arr[0];
  // 随机选择基准（pivot）
  const T pivot = arr[rand() % len];
  // i：当前操作的元素
  // j：第一个等于pivot的元素
  // k：第一个大于pivot的元素
  int i = 0, j = 0, k = len;
  // 完成一趟三路快排，将序列分为：小于pivot的元素 ｜ 等于pivot的元素 ｜
  // 大于pivot的元素
  while (i < k) {
    if (arr[i] < pivot)
      swap(arr[i++], arr[j++]);
    else if (pivot < arr[i])
      swap(arr[i], arr[--k]);
    else
      i++;
  }
  // 根据要找的排名与两条分界线的位置，去不同的区间递归查找第k大的数
  // 如果小于pivot的元素个数比k多，则第k大的元素一定是一个小于pivot的元素
  if (rk < j) return find_kth_element(arr, rk, j);
  // 否则，如果小于pivot和等于pivot的元素加起来也没有k多，则第k
  // 大的元素一定是一个大于pivot的元素
  else if (rk >= k)
    return find_kth_element(arr + k, rk - k, len - k);
  // 否则，pivot就是第k大的元素
  return pivot;
}
```
# 二分
```c++
int binary_search(int start, int end, int key) {
  int ret = -1;  // 未搜索到数据返回-1下标
  int mid;
  while (start <= end) {
    mid = start + ((end - start) >> 1);  // 直接平均可能会溢出，所以用这个算法
    if (arr[mid] < key)
      start = mid + 1;
    else if (arr[mid] > key)
      end = mid - 1;
    else {  // 最后检测相等是因为多数搜索情况不是大于就是小于
      ret = mid;
      break;
    }
  }
  return ret;  // 单一出口
}
```
# 倍增
## 应用：LCA
```c++
const int maxn=500005;
int n,m,rt;
int d[maxn],f[maxn][22];
int lg[maxn];
vector<int>G[maxn];

void dfs(int u,int fa){
	d[u]=d[fa]+1;
	f[u][0]=fa;
	//f[u][0]是父亲，然后是第二个祖先，第四个祖先，。。。以此类推 
	for(int i=0;i<lg[d[u]];++i) f[u][i+1]=f[f[u][i]][i];
	for(auto i:G[u]) if(i!=fa) dfs(i,u);
}

int lca(int x,int y){
	if(d[x]<d[y]) swap(x,y);            //x的深度大
	//倍增【1】：先跳到同一深度。只要前者深度大于后者，前者跳到第
	// 2^【log2(深度差)-1】（即：深度差 - 1）个祖先 
	while(d[x]>d[y]) x = f[x][lg[d[x]-d[y]]-1]; 
	if(x == y) return x;               //此时深度相同。两者相
	// 等的话说明原先的y是x的祖先。此时的x或y就是lca
	/*
	倍增【2】：找lca。
	【二进制思维】从大到小尝试往上跳 如果两者祖先不相等就往上跳：
	x跳到第2^【log2(x的深度)-1】（即：x的深度-1）个祖先 
	鉴于x和y的深度一直一样，因此如果祖先相等，说明可能超过了lca的位置。
	*/
	for(int k = lg[d[x]]-1;k>=0;--k)           
		if(f[x][k]!=f[y][k]) x = f[x][k],y = f[y][k]; 
	return f[x][0];                   
}

int main(){
	scanf("%d%d%d",&n,&m,&rt);
	for(int i=1;i<=n;++i) lg[i] = lg[i-1] + (1<<lg[i-1] == i);
	int x,y,k;
	for(int i=1;i<n;++i) scanf("%d%d",&x,&y),
	G[x].push_back(y),G[y].push_back(x);
	// 常数优化，求出 log2(x) + 1 
	dfs(rt,0);
	for(int i=1;i<=m;++i) scanf("%d%d",&x,&y),printf("%d\n",lca(x,y));
}
```
# 分治
题意：有一个长度为n的数组a, (a1,a2,…,an)，选一个整数x，让 (ai^x)中的最大值最小.
思路：看了题解，大意就是，因为最高位更小的话得到的数字一定更小，所以我们从最高位开始考虑。但是最高位的选择又会影响到下一位的选择，因此我们使用dfs的逻辑dp遍历。
只考虑当前的位。如果数组中所有数字在这一位都为0或都为1，那么答案的这一位（最小的最大值）可以为0。但如果这一位有的是1，有的是0，那么不管我们在这一位选择1还是0，答案的这一位必然都是1。但是这一位到底选择1还是0呢？这个决定会影响接下来跟随哪些数字寻找最大值。那么我们就用dfs的dp，找寻两个情况的最小值即可。
```c++
#include<bits/stdc++.h>
using namespace std;
int a[100004];
int dfs(vector<int> v,int bit){
	if(bit<0) return 0;
	vector<int> v0,v1;
	for(auto &i:v){
		if((i>>bit)&1) v1.push_back(i);
		else v0.push_back(i);
	}
	if(v0.empty()) return dfs(v1,bit-1);
	if(v1.empty()) return dfs(v0,bit-1);
	return min(dfs(v0,bit-1),dfs(v1,bit-1))+(1<<bit);
}
int main(){
	int n;
	vector<int> v;
	scanf("%d",&n);
	for(int i=1;i<=n;i++){
		scanf("%d",&a[i]);
		v.push_back(a[i]);
	}
	printf("%d\n",dfs(v,30));
}

```

# 分块
# 压缩图
注意：有时候可能要专门加进去一个边界点。
```c++
int main(){
	int n,m;
	scanf("%d%d",&n,&m);
	for(int i=0;i<n;i++){
		scanf("%d",&a[i]);
		scanf("%d",&b[i]);
		scanf("%d",&c[i]);
		xx.push_back(a[i]);
		xx.push_back(b[i]);
		yy.push_back(c[i]);
	}
	for(int i=0;i<m;i++){
		scanf("%d",&D[i]);
		scanf("%d",&e[i]);
		scanf("%d",&f[i]);
		xx.push_back(D[i]);
		yy.push_back(e[i]);
		yy.push_back(f[i]);
	}
	sort(xx.begin(),xx.end());//先按照从小到大sort
	sort(yy.begin(),yy.end());
	xx.resize(unique(xx.begin(),xx.end())-xx.begin());//然后把重复项去掉！（压缩）
	yy.resize(unique(yy.begin(),yy.end())-yy.begin());
	for(int i=0;i<n;i++){
		//接下来对于每条线段，先把线段查找出来，再用新的压缩坐标标记。
		//新的压缩坐标是vector中的下标。
		int aa=lower_bound(xx.begin(),xx.end(),a[i])-xx.begin();
		int bb=lower_bound(xx.begin(),xx.end(),b[i])-xx.begin();
		int cc=lower_bound(yy.begin(),yy.end(),c[i])-yy.begin();
		for(int x=aa;x<bb;x++){
			u[x][cc-1]=1;
			d[x][cc]=1;
		}
	}
	for(int i=0;i<m;i++){
		int dd=lower_bound(xx.begin(),xx.end(),D[i])-xx.begin();
		int ee=lower_bound(yy.begin(),yy.end(),e[i])-yy.begin();
		int ff=lower_bound(yy.begin(),yy.end(),f[i])-yy.begin();
		for(int y=ee;y<ff;y++){
			l[dd][y]=1;
			r[dd-1][y]=1;
		}
	}
	//找到牛位置
	int x=lower_bound(xx.begin(),xx.end(),0)-xx.begin();
	int y=lower_bound(yy.begin(),yy.end(),0)-yy.begin();
	//用左下角方格做起点
	bfs(x-1,y-1);
	if(vis[0][0]){
		printf("INF\n");return 0;
	}
	ll ans=0;
	for(int i=0;i<xx.size()-1;i++){
		for(int j=0;j<yy.size()-1;j++){
			if(vis[i][j]) ans+=(ll)(xx[i+1]-xx[i])*(yy[j+1]-yy[j]); 
			// 用后一个方块的最前端减去当前方块最前端代表当前方块的边长
		}
	}
	printf("%lld\n",ans);
}

```
# 二叉树
已知先序中序求后序
```c++
#include<bits/stdc++.h>
using namespace std;
void dfs(string xx,string zx){
	if(!xx.size()) return;
	int pos=zx.find(xx[0]);
	dfs(xx.substr(1,pos),zx.substr(0,pos));
	dfs(xx.substr(pos+1),zx.substr(pos+1));
	printf("%c",xx[0]);
}
int main(){
	string xx,zx;
	cin>>zx>>xx;
	dfs(xx,zx);
	printf("\n");
}

```
已知后序中序求前序
```c++
#include<bits/stdc++.h>
using namespace std;
void dfs(string zx,string hx){
	if(!zx.size())return;
	int pos=zx.find(hx[hx.size()-1]);
	printf("%c",zx[pos]);
	dfs(zx.substr(0,pos),hx.substr(0,pos));
	dfs(zx.substr(pos+1),hx.substr(pos,hx.size()-pos-1));
} 
int main(){
	string zx,hx;
	cin>>zx>>hx;
	dfs(zx,hx);
	printf("\n");
}

```
已知先序后序求中序数目
```c++
#include<bits/stdc++.h>
using namespace std;
int main(){
	char a[400],b[400];
	scanf("%s%s",a,b);
	long long n=0;
	for(int i=0;i<strlen(a)-1;i++)
		for(int j=1;j<strlen(b);j++)
			if(b[j]==a[i]&&a[i+1]==b[j-1]){
				n++;continue;
			}
	printf("%lld\n",1<<n);
}

```
# BST
```c++
#include<bits/stdc++.h>
using namespace std;
const int inf=0x7fffffff;
//BST
int cont=0,ans;
struct node{
	int ls,rs,val,siz,cnt;
}tree[500004];
//x 现在的节点 y 值 
void insert(int x,int y){
	tree[x].siz++;
	if(tree[x].val==y){
		tree[x].cnt++;
		return;
	}
	if(tree[x].val>y){
		if(tree[x].ls!=0)
			insert(tree[x].ls,y);
		else{
			cont++;
			tree[x].ls=cont;
			tree[cont].cnt=tree[cont].siz=1;
			tree[cont].val=y;
		}
	}
	else{
		if(tree[x].rs!=0)
			insert(tree[x].rs,y);
		else{
			cont++;
			tree[x].rs=cont;
			tree[cont].cnt=tree[cont].siz=1;
			tree[cont].val=y;
		}
	}
}
//x 现在的节点 val 值 
int search(int x,int val){
	if(x==0) return 0;
	if(val==tree[x].val)return tree[tree[x].ls].siz;
	if(val>tree[x].val)return tree[tree[x].ls].siz+tree[x].cnt+search(tree[x].rs,val);
	if(val<tree[x].val)return search(tree[x].ls,val);
}
//y排名 
int search2(int x,int y){
	if(x==0) return inf;
	if(tree[tree[x].ls].siz>=y)return search2(tree[x].ls,y);
	if(tree[tree[x].ls].siz+tree[x].cnt>=y)return tree[x].val;
	if(y>tree[tree[x].ls].siz)return search2(tree[x].rs,y-tree[tree[x].ls].siz-tree[x].cnt);
}
//找前驱
void searchf(int x,int val){
	if(x==0) return;
	if(tree[x].val<val){
		ans=max(ans,tree[x].val);
		searchf(tree[x].rs,val);
	}
	if(tree[x].val>=val){
		searchf(tree[x].ls,val);
	}
}
//找后继
void searchl(int x,int val){
	if(x==0) return;
	if(tree[x].val>val){
		ans=min(ans,tree[x].val);
		searchl(tree[x].ls,val);
	}
	if(tree[x].val<=val){
		searchl(tree[x].rs,val);
	}
}
int main(){
	int n;
	int type,num;
	scanf("%d",&n);
	while(n--){
		scanf("%d %d",&type,&num);
		if(type==5){
			if(cont==0){
				cont++;
				tree[1].cnt=tree[1].siz=1;
				tree[1].val=num;
			}
			else insert(1,num);
		}
		if(type==1)printf("%d\n",search(1,num)+1);//查找排名 
		if(type==2)printf("%d\n",search2(1,num));//依据排名查找数字 
		if(type==3)ans=-inf,searchf(1,num),printf("%d\n",ans);//查找前驱 
		if(type==4)ans=inf,searchl(1,num),printf("%d\n",ans);//查找后继 
	}
}

```
# 一元三次方程组
```c++
#include<bits/stdc++.h>
using namespace std;
double a,b,c,d;

double f(double x){
	return a*x*x*x+b*x*x+c*x+d;
}

int main(){
	int cnt=0;
	double ans;
	double l,r,m,tmp1,tmp2;
	scanf("%lf%lf%lf%lf",&a,&b,&c,&d);
	for(int i=-100;i<100;i++){
		l=i,r=i+1;
		tmp1=f(l),tmp2=f(r);
		if(!tmp1){
			printf("%.2lf ",l);
			cnt++;
			continue;
		}
		if(tmp1*tmp2<0){
			while(r>=l){
            	m=(l+r)/2;
				if(f(m)*f(l)<=0)r=m-0.0001,ans=m;
				else l=m+0.0001,ans=m;
			}
			printf("%.2lf ",ans);
			cnt++;
		}
		if(cnt==3)break;
	}
	l=100;
	if(cnt==2)printf("%.2lf ",l);
	return 0;
}

```
