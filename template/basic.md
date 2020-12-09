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
## CDQ分治

# 分块