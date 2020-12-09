# 数论
# 素数
## Miller-Rabin素性测试
时间复杂度O(klog^3(n))
```c++
bool millerRabbin(int n) {
  if (n < 3) return n == 2;
  int a = n - 1, b = 0;
  while (a % 2 == 0) a /= 2, ++b;
  // test_time 为测试次数,建议设为不小于 8
  // 的整数以保证正确率,但也不宜过大,否则会影响效率
  for (int i = 1, j; i <= test_time; ++i) {
    int x = rand() % (n - 2) + 2, v = quickPow(x, a, n);
    if (v == 1 || v == n - 1) continue;
    for (j = 0; j < b; ++j) {
      v = (long long)v * v % n;
      if (v == n - 1) break;
    }
    if (j >= b) return 0;
  }
  return 1;
}
```
## 求因子数一定的最小数
```c++
#include <stdio.h>
#define ULL unsigned long long
#define INF ~0ULL
ULL p[16] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53};

ULL ans;
ULL n;

// depth: 当前在枚举第几个素数。num: 当前因子数。
// temp: 当前因子数量为 num
// 的时候的数值。up：上一个素数的幂，这次应该小于等于这个幂次嘛
void dfs(ULL depth, ULL temp, ULL num, ULL up) {
  if (num > n || depth >= 16) return;
  if (num == n && ans > temp) {
    ans = temp;
    return;
  }
  for (int i = 1; i <= up; i++) {
    if (temp / p[depth] > ans) break;
    dfs(depth + 1, temp = temp * p[depth], num * (i + 1), i);
  }
}

int main() {
  while (scanf("%llu", &n) != EOF) {
    ans = INF;
    dfs(0, 1, 1, 64);
    printf("%llu\n", ans);
  }
  return 0;
}
```
## 求 n 以内因子数最多的数
```c++
#include <cstdio>
#include <iostream>
#define ULL unsigned long long

int p[16] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53};
ULL n;
ULL ans, ans_num;  // ans 为 n 以内的最大反素数（会持续更新），ans_sum 为 ans
                   // 的因子数。

void dfs(int depth, ULL temp, ULL num, int up) {
  if (depth >= 16 || temp > n) return;
  if (num > ans_num) {
    ans = temp;
    ans_num = num;
  }
  if (num == ans_num && ans > temp) ans = temp;
  for (int i = 1; i <= up; i++) {
    if (temp * p[depth] > n) break;
    dfs(depth + 1, temp *= p[depth], num * (i + 1), i);
  }
  return;
}

int main() {
  while (scanf("%llu", &n) != EOF) {
    ans_num = 0;
    dfs(0, 1, 1, 60);
    printf("%llu\n", ans);
  }
  return 0;
}
```
# exgcd
```c++
inline ll exgcd(ll a,ll b,ll &x,ll &y){
    if(!b) {x=1,y=0;return a;}
    ll g = exgcd(b,a%b,y,x);
    y-=x*(a/b);
    return g;
}
```
非递归版本
```c++
int gcd(int a, int b, int& x, int& y) {
    x = 1, y = 0;
    int x1 = 0, y1 = 1, a1 = a, b1 = b;
    while (b1) {
        int q = a1 / b1;
        tie(x, x1) = make_tuple(x1, x - q * x1);
        tie(y, y1) = make_tuple(y1, y - q * y1);
        tie(a1, b1) = make_tuple(b1, a1 - q * b1);
    }
    return a1;
}
```
# 欧拉函数
## 性质
1. $\varphi (n)=n \prod_{i=1}^{s}(1-\frac{1}{p_i})$
2. $n=\sum_{d|n}\varphi(d)$
3. 积性函数：if $gcd(a,b)=1$ then $\varphi(a\times b)=\varphi(a)\times \varphi(b)$
## 筛法
1. 求单个欧拉函数的值
   （可以用Pollard Rho算法优化）
   time complexity: $O(n)$
```cpp
int euler_phi(int n) {
  int m = int(sqrt(n + 0.5));
  int ans = n;
  for (int i = 2; i <= m; i++)
    if (n % i == 0) {
      ans = ans / i * (i - 1);
      while (n % i == 0) n /= i;
    }
  if (n > 1) ans = ans / n * (n - 1);
  return ans;
}
```
2. 筛一堆

(a) time complexity: $O(nloglogn)$，基于埃氏筛，不需要记录质数
```cpp
// 不用筛质数
void phi_table(int n, int* phi) {
  for (int i = 2; i <= n; i++) phi[i] = 0;
  phi[1] = 1;
  for (int i = 2; i <= n; i++)
    if (!phi[i])
      for (int j = i; j <= n; j += i) {
        if (!phi[j]) phi[j] = j;
        phi[j] = phi[j] / i * (i - 1);
      }
}
```
(b) time complexity: $O(n)$，基于欧拉筛（线性筛），需要记录质数

```cpp
void init() {
  phi[1] = 1;
  for (int i = 2; i < MAXN; ++i) {
    if (!vis[i]) {
      phi[i] = i - 1;
      pri[cnt++] = i;
    }
    for (int j = 0; j < cnt; ++j) {
      if (1ll * i * pri[j] >= MAXN) break;
      vis[i * pri[j]] = 1;
      if (i % pri[j]) {
        phi[i * pri[j]] = phi[i] * (pri[j] - 1);
      } else {
        phi[i * pri[j]] = phi[i] * pri[j];
        break;
      }
    }
  }
}
```
## 欧拉定理&费马小定理
- 欧拉定理：
  若$\gcd(a,m)=1$，则$a^{\varphi (m)} \equiv 1 \pmod m$.

- 费马小定理：
  欧拉定理的$m$为**质数**$p$，则$\varphi(p)=p-1$，
  故有$a^{p-1} \equiv 1 \pmod p$.

- 扩展欧拉定理：
  $a$为任意整数，$b$和$m$是正整数。
  $$a^b\equiv \begin{cases} a^{b\mod \varphi(m)}\, &\ \gcd(a,m)=1\\a^{b\mod \varphi(m)+\varphi(m)}&\ \gcd(a,m)\ne 1\end{cases} \pmod m$$

## 求欧拉函数的反函数

```cpp

#include<bits/stdc++.h>

using namespace std;
typedef long long ll;
const int N = 1e5+5;
ll ans=200000000000;
int pri[N],vis[N],cnt,m;

void get_prime(){
    for(int i=2;i<=50000;++i){
        if(!vis[i]) pri[++cnt] = i;
        for(int j=1;j<=cnt&&i*pri[j]<=50000;++j){
            vis[i*pri[j]]=1;
            if(i%pri[j]==0) break;
        }
    }
    // 因为空间问题，只能筛出sqrt n的质数. 然后我们在dfs过程中对huge prime进行特判。
}

bool isprime(int n){
    int k=(int)sqrt(n+0.5);
    for(int i=2;i<=k;++i) if(n%i==0) return 0;
    return 1;
}

void dfs(int k,int rem,ll pro){
    // now at k th prime, the previous product is pro, the remainder is rem
    if(pro>=ans) return;  // if the current product is already larger, cut the edge
    if(rem==1){ans = pro;return;} // if the product is not larger, and it can divide n, update
    if(rem>m&&isprime(rem+1)) ans = min(ans,pro*(rem+1)); // huge prime
    for(int i=k+1;pri[i]-1<=m;++i){
        if(pri[i]-1>rem) break;
        if(rem%(pri[i]-1)==0){
            int x = rem/(pri[i]-1);ll y = pro*pri[i];
            dfs(i,x,y);
            while(x%pri[i]==0){
                x/=pri[i];y*=pri[i];   // multiply the prime
                dfs(i,x,y);
            }
        }
    }
}

int main(){
    int k;scanf("%d",&k);
    m = (int)sqrt(k+0.5);
    get_prime();
    dfs(0,k,1);
    if(ans==200000000000) cout<<0<<endl;
    else cout<<ans<<endl;
}

```

# 逆元
- 单个求法：扩欧、费马小定理（要求p是质数），时间复杂度都是$O(log n)$

```cpp
// 扩欧
void exgcd(int a, int b, int& x, int& y) {
  if (b == 0) {
    x = 1, y = 0;
    return;
  }
  exgcd(b, a % b, y, x);
  y -= a / b * x;
}
// 快速幂（费马小定理）
inline int qpow(long long a, int b) {
  int ans = 1;
  a = (a % p + p) % p;
  for (; b; b >>= 1) {
    if (b & 1) ans = (a * ans) % p;
    a = (a * a) % p;
  }
  return ans;
}
```

- 批量求法：1~n线性求$O(n)$、n个独立数求$O(n+log p)$

```cpp
// 线性
inv[1] = 1;
for (int i = 2; i <= n; ++i) {
  inv[i] = (long long)(p - p / i) * inv[p % i] % p;
}

// 独立数
s[0] = 1;
for (int i = 1; i <= n; ++i) s[i] = s[i - 1] * a[i] % p;
sv[n] = qpow(s[n], p - 2);
// 当然这里也可以用 exgcd 来求逆元,视个人喜好而定.
for (int i = n; i >= 1; --i) sv[i - 1] = sv[i] * a[i] % p;
for (int i = 1; i <= n; ++i) inv[i] = sv[i] * s[i - 1] % p;
```
# 数论分块
### 理论：
数论分块用于快速处理形如$\sum_{i}^{n}\lfloor\frac{n}{i}\rfloor$的式子。复杂度为$O(\sqrt{n})$（大概）。
通过观察可以发现，对于同一个$n$，$\lfloor\frac{n}{i}\rfloor$的值随着$i$呈现**块状分布**。对于任意一个$i( i\leq n)$，能找到一个最大的$j(i\leq j\leq n)$，使得$\lfloor\frac{n}{i}\rfloor=\lfloor\frac{n}{j}\rfloor$，在这之间的$k(i\leq k\leq j)$都有$\lfloor\frac{n}{k}\rfloor=\lfloor\frac{n}{i}\rfloor$。
而这个$j=\lfloor\frac{n}{\lfloor\frac{n}{i}\rfloor}\rfloor$。
证明：
$\lfloor\frac{n}{i}\rfloor\leq \frac{n}{i}$
$\implies\lfloor\frac{n}{\lfloor\frac{n}{i}\rfloor}\rfloor\geq\lfloor\frac{n}{\frac{n}{i}}\rfloor=\lfloor i\rfloor=i$
$\implies j\geq i$
因此，每次以$[i,j]$为一块，分块求和即可。
### 实现：

```cpp
int ans=0;
for(int i=1,j;i<=n;i=j+1){
	j=n/(n/i);
	ans+=(j-(i-1))*(n/i);
}
```
# Mobius函数(莫比乌斯)
#### 定义：
$\mu(d)=\begin{cases} 1&d=1\\{(-1)^r}&{d=p_1...p_r,p_1,...,p_r}\text{是两两不同的素数}\\0&{\text{其他，即}d\text{有大于}1\text{的平方因数}}\end{cases}$
#### 性质：
（证明见课本）
（1）积性函数（$(d_1,d_2)=1$时$\mu(d_1d_2)=\mu(d_1)\mu(d_2)$）
（2）$\sum_{d|n}\mu(d)=[\frac{1}{n}]=\begin{cases} 1&{n=1}\\0&{n \neq 0}\end{cases}$
（3）$[gcd(i,j)=1]\iff\sum_{d|gcd(i,j)}\mu(d)$
*（3）的证明：通过Dirichlet卷积。
#### 筛法：
线性筛

```cpp
void getMu(){
	mu[1]=1;
	for(int i=2;i<=n;i++){
		if(!vis[i]) p[++tot]=i,mu[i]=-1;
		for(int j=1;j<=tot&&i*p[j]<=n;j++){
			vis[i*p[j]]=1;
			if(i%p[j]==0){
				mu[i*p[j]]=0;
				break;
			}
			mu[i*p[j]]=-mu[i];
		}
	}
}
```

# 杜教筛
被用来处理数论函数的前缀和问题。低于线性时间复杂度。
主要思路：通过寻找**递推关系**的方法来压缩计算步骤。

```cpp
ll get_large_mu(ll x){
	if(x<n) return sum_mu[x];//预处理（线性筛）计算的
	if(mp_mu.find(x)!=mp_mu.end()) return mp_mu[x];//get_large的时候计算过的
	ll res=1ll;
	for(ll i=2,j;i<=x;i=j+1){//从2开始
		j=x/(x/i);
		res-=get_large_mu(x/i)*(j-(i-1));//数论分块，递归调用
	}
	return mp_mu[x]=res;
}
```
# 线性同余方程
ax=c mod b
```c++
int ex_gcd(int a, int b, int& x, int& y) {
  if (b == 0) {
    x = 1;
    y = 0;
    return a;
  }
  int d = ex_gcd(b, a % b, x, y);
  int temp = x;
  x = y;
  y = temp - a / b * y;
  return d;
}
bool liEu(int a, int b, int c, int& x, int& y) {
  int d = ex_gcd(a, b, x, y);
  if (c % d != 0) return 0;
  int k = c / d;
  x *= k;
  y *= k;
  return 1;
}
```
# CRT
见doc
# 二次剩余
```c++
#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
int t;
ll n, p;
ll w;

struct num {  //建立一个复数域

  ll x, y;
};

num mul(num a, num b, ll p) {  //复数乘法
  num ans = {0, 0};
  ans.x = ((a.x * b.x % p + a.y * b.y % p * w % p) % p + p) % p;
  ans.y = ((a.x * b.y % p + a.y * b.x % p) % p + p) % p;
  return ans;
}

ll binpow_real(ll a, ll b, ll p) {  //实部快速幂
  ll ans = 1;
  while (b) {
    if (b & 1) ans = ans * a % p;
    a = a * a % p;
    b >>= 1;
  }
  return ans % p;
}

ll binpow_imag(num a, ll b, ll p) {  //虚部快速幂
  num ans = {1, 0};
  while (b) {
    if (b & 1) ans = mul(ans, a, p);
    a = mul(a, a, p);
    b >>= 1;
  }
  return ans.x % p;
}

ll cipolla(ll n, ll p) {
  n %= p;
  if (p == 2) return n;
  if (binpow_real(n, (p - 1) / 2, p) == p - 1) return -1;
  ll a;
  while (1) {  //生成随机数再检验找到满足非二次剩余的a
    a = rand() % p;
    w = ((a * a % p - n) % p + p) % p;
    if (binpow_real(w, (p - 1) / 2, p) == p - 1) break;
  }
  num x = {a, 1};
  return binpow_imag(x, (p + 1) / 2, p);
}
```
# Lucas 定理
```c++
long long Lucas(long long n, long long m, long long p) {
  if (m == 0) return 1;
  return (C(n % p, m % p, p) * Lucas(n / p, m / p, p)) % p;
}
```
# exLucas
```c++
LL CRT(int n, LL* a, LL* m) {
  LL M = 1, p = 0;
  for (int i = 1; i <= n; i++) M = M * m[i];
  for (int i = 1; i <= n; i++) {
    LL w = M / m[i], x, y;
    exgcd(w, m[i], x, y);
    p = (p + a[i] * w * x % mod) % mod;
  }
  return (p % mod + mod) % mod;
}
LL calc(LL n, LL x, LL P) {
  if (!n) return 1;
  LL s = 1;
  for (int i = 1; i <= P; i++)
    if (i % x) s = s * i % P;
  s = Pow(s, n / P, P);
  for (int i = n / P * P + 1; i <= n; i++)
    if (i % x) s = s * i % P;
  return s * calc(n / x, x, P) % P;
}
LL multilucas(LL m, LL n, LL x, LL P) {
  int cnt = 0;
  for (int i = m; i; i /= x) cnt += i / x;
  for (int i = n; i; i /= x) cnt -= i / x;
  for (int i = m - n; i; i /= x) cnt -= i / x;
  return Pow(x, cnt, P) % P * calc(m, x, P) % P * inverse(calc(n, x, P), P) %
         P * inverse(calc(m - n, x, P), P) % P;
}
LL exlucas(LL m, LL n, LL P) {
  int cnt = 0;
  LL p[20], a[20];
  for (LL i = 2; i * i <= P; i++) {
    if (P % i == 0) {
      p[++cnt] = 1;
      while (P % i == 0) p[cnt] = p[cnt] * i, P /= i;
      a[cnt] = multilucas(m, n, i, p[cnt]);
    }
  }
  if (P > 1) p[++cnt] = P, a[cnt] = multilucas(m, n, P, P);
  return CRT(cnt, a, p);
}
```
# Min_25
简单的函数
```c++
/* 「LOJ #6053」简单的函数 */
#include <algorithm>
#include <cmath>
#include <cstdio>

using i64 = long long;

constexpr int maxs = 200000;  // 2sqrt(n)
constexpr int mod = 1000000007;

template <typename x_t, typename y_t>
inline void inc(x_t &x, const y_t &y) {
  x += y;
  (mod <= x) && (x -= mod);
}
template <typename x_t, typename y_t>
inline void dec(x_t &x, const y_t &y) {
  x -= y;
  (x < 0) && (x += mod);
}
template <typename x_t, typename y_t>
inline int sum(const x_t &x, const y_t &y) {
  return x + y < mod ? x + y : (x + y - mod);
}
template <typename x_t, typename y_t>
inline int sub(const x_t &x, const y_t &y) {
  return x < y ? x - y + mod : (x - y);
}
template <typename _Tp>
inline int div2(const _Tp &x) {
  return ((x & 1) ? x + mod : x) >> 1;
}
template <typename _Tp>
inline i64 sqrll(const _Tp &x) {
  return (i64)x * x;
}

int pri[maxs / 7], lpf[maxs + 1], spri[maxs + 1], pcnt;

inline void sieve(const int &n) {
  for (int i = 2; i <= n; ++i) {
    if (lpf[i] == 0)
      pri[lpf[i] = ++pcnt] = i, spri[pcnt] = sum(spri[pcnt - 1], i);
    for (int j = 1, v; j <= lpf[i] && (v = i * pri[j]) <= n; ++j) lpf[v] = j;
  }
}

i64 global_n;
int lim;
int le[maxs + 1],  // x \le \sqrt{n}
    ge[maxs + 1];  // x > \sqrt{n}
#define idx(v) (v <= lim ? le[v] : ge[global_n / v])

int G[maxs + 1][2], Fprime[maxs + 1];
i64 lis[maxs + 1];
int cnt;

inline void init(const i64 &n) {
  for (i64 i = 1, j, v; i <= n; i = n / j + 1) {
    j = n / i;
    v = j % mod;
    lis[++cnt] = j;
    idx(j) = cnt;
    G[cnt][0] = sub(v, 1ll);
    G[cnt][1] = div2((i64)(v + 2ll) * (v - 1ll) % mod);
  }
}

inline void calcFprime() {
  for (int k = 1; k <= pcnt; ++k) {
    const int p = pri[k];
    const i64 sqrp = sqrll(p);
    for (int i = 1; lis[i] >= sqrp; ++i) {
      const i64 v = lis[i] / p;
      const int id = idx(v);
      dec(G[i][0], sub(G[id][0], k - 1));
      dec(G[i][1], (i64)p * sub(G[id][1], spri[k - 1]) % mod);
    }
  }
  /* F_prime = G_1 - G_0 */
  for (int i = 1; i <= cnt; ++i) Fprime[i] = sub(G[i][1], G[i][0]);
}

inline int f_p(const int &p, const int &c) {
  /* f(p^{c}) = p xor c */
  return p xor c;
}

int F(const int &k, const i64 &n) {
  if (n < pri[k] || n <= 1) return 0;
  const int id = idx(n);
  i64 ans = Fprime[id] - (spri[k - 1] - (k - 1));
  if (k == 1) ans += 2;
  for (int i = k; i <= pcnt && sqrll(pri[i]) <= n; ++i) {
    i64 pw = pri[i], pw2 = sqrll(pw);
    for (int c = 1; pw2 <= n; ++c, pw = pw2, pw2 *= pri[i])
      ans +=
          ((i64)f_p(pri[i], c) * F(i + 1, n / pw) + f_p(pri[i], c + 1)) % mod;
  }
  return ans % mod;
}

int main() {
  scanf("%lld", &global_n);
  lim = sqrt(global_n);

  sieve(lim + 1000);
  init(global_n);
  calcFprime();
  printf("%lld\n", (F(1, global_n) + 1ll + mod) % mod);

  return 0;
}
```
# Pollard-Rho
```c++
#include <bits/stdc++.h>

using namespace std;

typedef long long ll;
#define lll __int128

int t;
ll max_factor, n;

ll gcd(ll a, ll b) {
  if (b == 0) return a;
  return gcd(b, a % b);
}

ll quick_pow(ll x, ll p, ll mod) {
  ll ans = 1;
  while (p) {
    if (p & 1) ans = (lll)ans * x % mod;
    x = (lll)x * x % mod;
    p >>= 1;
  }
  return ans;
}

bool Miller_Rabin(ll p) {
  if (p < 2) return 0;
  if (p == 2) return 1;
  if (p == 3) return 1;
  ll d = p - 1, r = 0;
  while (!(d & 1)) ++r, d >>= 1;
  for (ll k = 0; k < 10; ++k) {
    ll a = rand() % (p - 2) + 2;
    ll x = quick_pow(a, d, p);
    if (x == 1 || x == p - 1) continue;
    for (int i = 0; i < r - 1; ++i) {
      x = (lll)x * x % p;
      if (x == p - 1) break;
    }
    if (x != p - 1) return 0;
  }
  return 1;
}

ll f(ll x, ll c, ll n) { return ((lll)x * x + c) % n; }

ll Pollard_Rho(ll x) {
  ll s = 0, t = 0;
  ll c = (ll)rand() % (x - 1) + 1;
  int step = 0, goal = 1;
  ll val = 1;
  for (goal = 1;; goal <<= 1, s = t, val = 1) {
    for (step = 1; step <= goal; ++step) {
      t = f(t, c, x);
      val = (lll)val * abs(t - s) % x;
      if ((step % 127) == 0) {
        ll d = gcd(val, x);
        if (d > 1) return d;
      }
    }
    ll d = gcd(val, x);
    if (d > 1) return d;
  }
}

void fac(ll x) {
  if (x <= max_factor || x < 2) return;
  if (Miller_Rabin(x)) {
    max_factor = max(max_factor, x);
    return;
  }
  ll p = x;
  while (p >= x) p = Pollard_Rho(x);
  while ((x % p) == 0) x /= p;
  fac(x), fac(p);
}

int main() {
  scanf("%d", &t);
  while (t--) {
    srand((unsigned)time(NULL));
    max_factor = 0;
    scanf("%lld", &n);
    fac(n);
    if (max_factor == n)
      printf("Prime\n");
    else
      printf("%lld\n", max_factor);
  }
  return 0;
}
```