# 树
## 基础
### 树的大小
```c++
void dfs(int u){
	sz[u] = 1;
	for(auto i:G[u]) if(!sz[i]) dfs(i),sz[u] += sz[i];
}
```
### 树的直径
```c++
vector<int>G[maxn];
int dep[maxn],f[maxn],in[maxn];
int dia,d;

void dfs(int u,int fa){
	dep[u]=dep[f[u]=fa]+1;
	if(dep[u]>dia) dia = dep[u],d = u;
	for(auto i:G[u]) if(i-fa) dfs(i,u);
}

int main(){
	dfs(1,0);
	x = d;
	dia = 0,dfs(x,0);
	y = d;
	// dia就是答案
}
```
### 