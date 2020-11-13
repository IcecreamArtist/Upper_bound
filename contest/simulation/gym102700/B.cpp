# include <bits/stdc++.h>
# define cl(x) memset(x,0,sizeof x)
using namespace std;
const int maxs = 2e5+10;
namespace gsa{
	int sa[maxs],rank[maxs],qzh[maxs];
	int tmpsa[maxs],tmpr[maxs],t[maxs];
	char s[maxs];
	inline bool same(int a,int b,int p){return t[a]==t[b]&&t[a+p]==t[b+p];}
	void getsa(int n,int m=233){	//m = |Sigma|, the size of sets of characters.
		cl(sa); cl(rank); cl(qzh); cl(tmpsa); cl(tmpr); cl(t);
		for(int i=0;i<n;++i) rank[i] = s[i],++qzh[rank[i]];
		for(int i=1;i<m;++i) qzh[i] += qzh[i-1];
		for(int i=n-1;i>=0;--i) sa[--qzh[s[i]]] = i;
		for(int j=1;j<=n;j<<=1){
			int cur = -1;
			for(int i=n-j;i<n;++i) tmpsa[++cur] = i;
			for(int i=0;i<n;++i) if(sa[i] >= j) tmpsa[++cur] = sa[i] - j;
			for(int i=0;i<n;++i) tmpr[i] = rank[tmpsa[i]];
			
			for(int i=0;i<m;++i) qzh[i] = 0;
			for(int i=0;i<n;++i)++qzh[tmpr[i]];
			for(int i=1;i<m;++i) qzh[i] += qzh[i-1];
			for(int i=n-1;i>=0;--i) t[i]=rank[i],sa[--qzh[tmpr[i]]]=tmpsa[i];
			m=0;
			for(int i=0;i<n;++i)
				rank[sa[i]]=(i&&same(sa[i],sa[i-1],j)?m:++m);
			++m;
		}
		for(int i=0;i<n;++i) rank[sa[i]] = i;
	}
}
char s1[maxs],s2[maxs];
int l1=0,l2=0;
int n;
int main(){
	scanf("%s",gsa::s);
	n = strlen(gsa::s); gsa::getsa(n+1);
	
	for(int i=gsa::sa[n];i<n;++i){
		s1[l1++] = gsa::s[i];
	}
	
	scanf("%s",gsa::s);
	n = strlen(gsa::s); gsa::getsa(n+1);
	for(int i=gsa::sa[n];i<n;++i){
		s2[l2++] = gsa::s[i];
	}
	
	int i;
	printf("%c",s1[0]);
	for(i=1;i<l1;++i){
		if(s1[i] < s2[0]) break;
		else printf("%c",s1[i]);
	}
	
	printf("%s\n",s2);
	return 0;
}
