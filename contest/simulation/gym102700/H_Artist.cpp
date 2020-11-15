/*
 * 题意：求一个字符串中有多少个本质不同回文子串，且这些子串的长度为奇数
 * 思路：用manacher求出所有回文子串然后插入hash统计
 */
#include <bits/stdc++.h>
using namespace std;
const int N = 1e5+6;
typedef unsigned long long ull;
int n;
ull seed = 31,base[N],_hash[N];
char s[N];
map<ull,int>mp;
int ans = 0;

void init(){
    base[0] = 1;
    for(int i=1;i<N;++i) base[i] = base[i-1]*seed;
}

void makehash(int len,char str[]){
    for(int i=1;i<=len;++i) _hash[i] = _hash[i-1]*seed+(str[i]-'a'+1);
}

// 从i开始，长度为l
ull gethash(int i,int l){
    return _hash[i+l-1]-_hash[i-1]*base[l];
}

void insert(int i,int l){
    if(l==1) return;
    ull tmp = gethash(i,l);
    if(mp.find(tmp)==mp.end()) mp[tmp]=1,ans++;
}

int d1[N];

void manacher_odd(char str[]){
    for(int i=1,l=1,r=0;i<=n;++i){
        int k = (i>r)?1:min(d1[l+r-i],r-i);
        while(1<=i-k&&i+k<=n&&str[i-k]==str[i+k]) k++;
        d1[i] = k--;
        cout<<d1[i]<<endl;
        if(i+k>r) l=i-k,r=i+k;
    }
}

int main(){
    scanf("%d",&n);
    scanf("%s",s+1);
    init();
    makehash(n,s);
    manacher_odd(s);
    printf("%d\n",ans);
}