# 高精度
// 原高精度算法的乘法有问题，那位作者两个板子都有问题，汗
# 快速乘法、快速幂
```c++
inline ll mul(ll a,ll b,ll p){  // 大数乘（防止乘的时候爆ll）
    ll ans = 0;
    while(b){
        if(b&1) ans = (ans+a)%p;
        a = (a+a)%p;
        b>>=1;
    }
    return ans;
}

ll qpow(ll a,ll n,ll p){
    ll ans = 1;
    while(n){if(n&1) ans=mul(ans,a,p);a=mul(a,a,p);n>>=1;}
    return ans;
}
```
