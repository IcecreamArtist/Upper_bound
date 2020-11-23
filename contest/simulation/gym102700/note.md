## Problem C
设$f(x)$表示，长度为$x$的，不是“循环串”的串的个数，那么设$F(x)$表示长度为$x$的所有串的个数，那么
$$|\Sigma|^x=F(x)=\sum_{d|x} f(d)$$
$$f(x)=\sum_{d|x}\mu(\frac{x}{d})F(d)=\sum_{d|x}\mu(\frac{x}{d})|\Sigma|^d$$
然后我们要求$\sum_{i=1}^nf(i)$，可以发现就是$\sum_{d=1}^n|\Sigma|^d\sum_{x=1}^n[d|x]\mu(\frac{x}{d})=\sum_{d=1}^n|\Sigma|^d\sum_{j=1}^{[n/d]}\mu(j)$