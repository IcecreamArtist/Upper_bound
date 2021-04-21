# include <bits/stdc++.h>
# define fo(i,a,b) for(int i=(a);i<=(b);++i)
using namespace std;
struct _complex{
    double a,b;
    _complex(){}
    _complex(double a,double b): a(a),b(b){}
    _complex operator+(const _complex& y)const {
        return _complex(a+y.a,b+y.b);
    }
    _complex operator-(const _complex& y)const {
        return _complex(a-y.a,b-y.b);
    }
    _complex operator*(const _complex& y)const {
        return _complex(a*y.a-b*y.b,a*y.b-b*y.a);
    }
};
int main(){
    
    return 0;
}