# include <bits/stdc++.h>
# define eps 1e-6
using namespace std;
struct Point{
	//It's a point and a vector!!!
	double x,y;
	Point(){}
	Point(double x,double y):x(x),y(y){}
	Point operator+(Point t)const{return Point(x+t.x,y+t.y);}
	Point operator-(Point t)const{return Point(x-t.x,y-t.y);}
	Point operator*(double v)const{return Point(v*x,v*y);}
};
inline double q2(double x){return x*x;}
inline double dist2(Point a,Point b){return sqrt(q2(a.x-b.x)+q2(a.y-b.y)); }
Point A,B,C,D;
double d1,d2;
inline double calc(double t){
	Point a=(t>=d1?B:A+(B-A)*(t/d1));
	Point b=(t>=d2?D:C+(D-C)*(t/d2));
	return dist2(a,b);
}
//Line AB & CD
int main(){
//	freopen("A.in","r",stdin);
	cin >> A.x >> A.y;
	cin >> B.x >> B.y;
	cin >> C.x >> C.y;
	cin >> D.x >> D.y;
	
	d1=dist2(A,B);
	d2=dist2(C,D);
	//now use trinary search to get the minimum dist of all
	double res1,res2;
	double l=0,r=min(d1,d2);
	while(r-l>eps){
		double lm=(l+l+r)/3,rm=(l+r+r)/3;
		double lv = calc(lm),rv=calc(rm);
	//	cout << l << ':' << lv << ' ' << r << ':' << rv << endl;
		if(lv<rv){
			r = rm;
		}
		else{
			l = lm;
		}
	}
	res1 = l;
	
	l=min(d1,d2),r=max(d1,d2);
	while(r-l>eps){
		double lm=(l+l+r)/3,rm=(l+r+r)/3;
		double lv = calc(lm),rv=calc(rm);
	//	cout << l << ':' << lv << ' ' << r << ':' << rv << endl;
		if(lv<rv){
			r = rm;
		}
		else{
			l = lm;
		}
	}
	res2 = l;
	printf("%.12lf\n",min(calc(res1),calc(res2)));
	return 0;
}
