# include <bits/stdc++.h>
using namespace std;
struct point{
	double x,y;
	point(){}
	point(double x,double y):x(x),y(y){}
	point operator+(const point& a)const{return point(x+a.x,y+a.y);}
	point operator-(const point& a)const{return point(x-a.x,y-a.y);}
	bool operator<(const point& a)const{return x==a.x?y<a.y:x<a.x;}
	double cross(const point& a,const point& b){
		return a.x*b.y-a.y*b.x;
	}
	double dot(const point& a,const point& b){
		return a.x*b.x+a.y*b.y;
	}
} p[maxn];
int n;
int convex[maxn],tot1,tot2;
int main(){
	scanf("%d",&n);
	for(int i=1;i<=n;++i) scanf("%lf%lf",&p[i].x,&p[i].y);
	
	//Convex Hull
	sort(p+1,p+n+1);
	//1...tot1 is the lower hull
	//tot1 ... tot2 is the upper hull
	tot1 = 0;
	for(int i=1;i<=n;++i){
		while(tot>2&&cross(p[convex[tot1]]-p[convex[tot1-1]],p[i]-p[convex[tot1]]) <= 0)--tot1;
		convex[++tot1]=i;
	}
	
	tot2 = tot1;
	for(int i=n-1;i>0;--i){
		while((tot2 > tot1)&&dot(p[convex[tot2]-p[convex[tot2-1]],p[i]-p[convex[tot2]]) <= 0) --tot2;
		convex[++tot2]=i;
	}
	
}