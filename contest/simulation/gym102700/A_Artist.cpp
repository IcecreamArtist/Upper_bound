#include<bits/stdc++.h>
using namespace std;
double ax,ay,bx,by,cx,cy,dx,dy,len1,len2;
const double eps = 1e-12;     // 精度
/*
 * 思路：三分
 * 模拟过程，分两段，先是三分从A点走到B点（要使得前面的路径更短）。然后再站定在A点，三分另一个人继续走。
 * 其实第一段可以直接判断（1）是否交叉=0（2）只考虑两个起点的连线与两个终点的连线。但由于这样也符合为一个抛物线。
 * 为减少讨论，用三分处理。注意精度问题。
 */

double get_dis(double x1,double y1,double x2,double y2){return (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1);}

double check(double dis){
    // 计算长边上的点距离起点走了dis，两点的距离
    double x1,y1,x2,y2;  // 两点的坐标
    if(dis>=len1-eps) x1=bx,y1=by;
    else x1 = ax+(bx-ax)/len1*dis,y1 = ay+(by-ay)/len1*dis;
    if(dis>=len2-eps) x2=dx,y2=dy;
    else x2 = cx+(dx-cx)/len2*dis,y2 = cy+(dy-cy)/len2*dis;
    return  get_dis(x1,y1,x2,y2);
}
// 三分板子
double three_div(double left,double right){
    // left and right represent the distance gone
    double ans=1e25;
    for(int i=1;i<=1e5;++i){
        double lmid = (2*left+right)/3;
        double rmid = (2*right+left)/3;
        double resl = check(lmid),resr = check(rmid);
        if(resl+eps<=resr) right = rmid;
        else left = lmid;
        ans = min({ans,resl,resr});
    }
    return ans;
}

int main(){
    scanf("%lf%lf%lf%lf",&ax,&ay,&bx,&by);
    scanf("%lf%lf%lf%lf",&cx,&cy,&dx,&dy);
    len1 = sqrt(get_dis(ax,ay,bx,by)),len2 = sqrt(get_dis(cx,cy,dx,dy));
    double ans = min(three_div(0,min(len1,len2)),three_div(min(len1,len2),max(len1,len2)));
    printf("%.12lf\n",sqrt(ans));
}