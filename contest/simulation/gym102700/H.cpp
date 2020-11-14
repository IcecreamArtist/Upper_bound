#include<bits/stdc++.h>
using namespace std;
const int MAX=2e6+10;
const int Hashsize=2000003;
const unsigned long long p=131;
typedef long long ll;
typedef unsigned long long ull;
struct lenka
{
    int next;
    ull val;
}ed[MAX];
int head[MAX],cnt=0;
ull f[MAX],sum[MAX];
int ans=0;
void Insert(int x,int y)
{
    ull tot=sum[y]-sum[x-1]*f[y-x+1];
    for(int i=head[tot%Hashsize];i!=-1;i=ed[i].next)
    {
        if(tot==ed[i].val)return;
    }
    ans++;
    ed[cnt].next=head[tot%Hashsize];
    ed[cnt].val=tot;
    head[tot%Hashsize]=cnt++;
}
char s[MAX];
int len[MAX];
int main()
{
	int n; scanf("%d",&n);
    scanf("%s",s+1);
    n=strlen(s+1);
    f[0]=1;
    for(int i=1;i<=n;i++)
    {
        f[i]=f[i-1]*p;
        sum[i]=sum[i-1]*p+s[i];
    }
    memset(head,-1,sizeof head);
    cnt=0;
    int mx=0,x=0;
    for(int i=1;i<=n;i++)
    {
    //    Insert(i,i);
        if(mx>i)len[i]=min(mx-i,len[2*x-i]);
        while(i+len[i]+1<=n&&s[i+len[i]+1]==s[i-len[i]-1])
        {
            Insert(i-len[i]-1,i+len[i]+1);
            len[i]++;
        }
        if(i+len[i]>mx)
        {
            mx=i+len[i];
            x=i;
        }
    }
    mx=x=0;
    memset(len,0,sizeof len);
    memset(head,-1,sizeof head);
    cnt=0;
    /*
    for(int i=2;i<=n;i++)
    {
        if(mx>i)len[i]=min(mx-i+1,len[2*x-i]);
        while(i+len[i]<=n&&s[i+len[i]]==s[i-len[i]-1])
        {
            Insert(i-len[i]-1,i+len[i]);
            len[i]++;
        }
        if(i+len[i]-1>mx)
        {
            mx=i+len[i]-1;
            x=i;
        }
    }
    */
    printf("%d\n",ans);
    return 0;
}
