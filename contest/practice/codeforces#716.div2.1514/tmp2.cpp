#include <cstdio>
#include <algorithm>
#include <cmath>
#include <cstring>
#define LL long long
using namespace std;
const int MAXN = 1e5 + 10, INF = 1e9 + 7;
inline int read()
{
    char c = getchar();
    int x = 0, f = 1;
    while (c < '0' || c > '9')
    {
        if (c == '-')
            f = -1;
        c = getchar();
    }
    while (c >= '0' && c <= '9')
        x = x * 10 + c - '0', c = getchar();
    return x * f;
}
int N, Q;
LL out[MAXN], ans;
int be[MAXN], date[MAXN], cnt[MAXN], a[MAXN], tot, base, num;
struct Query
{
    int l, r, id;
    bool operator<(const Query &rhs) const
    {
        return be[l] == be[rhs.l] ? r < rhs.r : be[l] < be[rhs.l];
    }
} q[MAXN];
LL solve(int l, int r)
{
    static int tim[MAXN];
    LL ans = 0;
    for (int i = l; i <= r; i++)
        tim[a[i]] = 0;
    for (int i = l; i <= r; i++)
        tim[a[i]]++, ans = max(ans, 1ll * tim[a[i]] * date[a[i]]);
    return ans;
}
void Add(int x)
{
    cnt[a[x]]++;
    ans = max(ans, 1ll * cnt[a[x]] * date[a[x]]);
}
void Del(int x)
{
    cnt[a[x]]--;
}
int Get(int i, int id)
{
    int R = min(N, id * base), ql = R + 1, qr = ql - 1;
    ans = 0;
    memset(cnt, 0, sizeof(cnt));
    for (; be[q[i].l] == id; i++)
    {
        if (be[q[i].l] == be[q[i].r])
        {
            out[q[i].id] = solve(q[i].l, q[i].r);
            continue;
        }
        while (qr < q[i].r)
            Add(++qr);
        LL cur = ans;
        while (ql > q[i].l)
            Add(--ql);
        out[q[i].id] = ans;
        while (ql < R + 1)
            Del(ql++); //每次询问完之后重新统计答案
        ans = cur;
    }
    return i;
}
main()
{
    //freopen("4241.in", "r", stdin);
    //freopen("4241.out", "w", stdout);
    N = read();
    Q = read();
    base = sqrt(N);
    for (int i = 1; i <= N; i++)
    {
        a[i] = date[i] = read();
        be[i] = (i - 1) / base + 1;
        num = max(num, be[i]);
    }

    sort(date + 1, date + N + 1);
    int tot = unique(date + 1, date + N + 1) - date - 1;
    for (int i = 1; i <= N; i++)
        a[i] = lower_bound(date + 1, date + N + 1, a[i]) - date;

    for (int i = 1; i <= Q; i++)
        q[i].l = read(), q[i].r = read(), q[i].id = i;
    sort(q + 1, q + Q + 1);

    for (int i = 1, id = 1; id <= num; id++)
        i = Get(i, id);

    for (int i = 1; i <= Q; i++)
        printf("%lldn", out[i]);
    return 0;
}