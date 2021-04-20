#include <bits/stdc++.h>
#define pb push_back
#define fi first
#define se second
#define ll long long
#define fo(i, a, b) for (int i = (a); i <= (b); ++i)
using namespace std;
const int maxn = 3e5 + 10;
const int maxs = 600;
int n, q;
int a[maxn];
int sum[maxs][maxn], cc[maxs][maxs];
int tmp[maxn];
int main()
{
    scanf("%d%d", &n, &q);
    for (int i = 0; i < n; ++i)
        scanf("%d", a + i);
    int bl = sqrt(n) + 1;
    int cnt = (n + bl - 1) / bl;
    for (int i = 1; i <= bl; ++i)
    {
        for (int j = 0; j < n; ++j)
            sum[i][j] = sum[i - 1][j];
        for (int j = (i - 1) * cnt; j < n && j < i * cnt; ++j)
            sum[i][a[j]]++;
    }
    for (int i = 1; i <= bl; ++i)
    {
        for (int j = i; j <= bl; ++j)
        {
            cc[i][j] = cc[i][j - 1];
            for (int k = (j - 1) * cnt; k < n && k < j * cnt; ++k)
            {
                int x = sum[j][a[k]] - sum[i - 1][a[k]];
                if (x > cc[i][j])
                {
                    cc[i][j] = x;
                }
            }
        }
    }
    while (q--)
    {
        int L, R;
        scanf("%d%d", &L, &R);
        --L;
        --R;
        int x = L / cnt + 1, y = R / cnt + 1;
        int ans = cc[x + 1][y - 1];
        for (int i = L; i <= R && i < x * cnt; i++)
            tmp[a[i]] = max(0, sum[y - 1][a[i]] - sum[x][a[i]]);
        for (int i = max(L, (y - 1) * cnt); i <= R; i++)
            tmp[a[i]] = max(0, sum[y - 1][a[i]] - sum[x][a[i]]);
        for (int i = L; i < x * cnt && i <= R; i++)
        {
            tmp[a[i]]++;
            if (tmp[a[i]] > ans)
            {
                ans = tmp[a[i]];
            }
        }
        if (x != y)
        {
            for (int i = (y - 1) * cnt; i <= R; i++)
            {
                tmp[a[i]]++;
                if (tmp[a[i]] > ans)
                {
                    ans = tmp[a[i]];
                }
            }
        }
        if (ans <= (R - L + 1 + 1) / 2)
            printf("1\n");
        else
            printf("%d\n", ans * 2 - (R - L + 1));
    }
    return 0;
}