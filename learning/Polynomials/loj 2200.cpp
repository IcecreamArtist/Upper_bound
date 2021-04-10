#include <bits/stdc++.h>
using namespace std;
const double PI = acos(-1.0);
struct Complex {
    double x, y;
    Complex(double _x = 0.0, double _y = 0.0) {
        x = _x;
        y = _y;
    }
    Complex operator-(const Complex &b) const {
        return Complex(x - b.x, y - b.y);
    }
    Complex operator+(const Complex &b) const {
        return Complex(x + b.x, y + b.y);
    }
    Complex operator*(const Complex &b) const {
        return Complex(x * b.x - y * b.y, x * b.y + y * b.x);
    }
};

void change(Complex y[], int len) {
    int i, j, k;
    for (int i = 1, j = len / 2; i < len - 1; i++) {
        if (i < j)
            swap(y[i], y[j]);
        k = len / 2;
        while (j >= k) {
            j -= k;
            k /= 2;
        }
        if (j < k)
            j += k;
    }
}
void fft(Complex y[], int len, int on) {
    change(y, len);
    for (int h = 2; h <= len; h <<= 1) {
        Complex wn(cos(2 * PI / h), sin(on * 2 * PI / h));
        for (int j = 0; j < len; j += h) {
            Complex w(1, 0);
            for (int k = j; k < j + h / 2; k++) {
                Complex u = y[k];
                Complex t = w * y[k + h / 2];
                y[k] = u + t;
                y[k + h / 2] = u - t;
                w = w * wn;
            }
        }
    }
    if (on == -1) {
        for (int i = 0; i < len; i++) {
            y[i].x /= len;
        }
    }
}
const int maxn = 400020;
Complex x1[maxn], x2[maxn], x3[maxn];

int main() {
    int n;
    scanf("%d", &n);
    for (int i = 0; i < n; i++) {
        scanf("%lf", &x1[i].x);
        x3[n - 1 - i].x = x1[i].x;
    }
    for (int i = 1; i < n; i++)
        x2[i] = 1.0 / i / i;
    for (int i = 1; i < n; i++)
        cout << "x:" << x2[i].x << " y:" << x2[i].y << endl;
    int len = 1;
    while (len < n * 2)
        len <<= 1;
    fft(x1, len, 1);
    fft(x2, len, 1);
    fft(x3, len, 1);
    for (int i = 0; i < len; i++)
        x1[i] = x1[i] * x2[i], x3[i] = x3[i] * x2[i];
    fft(x1, len, -1);
    fft(x3, len, -1);
    for (int i = 0; i < n; i++)
        printf("%lf\n", x1[i].x - x3[n - i - 1].x);
    return 0;
}
