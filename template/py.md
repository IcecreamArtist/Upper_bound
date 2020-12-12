```python
N=int(input())
s=input().split()
S=int(s[0])
T=int(s[1])
a=[]#一个列表
for i in range(1,N+1):
    k=input().split()
    a.append((int(k[0]),int(k[1])))
a.sort(key=lambda x:x[0]*x[1])
ans=0
for i in range(0,N):
    if(S//(a[i])[1]>ans):
        ans=S//(a[i])[1]
    S*=(a[i])[0]
print(ans)
```
```python
for _ in range(int(input())):
    a,b=map(int,input().split())
    print(a^b)
```
读到end of file：
```python
if __name__ == '__main__':
    while True:
        try:
            # do stuffs of the problem
        except EOFError:   #end of file退出
            break


```
