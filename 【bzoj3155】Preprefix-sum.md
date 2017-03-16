---
title: 【bzoj3155】Preprefix sum
date: 2017-02-14 14:06:48
tags:
  - 线段树
  - 树状数组
categories:
  - 数据结构
---

> 这是一道没什么与意思的板子题。

<!--more-->

# 题目大意
两个操作：① 修改一个位置的权值； ② 查询前缀和的前缀和；
# 解题报告
做的一道**弱智题**；

1. 可以用**线段树**维护前缀和序列，$a_i$ 对前缀和序列的影响是在$[i, n]$加上$a_i$, 就是一个区间加，查询时候就是一个区间和；
2. 查询$S_{S_x}$, 展开看就是 $\Sigma_{i=1}^{x} a_i * (n - i + 1)$,也就是$(x + 1) * \Sigma_{i=1}^{x} a_i - \Sigma_{i=1}^{x} a_i * i$;
用**bit**分别维护$a_i$和$a_i * i$的前缀和；

# 代码
1. 线段树
```c++
#include <iostream> 
#include <cstdio> 
#include <cstring> 
#include <cstdlib> 
#include <vector> 
using namespace std; 
#define for1(i, n) for (i = 1; i <= n; ++i) 
#define for0(i, n) for (i = 0; i < n; ++i) 
#define mp make_pair
#define pb push_back 
#define F first
#define S second
typedef pair<int, int> pii; 
typedef vector<int> vi; 
typedef long long ll; 
const int N = 101000; 
//string s; 
char s[20];
int n, x, m, a[N]; 
ll sm[N << 2], tag[N << 2]; 
inline void in(int &x) { 
    char ch = getchar(); 
    for (; ch < '0' || ch > '9'; ch = getchar()); 
    for (x = 0; ch >= '0' && ch <= '9'; ch = getchar())
        x = x * 10 + ch - 48; 
}
inline void down(int x, int l, int r) { 
    if (tag[x]) { 
        int mid = l + r >> 1; 
        sm[x << 1] += 1LL * (mid - l + 1) * tag[x]; 
        sm[x << 1 | 1] += 1LL * (r - mid) * tag[x]; 
        tag[x << 1] += tag[x], tag[x << 1 | 1] += tag[x]; 
        tag[x] = 0; 
    }
}
inline void up(int x) { 
    sm[x] = sm[x << 1] + sm[x << 1 | 1]; 
}    
void add(int x, int l, int r, int L, int R, int delta) { 
    if (L <= l && r <= R) {
        sm[x] += 1LL * (r-l+1) * delta; 
        tag[x] += delta; 
    } else {
        int mid = (l + r) >> 1; down(x, l, r); 
        if (L <= mid) add(x << 1, l, mid, L, R, delta); 
        if (R > mid) add(x << 1 | 1, mid + 1, r, L, R, delta); 
        up(x); 
    }
}
ll sum(int x, int l, int r,int L, int R) { 
    if (L <= l && r <= R) return sm[x]; 
    else { 
        int mid = (l + r) >> 1; ll tmp = 0; down(x, l, r); 
        if (L <= mid) tmp += sum (x << 1, l, mid, L, R); 
        if (R > mid) tmp += sum(x << 1 | 1, mid + 1, r, L, R); 
        return tmp ; 
    }
}
inline void modify(int w, int x) { 
    add(1, 1, n, w, n, x - a[w]), a[w] = x; 
}
inline void query(int w) { 
    printf("%lld\n", sum(1, 1, n, 1, w)); 
}
int main() { 
    in(n), in(m); int i, j; 
    for1(i, n) in(x), modify(i, x); 
    for1(i, m) { 
        scanf("%s", s); 
        if (s[0] == 'Q') in(x), query(x); 
        if (s[0] == 'M') in(j), in(x), modify(j, x); 
    }
    return 0; 
}
    
```
2. bit
```c++
#include <iostream> 
#include <cstdio> 
#include <cstdlib> 
#include <cstring> 
#include <algorithm> 
using namespace std; 
const int N = 100100; 
typedef long long ll; 
ll sm[N], bit[N]; 
int  n, a[N], m ;
char s[20]; 
inline void in(int &x) { 
    char ch = getchar(); 
    for (; ch < '0' || ch > '9'; ch = getchar());   
    for (x = 0; ch >= '0' && ch <= '9'; ch = getchar())
        x =  x * 10 + ch - 48; 
}
inline void add(ll *a, int x, ll val) { 
    for (; x <= n; x += x & -x) a[x] += val; 
}
inline ll sum(ll *a, int x) { 
    ll tmp = 0; 
    for (; x; x -= x & -x) tmp += a[x]; 
    return tmp; 
}
inline void modify(int place, int x) {
    add(sm, place, x - a[place]), add(bit, place, 1LL * (x - a[place]) * place);
    a[place] = x; 
}
inline void query(int x) { 
//    cout << x << ": " << sum(sm, x) << ' ' << sum(bit, x) << endl; 
    ll tmp = 1LL * sum(sm, x) * (x + 1) - sum(bit, x); 
    printf("%lld\n", tmp); 
}
int main() { 
//    freopen("A.in", "r", stdin); 
//    freopen("A.out", "w", stdout); 
    in(n), in(m); int i, j, x; 
    for (i = 1; i <= n; ++i) 
        in(x), modify(i, x); 
    for (i = 1; i <= m; ++i) {
        scanf("%s", s); 
        if (s[0] == 'Q') in(x), query(x); 
        if (s[0] == 'M') in(j), in(x), modify(j, x); 
    }
    return 0;  
}
```
# 题目
[题目链接][1]

  [1]: http://www.lydsy.com/JudgeOnline/problem.php?id=3155

<center><h2>3155: Preprefix sum</h2><span class="green">Time Limit: </span>1 Sec&nbsp;&nbsp;<span class="green">Memory Limit: </span>512 MB<br></center><h2>Description</h2><div class="content"><p><img src="http://www.lydsy.com/JudgeOnline/upload/201503/222(2).png" width="719" height="453" alt=""></p></div><h2>Input</h2><div class="content"><p class="MsoNormal" style="margin: 0cm 0cm 10pt; layout-grid-mode: both; mso-margin-top-alt: auto; mso-margin-bottom-alt: auto; mso-layout-grid-align: auto"></p>
<div>第一行给出两个整数N，M。分别表示序列长度和操作个数</div>
<div>接下来一行有N个数，即给定的序列a1,a2,....an</div>
<div>接下来M行，每行对应一个操作，格式见题目描述</div>
<p></p></div><h2>Output</h2><div class="content"><p class="MsoNormal" style="margin: 0cm 0cm 10pt; layout-grid-mode: both; mso-margin-top-alt: auto; mso-margin-bottom-alt: auto; mso-layout-grid-align: auto"><span style="font-size: 12pt; font-family: 宋体; mso-ascii-font-family: 'Times New Roman'">对于每个询问操作，输出一行，表示所询问的</span><span lang="EN-US" style="font-size: 12pt; font-family: 'Times New Roman'; mso-fareast-font-family: 宋体; mso-no-proof: yes">SSi</span><span style="font-size: 12pt; font-family: 宋体; mso-ascii-font-family: 'Times New Roman'">的值。</span><span lang="EN-US" style="font-size: 12pt; font-family: 'Times New Roman'; mso-fareast-font-family: 宋体"><o:p></o:p></span></p>
<p class="MsoNormal" style="margin: 0cm 0cm 0pt; layout-grid-mode: both; mso-layout-grid-align: auto; tab-stops: 45.8pt 91.6pt 137.4pt 183.2pt 229.0pt 274.8pt 320.6pt 366.4pt 412.2pt 458.0pt 503.8pt 549.6pt 595.4pt 641.2pt 687.0pt 732.8pt"><span lang="EN-US" style="font-size: 12pt; font-family: 'Times New Roman'; mso-fareast-font-family: 宋体"><o:p></o:p></span></p></div><h2>Sample Input</h2>
			<div class="content"><span class="sampledata">5 3<br>
1 2 3 4 5<br>
Query 5<br>
Modify 3 2<br>
Query 5</span></div><h2>Sample Output</h2>
			<div class="content"><span class="sampledata">35<br>
32</span></div><h2>HINT</h2>
