---
title: 【bzoj3309】DZY Loves Math
date: 2017-02-14 14:17:46
tags:
  - 筛法
  - 反演
categories:
  - 数学题
---

> 龙队长和yveh屠杀无辜数论题的帮凶...

<!--more-->

# 题目大意
定义 $f(n)$ 为 $n$ 所含质因子的最大幂指数, 给定正整数 $a,b$ ，求 $\displaystyle \sum_{i=1}^{a} \sum_{j=1}^{b} f(gcd(i, j))$
# 解题报告
首先化简式子...设 $a \leqslant  b$ ,
$$\displaystyle \sum_{i=1}^{a} \sum_{j=1}^{b} f(gcd(i, j)) = \sum_{i=1}^{a}f[i] \sum_{x=1}^{\lfloor \frac{a}{i} \rfloor} \sum_{y=1}^{\lfloor \frac{b}{i} \rfloor}[(x, y)=1]$$
$$= \sum_{i=1}^{a}f[i] \sum_{d=1}^{\lfloor \frac{a}{i} \rfloor} \mu(d) \lfloor \frac{a}{id} \rfloor \lfloor \frac{b}{id} \rfloor$$
$$= \sum_{k=1}^{n} \sum_{d|k} \mu(d) f(\frac{k}{d}) \lfloor \frac{a}{k} \rfloor \lfloor \frac{b}{k} \rfloor$$
上面这个式子, $\lfloor \frac{a}{k} \rfloor \lfloor \frac{b}{k} \rfloor$ 的取值只有 $4 \sqrt{n}$ 个, 所以需要筛出 $\mu \times  f$ , 就可以得到 $O(q \sqrt{a})$ 的做法; <br>
首先, 利用 $\mu$ 的性质, $\mu(p^k, k > 1) = 0$ , 所以, 产生贡献的 $\mu(d)f(n/d)$ 中, $d$ 的每个质因子的幂指数为 $0, 1$ ;   
考虑 $f(d)$ 的意义是最大幂指数, 所以将 $n$ 质因数分解得到 $\prod p_i^{k_i}$ , 设最大幂指数为 $k_m$ , 则质因数可以分为幂指数为 $k_m$ 的, 和幂指数小于 $k_m$ 的...<br>
对 $\mu \times f$ 产生贡献的 $k$ , 显然, 是对每个质因数进行不选或者选一个的抉择, 不妨分为在幂指数为 $k_m$ 的质因数中做抉择和在其他质因数中做抉择两部分.<br>
利用二项式定理, $\displaystyle \sum_{i=0}^{n} (-1)^{i}C_n^i = 0$, 不论在幂指数为 $k_m$ 的质因数中做怎样的抉择, 在其他质因数中作抉择形成的 $\mu$ 会使得贡献为 $0$ , 所以得到重要结论,  $u \times f$ 不为 $0$ 的数, 所有质因子的幂指数相同;<br>
所以现在只关心幂指数相同的质因子组成的数, 仍然根据二项式定理, 如果不论怎么选 $d$ , $f(n/d) = k_m$ 的话, $(\mu \times f)(n) = 0$ , 但是, 如果 $d$ 选择所有的质因子, $f(n/d) = k_m - 1$ , 根据二项式系数的正负, 如果质因子个数为奇数, $(\mu \times f)(n) = 1$ , 否则 $(\mu \times f)(n) = -1$ ;<br>
现在可以筛 $(\mu \times f)$ 了, yveh 和 龙队分别筛了 $5$ 和 $3$ 个量, 我不是很会他们的做法...我的做法是先只筛 $k_m = 1$ 的数, 然后它的幂的 $f$ 与它相同, 直接赋值.. 复杂度是 $O(n + sq(n))$ , $sq(n)$ 是 $[1, n]$ 中的完全 $k$ 次方数;<br>
有了 $\mu \times f$ , 可以暴力分块求答案了..
![激烈鼓掌](1.jpg)

# 代码

```c++
#include <bits/stdc++.h>

using namespace std;

#define fi first
#define se second
#define mp make_pair
#define pb push_back
#define FORU(i, a, b) for (int i = a, nn = int(b); i <= nn; ++i)
#define FORD(i, a, b) for (int i = a, nn = int(b); i >= nn; --i)
#define REP(i, b) for (int i = 0, nn = int(b); i < nn; ++i)

typedef long long ll;

const int N = 10000001;

int f[N], n, prime[2000100], test, pn;
bool no_prime[N];
char *cp = (char *) malloc(1000000);
char *os = (char *) malloc(1000000), *ot = os;
inline void in(int &x) {
    int f = 1;
    for (; *cp < '0' || *cp > '9'; cp ++)
        if (*cp == '-') f = -1;
    for (x = 0; *cp >= '0' && *cp <= '9'; cp ++)
        x = x * 10 + *cp - 48;
    x *= f;
}
inline void out(ll x) {
    if (x) out(x / 10), *ot ++ = x % 10 + '0';
}

inline void print(ll &x) {
    if (x == 0) *ot++ = '0';
    else out(x);
    *ot ++ = '\n';
}
inline void euler() {
    f[1] = 0, no_prime[1] = 1;
    FORU(i, 2, N - 1) {
        if (!no_prime[i]) {
            prime[++pn] = i;
            f[i] = 1;
        }
        for (int j = 1; prime[j] * i < N && j <= pn; ++j) {
            no_prime[i * prime[j]] = 1;
            if (i % prime[j]) {
                f[i * prime[j]] = -f[i];
            } else break;
        }
    }
    FORU(i, 1, 3162)
        if (f[i]) for (ll j = 1LL * i * i; j < N; j *= i)
            f[j] = f[i];
    FORU(i, 2, N - 1)
        f[i] += f[i - 1];
}
inline void query(int a, int b) {
    if (a > b) swap(a, b);
    ll ans = 0;
    for (int i = 1, j; i <= a; i = j + 1) {
        j = min(a/(a/i), b/(b/i));
        ans += 1LL * (f[j] - f[i-1]) * (a/i) * (b/i);
    }
    print(ans);
}

int main() {
    fread(cp, 1, 1000000, stdin);
    euler();
    in(test);  int a, b;
    while (test --) {
        in(a), in(b);
        query(a, b) ;
    }
    fwrite(os, 1, ot - os, stdout);
    return 0;
}
```

# 原题信息

[题目链接](http://www.lydsy.com/JudgeOnline/problem.php?id=3309)

<body>


<title>Problem 3309. -- DZY Loves Math</title><center><h2>3309: DZY Loves Math</h2><span class="green">Time Limit: </span>20 Sec&nbsp;&nbsp;<span class="green">Memory Limit: </span>512 MB<br><span class="green">Submit: </span>782&nbsp;&nbsp;<span class="green">Solved: </span>418<br>[<a href="submitpage.php?id=3309">Submit</a>][<a href="problemstatus.php?id=3309">Status</a>][<a href="bbs.php?id=3309">Discuss</a>]</center><h2>Description</h2><div class="content"><p><span style="font-size: medium">对于正整数n，定义f(n)为n所含质因子的最大幂指数。例如f(1960)=f(2^3 * 5^1 * 7^2)=3, f(10007)=1, f(1)=0。<br>
给定正整数a,b，求sigma(sigma(f(gcd(i,j)))) (i=1..a, j=1..b)。</span></p>
<p></p></div><h2>Input</h2><div class="content"><p><span style="font-size: medium">第一行一个数T，表示询问数。<br>
接下来T行，每行两个数a,b，表示一个询问。</span></p>
<p></p></div><h2>Output</h2><div class="content"><p><span style="font-size: medium">对于每一个询问，输出一行一个非负整数作为回答。</span></p>
<p></p></div><h2>Sample Input</h2>
			<div class="content"><span class="sampledata">4<br>
10000<br>
7558588 9653114<br>
6514903 4451211<br>
7425644 1189442<br>
6335198 4957<br>
</span></div><h2>Sample Output</h2>
			<div class="content"><span class="sampledata">35793453939901<br>
14225956593420<br>
4332838845846<br>
15400094813</span></div><h2>HINT</h2>
			<div class="content"><p></p><p>【数据规模】<br><br>
T&lt;=10000<br><br>
1&lt;=a,b&lt;=10^7</p><br>
<p></p><p></p></div><h2>Source</h2>
</body>

