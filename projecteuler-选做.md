---
title: projecteuler-选做
date: 2017-03-14 10:05:51
tags:
  - 筛法
  - 数论
  - 枚举
categories:
  - 题目集锦
---

> 通过做project euler上的题目， 做一些数学方面的练习

<!--more-->

# [Project Euler #439](https://projecteuler.net/problem=439)
## 题目大意
$$\displaystyle S(N) = \sum_{1 \leqslant i \leqslant N} \sum_{1 \leqslant j \leqslant N} d(i * j)$$
$$d(n) = \sum _{k|n} k$$
## 解题报告
本来想新年第一天搞一发, 但是可怜的是当天没调出来; <br>
首先将式子化简:
$$d(i*j) = \sum \frac{j * p}{q} (p|i,q|j,(p,q)=1) $$
$$S(n) = \sum_{i=1}^{n} \sum_{j=1}^{n}[(p,q)=1] \ p * \sum _{i*p \leqslant n} \sum_{j * q \leqslant n} j$$
$$S(n) = \sum_{p=1}^{n} \sum_{q=1}^{n} \sum_{d|(p,q)} \mu (d)\ p * \sum_{i * p \leqslant n} \sum_{j * q \leqslant n} j$$
$$S(n) = \sum_{d=1}^{n} \mu(d)d \sum_{pd \leqslant n} p\sum_{i * pd \leqslant n} \  \sum_{qd \leqslant n}\sum_{j * qd \leqslant n}j$$
如果把式子中的 $\lfloor \frac{n}{d} \rfloor$ 设为 $I$ , 则
$$S(n) = \sum_{d=1}^{n} \mu(d)d (\sum_{i}^I\sum_{d|i}d)^2$$
现在就可以做了, 前半部分是我刚学习的杜教筛, 后半部分只有 $\sqrt{n}$ 个取值, 是约数和可以直接分块求; 

预处理 $O(n^{\frac{2}{3}})$ , 总复杂度变成 $O(n^{\frac{2}{3}})$ ; 
目前没有解决的是复杂度的分析证明...

# [Project Euler #310](https://projecteuler.net/problem=310)
## 题目大意
只能取完全平方数的nim游戏, 只有三堆, 问必败局面的数量;
## 解题报告
直接暴力sg函数+暴力dp = $5s$, 跑出了答案; 

但是发现sg函数的值非常密集, 最大的不超过在$\sqrt{n}$  所以只需要枚举三堆的sg值, 就可以$n^{\frac{3}{2}}$地得到答案, 这个做法秒出答案; 

## 代码
```cpp
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
typedef vector<int> vi;
typedef pair<int, int> pii;

const int N = 100010;
int f[N];
ll g[2][4][130000];
vi square;
bitset<130000> ap;
int main() {
    f[0] = 0;
    for (int i = 1; i * i <= 100000; ++i)
        square.pb(i * i) ;
    FORU(i, 1, 100000) {
        ap.reset();
        REP(j, sz(square)) {
            if (square[j] > i) break;
            ap.set(f[i - square[j]]);
        }
        REP(j, i + 1) if (!ap[j]) {
            f[i] = j; break;
        }
    }
    mmst(g, -1);
    g[1][0][0] = 1;
    int pre, now;
    for (int i = 0; i <= 100000; ++i) {
        now = i & 1, pre = now ^ 1;
        for (int j = 0; j <= 3; ++j)
            for (int k = 0; k < 130000; ++k)
                for (int l = 0, tmp = 0; l <= j; ++l, tmp ^= f[i])
                    if (g[pre][j - l][k ^ tmp] != -1) {
                        if (g[now][j][k] == -1) g[now][j][k] = 0;
                        g[now][j][k] += g[pre][j-l][k ^ tmp];
                    }
        mmst(g[pre], -1);
    }
    cout << g[0][3][0] << endl;
    return 0;
}
```



# [Project Euler #265](https://projecteuler.net/problem=265)
## 题目大意
将所有长度为$n$的二进制数放在一个环中, 使得顺时钟读取的$2^n$个数不重不漏, 输出$n = 5$时, 所有合法环对应十进制数值的和;
## 解题报告
不小心开启了刷水模式, 但我也不知道这个题这么简单; 

直接枚举环上的二进制位肯定会T, 但是只要预处理出每个二进制数的后继数字(在后4位的末尾加一个0/1), 就可以有效的剪枝; 

另一个比较优越的优化, 就是不需要深搜$2^5$层, 只需要$2^5-4$层就可以, 这是因为为了避免环状置换, 可以强制让第一个数为$0$; 

## 代码
```cpp
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
typedef vector<int> vi;
typedef pair<int, int> pii;

const int N = 34;

vi nxt[N], have;
ll ans = 0;
bool vis[N];
inline int last(int x, int w) {
    int tmp = 0;
    REP(i, w) tmp |= (1 << i) & x;
    return tmp;
}

inline ll got(vi & v) {
    ll x = 0;
    REP(i, sz(v)) x = x << 1 | v[i];
    return x;
}
inline bool legal(int x) {
    REP(i, 4) {
        x = last(x, 4) << 1;
        if (vis[x]) return 0;
    }
    return 1;
}
void dfs(int pas, int x) {
    cout << pas << ' ' << x << endl;
    have.pb(x & 1);
    vis[x] = 1;
    if (pas == 27 && legal(x)) {
        LOOK(have, 0, 31);
        ans += got(have);
    } else
        REP(i, sz(nxt[x]))
            if (!vis[nxt[x][i]])
                dfs(pas + 1, nxt[x][i]);
    have.pop_back();
    vis[x] = 0;
}
int main() {
    REP(i, 32) {
        int tmp = last(i, 4) << 1;
        if (tmp ^ i) nxt[i].pb(tmp);
        tmp = last(i, 4) << 1 | 1;
        if (tmp ^ i) nxt[i].pb(tmp);
    }
    REP(i, 4) have.pb(0);
    dfs(0, 0);
    cout << ans <<  endl;
    return 0;
}

```

# [Project Euler #565](https://projecteuler.net/problem=565)
## Date: 2016/12/25
## 题目大意
$[1, 10^{11}]$ 中， 约数和能被$2017$整除的数的和；
## 解题报告
$$\sum(n) = \prod \sum p^i = \prod \frac{1-p^{k+1}}{1-p}$$

通过打表可以发现，对于一个质数$p$， 只有一个幂指数$k$满足$2017| \frac{1-p^{k+1}}{1-p}$;

并且，通过打表，可以发现，最小的三个$p^k$ 的乘积大于$10^{11}$, 且$10^{11}$内， 满足条件的$p^k$数量很少， 只有$10^4$左右；

但是问题来了，线性筛法， 不能做$10^{11}$内的质数；

继续打表， 可以发现， 对于大于$10^6$的所有合法的$p^k$, $k = 1$， 那么可以枚举$10^{11}$内， $2017$的倍数$y$， 然后$\sqrt{n}$检验$y-1$是否是质数；

当得到所有的合法的$p^k$, 他的所有倍数$x(gcd(x, p^k) = p^k)$添加到答案中，再$O(n^2)$枚举容斥， 就可以AC;

但问题是， 搞不懂为什么自己的程序在windows10上跑到一半就被当做**病毒**删了！

F××K


# [Project Euler #169](https://www.hackerrank.com/contests/projecteuler/challenges/euler169)
## 题目大意
一个数表示成$2^i$形式的方案数， 每个$2^i$可以使用两次；
## 解题报告
通过打表， 可以发现， 设方案数为$f[x]$, 则$f[x] = f[y] (x = 2*y-1)$, $f[x] = f[y] + f[x-1](x = 2 * y)$;

简单证明:
1. 对于$x = 2*y+1$, 可以将$y$的方案中， 所有数全部乘二， 然后加一个单独的1，得到合法方案， 容易证明， 任何将方案中的每个$2^i$拆解， 都会导致与另一方案重复；
2. 对于$x = 2 * y$, 首先可以将$y$的方案中， 所有数乘二，得到若干合法的方案。 然后考虑合法的不与其他方案重复的拆解， 发现只能就某个2拆成2个1，是不会冲突的； 那么去掉拆出的一个1， 就可以等价于f[x-1];

程序实现中，定义了一个30位大整数， 利用两个`long long`实现；

## UPDATE
mrazer 大爷 轻松的搞了一个DP的做法, 比我的科学多了, 复杂度更优越, 关键是, 人家并不是打表看的规律;
链接: http://blog.csdn.net/mrazer/article/details/53898297

## 代码

```cpp
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

const int N = 50;
const ll ba = 1000000000000000LL;
/*
int f[12][1000000];
int main() {
	f[0][0] = 1;
	FORU(i, 0, 10) {
		FORU(j, 0, 1024) f[i+1][j] = f[i][j];
		FORU(j, (1<<i), 1024)
			f[i+1][j] += f[i][j-(1<<i)];
		FORU(j, (2 * (1<<i)), 1024)
			f[i+1][j] += f[i][j-2*(1<<i)];
	}
	FORU(i, 1, 50)
		if ( !(i & 1) )
		cout << i << ": " << f[11][i] << endl;
	return 0;
}
*/
struct bigInteger {
	ll a, b;
	bigInteger() {a = b = 0;}
	bigInteger(ll a, ll b) :a(a), b(b){}
	bool operator<(const bigInteger c) const{
		return a == c.a ? b < c.b: a < c.a;
	}
	bigInteger div2() {
		bigInteger c;
		c.b = b / 2;
		if (a & 1) c.b += ba/2;
		c.a = a / 2;
		return c;
	}
	bigInteger operator - (int x) {
		bigInteger c = *this;
		c.b -= x;
		if (c.b < 0) c.a -= 1, c.b += ba;  
		return c;
	}

};// a*10e15 + b
map<bigInteger, ll> f;
string s;
ll dfs(bigInteger x) {
	if (f[x]) return f[x];
	if (x.b&1) f[x] = dfs(x.div2());
	else f[x] = dfs(x.div2()) + dfs(x-1);
	return f[x];
}
int main() {
	cin >> s;
	FORU(i, 0, (s.length()-1)/2)
		swap(s[i], s[ s.length() - i - 1]);
	ll a = 0, b = 0;
	FORD(i, 14, 0) {
		if (i + 1 > s.length()) continue;
		b = b * 10 + s[i] - '0';
	}
	FORD(i, 28, 15) {
		if (i + 1 > s.length()) continue;
		a =  a * 10 + s[i] - '0';
	}
	f[bigInteger()] = 1; 		
	cout << dfs(bigInteger(a, b)) << endl;
	return 0;
}
```


# [Project Euler #168](https://www.hackerrank.com/contests/projecteuler/challenges/euler168)
## 题目大意
定义一个数$x = 10  * a + b$, 一共有$n$位，$y=b * 10^{n-1}+a$, 如果$\displaystyle y=kx(k \in \mathbb{N})$, 则$x$是合法的， 求$(10, 10^m)$中，合法数的和的后五位；
## 解题报告

比较简单的一道题， 很容易发现， 只要固定个位数字和$k$， 利用$y$ 和$x$的$k$倍关系以及错位关系， 可以得到长度为$l (l \in [2, m])$的有可能合法的数；

这些数合法的条件是最高位和个位在处理完进位后满足$k$的关系；

所以可以通过$O(9 * 9 * m)$得到所有的合法数， 求和即可；

## 代码

```cpp
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
typedef double ff;
typedef vector<int> vi;
typedef pair<int, int> pii;

const int N = 110;

int arr[N], m, ans[N*N], answ;
inline void add_ans(int *a, int w) {
	FORU(i, 1, w) {
		ans[i] += a[i];
		if (ans[i] >= 10) ans[i] -= 10, ans[i + 1] += 1;
	}
	while (ans[w + 1] >= 10)
		++ w, ans[w] -= 10, ans[w + 1] += 1;  
}
int main() {
	cin >> m; mmst(ans, 0);
	FORU(k, 1, 9) FORU(i, 1, 9) {
		mmst(arr, 0); arr[1] = i;
		FORU(j,	 2, m) {
			arr[j] += arr[j - 1] * k;
			arr[j+1] += arr[j] / 10;
			arr[j] %= 10;
			if (arr[j + 1] + arr[j] * k == arr[1] && arr[j] != 0)
				add_ans(arr, j);
		}
	}
	answ = 5;
	while (ans[answ] == 0) answ --;
	if (!answ) ++ answ;
	FORD(i, answ, 1) cout << ans[i]; cout << endl;
	return 0;
} 	
```


# [【PE 153】Investigating Gaussian Integers](https://projecteuler.net/problem=153)
## 题目大意
$\frac{n}{a+bi}$ 若能表示成 $k+bi(k, b \in \mathbb{Z})$ 的形式, 则对于 $n$ , $(a, b)$ 是一个合法的整数对, $s_n = \sum a[(a, b) \ is \ legal \ for \ n]$ , 求 $\displaystyle \sum_{i=1}^{100000000} s_i$ ;
## 解题报告
合法的等价条件:

$$a^2 + b^2 | gcd(a, b)n$$

如果只考虑 $gcd(a, b) = 1$ 的情况, 则 $a^2 + b^2 | n$ ;

如果 $a' = ak, b' = bk$ , 则 $k(a^2 + b^2)|n$ ;

也就是说, 对于 $n$, 如果互质数对 $(a,b), a^2 + b^2 |n$ 那么 $(a, b)$ 可以对 $n$ , 造成贡献 $\displaystyle \sum_{k \times (a^2 + b^2)|n}ka$

所以, 对于 $n$ ,

$$s_n = \sum_{d|n}A(d) \times B(n/d)$$

其中, $\displaystyle A(n) = \sum_{a^2 + b^2 = n, (a, b) = 1} a$ , $\displaystyle B(n) = \sum_{d|n} d$ ;

求 $B$ ? 直接调和筛法, $O(NlnN)$ ;   
求 $A$ ? 把 $a^2 <= 10^8$ 的完全平方数 $n^2$ 组合, 判断互质, 统计, $O(NlnN)$ ;   
求 $s = A \times B$ ? 调和筛, $O(NlnN)$ ;

然后 $O(N)$ , 求 $s$ 的和.


# [【PE 184】Triangles containing the origin](https://projecteuler.net/problem=184)
## 题目大意
在半径为 $r$ 的圆中, 包含原点的整点三角形个数;
## 解题报告
考虑一个三角形 $ABC$ 包含原点的条件, 显然是, $\overrightarrow{OA}, \  \overrightarrow{OB}, \  \overrightarrow{OC}$ 同号;

其实就是固定原点, 按照某个顺序扫过 $ABC$ , 相邻两边转过的角度不超过 $\pi$ ; 
这样就有一个 $O(n^2logn)$ 的做法, 也就是按照极角排序后, 枚举两个点, 二分第三个点的范围; 

直接不进行排序+二分, 考虑合法三角形点的分布情况:
1. 三个点分别位于一个象限, 这样会空出一个象限: 利用圆和坐标轴的对称性, 可以固定空出的象限, 得到答案 $* 4$, 具体的, 三个有点的象限中, 相对的两个象限的点形成的线段与第三点在原点的两侧, 可以 $O(n^2)$ 枚举相对两点 + 叉积判断;
2. 三个点有两个点在一个象限, 另一个点在相对的象限: 枚举单独在一个象限的点, 它与原点的连线将相对的象限分为两部分, 分别统计两部分的点的数量并相乘就可以得到枚举的点对答案的贡献, 同样利用对称性, 答案需要乘 $4$ , 复杂度 $O(n^2)$ ; 
3. 三个点有一个在坐标轴上, 另外两个点在对侧相邻的两个象限: 这个情况最简单, 答案就是坐标轴本轴点的数量乘象限内点数量的平方, 答案需要乘 $4$ ;
4. 三个点有一个在坐标轴上, 另外两个在相对的两个象限内: 使用坐标半轴上的点的数量乘第一种情况中相对两点合法个数, 答案乘 $8$ ; 

5. 两个点在不同的坐标轴上, 第三点在相对的象限内: 用坐标半轴的点数乘一个象限内的点数, 答案乘 $4$ ; 

## 代码

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
const int N = 10100;
int x;
ll s;
vector< pii > t;
ll ans;
int main() {
    int n = 105;
    x = n - 1;
    FORU(i, 1, n - 1) {
        int haf = sqrt(n * n - i * i);
        if (haf * haf + i * i == n * n)
            -- haf;
        FORU(j, 1, haf) t.pb(mp(i, j));
    }
    int _x, _y, __x, __y;
    REP(i, sz(t)) REP(j, sz(t)) {
        _x = t[i].fi, _y = t[i].se;
        __x = -t[j].fi, __y = -t[j].se;
        if (_x * __y - __x * _y > 0)
            ++ s; //cout << _x << ' ' << _y << ' ' << __x << ' ' << __y << endl;
    }
    // cout << s << endl;
    REP(i, sz(t)) {
        int a = 0, b = 0;
        REP(j, sz(t)) {
            _x = t[i].fi, _y = t[i].se;
            __x = -t[j].fi, __y = -t[j].se;
            if (_x * __y - __x * _y > 0)
                ++ a;
            else if (_x * __y - __x * _y < 0)
                ++ b;
        }
//        cout << _x << ' ' << _y << " : " << a << ' ' <<b << endl;
        ans += a * b;
    }
    ans = ans * 4;
    ans += (s + x*x) * sz(t) * 4;
    ans += 4LL * x * sz(t) * sz(t);
    ans += 8 * x * s;
    cout << ans << endl;
    return 0;
}
```


# [【PE 156】Counting Digits](https://projecteuler.net/problem=156)
## 题目大意
$f[n][d]$ 表示前 $n$ 个数中， 数字 $d$ 出现的次数， 求 $\displaystyle \sum_n \sum_{d = 1}^{9} [f[n][d] = n]$
## 解题报告
利用数位dp, 可以得到第 $f[n][d]$ , 时间复杂度是 $O(logn)$ ;   
但是怎么求出所有满足条件的 $(n, d)$ ?  
对于不同的 $d$ , 寻找满足条件的 $n$ , 显然, 如果 $n = n + 1$ , $f[n][d]$ 可能不变, 最多增加 $w$ , $w$ 是 $n$ 的位数;   
那么如果 $abs(n - f[n][d]) = delta$ , 下一个满足条件的 $n$ 最少需要经过 $delta/w$ , 不妨令 $w = 10$ , 因为通过尝试, 发现最大的满足条件的 $n$ 在 `int` 范围内...  
复杂度? 不是很会证, 但可以在 $1s$ 内得到答案;

## 代码

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
typedef unsigned long long ull; 
const int N = 20;
ull ans, lim = 1LL << 50, step;
ull pre[N], uu[N];
inline ull got(ull n, ull d) {
    vi w;
    while (n) w.pb(n % 10), n /= 10;
    ull ca = 0, ans = 0; int x;
    FORD(i, sz(w) - 1, 0) {
        x = w[i];
        ans += ca * uu[i] * x;
        ans += pre[i] * x;
        if (x > d) ans += uu[i];
        if (x == d) ca ++;
    }
    ans += ca;
    return ans;
}

inline ull query(int d) {
    ull n = 1,  sm = 0; ll step = 1;
    while (n < lim) {
        ull tmp = got(n, d);
        if (tmp - n == 0) sm += n;// cout << d << ": " << n << endl;
        if (tmp > n) step = (tmp - n) / 10;
        else step = (n - tmp) / 10;
        if (!step) step = 1;
        n += step;
    }
    if (d == 1) cout << sm << endl;
    return sm;
}

int main() {
    ull tot = 0,  tmp = 1 ;
    uu[0] = 1;
    FORU(i, 1, 16) {
        tmp = tmp * 10;
        pre[i] = tmp * i / 10 ;
        uu[i] = tmp;
    }
    FORU(i, 1, 9) ans += query(i);
    cout << ans << endl;
    return 0;
}
```

# [【PE 157】Solving the diophantine equation](https://projecteuler.net/problem=157)
## 题目大意

如果 $\frac{1}{a} + \frac{1}{b} = \frac{p}{10^n}$ , 其中 $a, b, p \in \mathbb{Z}, n \in [1, 9]$ , 那么 $(a, b, p, n)$ 向答案贡献 $1$ , 求答案。

## 解题报告
对题目中的条件进行翻译:
$$\frac{1}{a} + \frac{1}{b} = \frac{p}{10^n} \Leftrightarrow ab|10^n(a+b)$$
如果 $(a, b) = 1$ , 那么 $ab|10^n$ ;  
如果 $(a, b) = k$ , 那么 $kab|10^n$ ;   
所以, 可以枚举 $(a, b) = 1$ , 使得 $ab|10^n$ , 贡献是 $\frac{10^n(a+b)}{ab}$ 的约数;

可以发现, 能乘除 $10^n$ 的 $a, b((a, b) = 1)$ 导致两种情况: ①$a = 1, b|10^n$ ② $a=2^l, b=5^k, ab|10^n$ , 第二种情况最多 $log_210^n$ 种, 第一种情况大约 $ln10^n$ 种, 所以可以枚举;

## 代码

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
typedef vector<int> vi;
typedef pair<int, int> pii;

int lim = 1e9, ans;
vi xy, x, y;
vi w[1000000];
inline int query(int n, int i) {
    ll _n = 1LL * n * (x[i] + y[i]) / xy[i];
    int tmp = 1;
    for (ll j = 2; j * j <= _n; ++j)
        if (_n % j == 0) {
            int k = 0;
            while (_n % j == 0) ++k, _n/=j;
            tmp = tmp * (k + 1);
        }
    if (_n > 1) tmp = tmp * 2;
    if (x[i] != 1 && y[i] != 1) {
        int _tmp = 1;
        _n = 1LL * n * (1 + x[i]*y[i]) / xy[i];
        for (ll j = 2; j * j <= _n; ++j)
            if (_n % j == 0) {
                int k = 0;
                while (_n % j == 0) ++k, _n/=j;
                _tmp = _tmp * (k + 1);
            }
        if (_n > 1) _tmp = _tmp * 2;
        tmp += _tmp;
    }
    if (xy[i] < 1000000) w[xy[i]].pb(tmp);
    return tmp;
}

int main() {
    ans = 0;
    for (int i = 1; i < lim; i *= 2)
        for (int j = 1; j < lim; j *= 5) {
            if (1LL * i * j > lim) break;
            xy.pb(i * j), x.pb(i),  y.pb(j);
        }
    for (int i = 10; i <= lim; i *= 10) {
        int tmp = 0;
        for (int j = 0; j < sz(xy); ++j)
            if (i % xy[j] == 0)
                tmp += query(i, j);
        ans += tmp;

    }
    cout << ans << endl;
    return 0;
}
```

# [【PE 160】Factorial trailing digits](https://projecteuler.net/problem=160)
## 题目大意
$10^{12}!$ 去零的后五位;
## 解题报告
1. 首先确定 $10^{12}!$ 次方有几个 $0$ , 其实就是 $10^{12}!$ 中 $10$ 的幂指数 $k_{10}$ , 而 $k_{10} = min(k_{2}, k_{5}) , k_{2} > k_{5}$ , 所以 $k_{10} = k_{5}$ ;  
计算 $k_{5}$ ? 令 $n = 10^{12}$ , $5$ 的倍数会分别贡献一个 $5$ , $25$ 的倍数会额外贡献一个 $5$ ...   
所以 $\displaystyle k_{5} = \sum_{i} \frac{n}{5^i}$   
记 $k = k_5$ ;
2. 确定 $\frac{n!}{10^k}(mod \ 10^5)$  
这个比较有难度, 因为 $n!$ 中, 有 $5$ , $2$ , 无法直接取模, 但是因为 $n$ 足够的大,所以可以得到, $k_2 - k_5 > 5$ , 也就是 $\frac{n!}{10^k}(mod \ 2^5) = 0$ ,所以如果能求出 $\frac{n!}{10^k}(mod \ 5^5)$ , 使用CRT就可以得到答案(因为 $(2^5, 5^5) = 1$ );  
将 $n!(mod \ 5^5)$ 中的 $5$ 进行递归提取, 则 $n_i! = 5^{k_5(n_i)} A(n_i) n_{i+1}!$ , 其中 $n_{i+1} = \lfloor \frac{n_i}{5} \rfloor$ .
展开递归式, $n! = 5^{k}A(n!)A(n!/5)A(n!/25)...A(1)$ , 其中 $A(n)$ 表示 $[1, n]$ 中与 $5$ 互质的数的乘积, 也就是 $A(n) = \prod_{i=1}^{n} [gcd(i, 5)=1]i$ ;  
这个 $A(n)$ 显然存在循环节, 且 $5^5$ 是一个已知的循环节 , 那么可以预处理 $A(1-5^5)$ ;  
$$\frac{5!}{10^k}(mod \ 5^5) = \frac{5^{k}A(n!)A(n!/5)A(n!/25)...A(1)}{5^k2^k}(mod \ 5^5) = 2^{-k}A(n!)A(n!/5)A(n!/25)...A(1)(mod \ 5^5)$$
3. 使用中国剩余定理 .

过程中的逆元需要用拓展欧几里得求.
## 代码

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
const int N = 3126;
int A[N], pre[N] ;
ll n = 1000000000000LL, d;
inline ll quick(ll x, ll k) {
     ll ans = 1;
     for (; k; k >>= 1) {
         if (k & 1) ans = ans * x % 3125;
         x = x * x % 3125;
     }
     return ans;
}
inline void exd_ou(ll a, ll b, ll &x, ll &y) {
    if (!b) x = 1, y = 0;
    else {
        exd_ou(b, a % b, y , x);
        y -= (a/b) * x;
    }
}
inline ll _A(ll n) {
    return ((int)pow(-1, n/3125) * A[n % 3125] + 3125) % 3125;
}

int main() {
    ll fiv = 5;

    for (; fiv <= n; fiv *= 5)
        d += n / fiv;
    A[0] = 1;
    for (int i = 1; i < N; ++i)
        if (i % 5) A[i] = A[i-1] * i % 3125;
        else A[i] = A[i-1];
    ll two, rev, tmp ;
    two = quick(2LL,  d);
    exd_ou(two, 3125, rev, tmp);
    if (rev < 0) rev += 3125;
    fiv = 1; ll _a = 1;
    for (; fiv <= n; fiv *= 5)
        _a = _a * _A(n / fiv) % 3125;
    _a = rev * _a % 3125;
    two = quick(2, 5);
    exd_ou(two, 3125, rev, tmp);
    ll ans = _a * two * rev;
    cout << ans % 100000<< endl;
    return 0;
}
```


