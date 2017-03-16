---
title: 【bzoj 2599】[IOI2011]Race
date: 2017-02-14 14:14:28
tags:
  - 树链剖分
  - 树分治
categories:
  - 数据结构
---

> 聪哥和TA爷在做的一道题, 之前好像所有人都当做点分治模板题做的, TA爷用dsu on a tree艹了, 我也更欣赏这个常数优越的做法;

<!--more-->

# 题目大意
给一棵树,每条边有权.求一条简单路径,权值和等于$K$,且边的数量最小.$N <= 200000$, $K <= 1000000$

# 解题报告
求满足条件的树的路径, 首先想到的是点分治, 在每个分治块, 从分治重心出发, 可以得到每个点到重心的距离, 并更新答案;

在点分治的过程中可以剪枝, 有两点:①已经经过的距离大于$k$, 停止`dfs`,  ②已经经过的边数大于`ans`, 停止`dfs`;

因为点分治本身不优越的常数, 所以这个做法比较慢(加了剪枝后还比较理想);

给出几个同学使用点分治的解题报告, [dada's](http://www.cnblogs.com/DaD3zZ-Beyonder/p/5618633.html) [abclzr's](http://www.cnblogs.com/abclzr/p/5337088.html)

再次扔出[yveh的dsu on the tree](http://blog.csdn.net/qaq__qaq/article/details/53455462), 是这个题常数比较小的一个做法; 

极其模板, 路径的长度可以表示为$dis[x] + dis[y] - dis[lca] * 2$, 所以可以利用记录 $dis$ 来在 $lca$ 处得更新 `ans`; 

每次先将轻儿子子树中的事情处理好, 清空数组, 然后, 处理重儿子的信息, 不清空数组, 将当前点和轻儿子子树的点重新加入数组得到答案; 

但是 $dis[x]$ 可以很大, 数组并不能很好存下, 那么我们可以利用模意义, 在 $mod \ k$ 的意义下记录信息, 同时记录原数大小, 更新答案的时候, 需要判断原数是否满足题意; 

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

const int N = 200100;
const int inf = N;

ll dis[N];
int n, k, dep[N], sz[N], son[N], ans = inf;
int sho[1000000], id[N], cnt, ent[N], out[N];
vi sta, to[N], v[N];

inline void in(int &x) {
    char ch = getchar(); int f = 1;
    for (; ch < '0' || ch > '9'; ch = getchar())
        if (ch == '-') f = -1;
    for (x = 0; ch >= '0' && ch <= '9'; ch = getchar())
        x = x * 10 + ch - 48;
    x = x * f;
}

void dfs(int x, int from) {
    sz[x] = 1, son[x] = n, id[++cnt] = x;
    ent[x] = cnt;
    REP(i, sz(to[x]))
        if (to[x][i] != from) {
            dis[to[x][i]] = dis[x] + v[x][i];
            dep[to[x][i]] = dep[x] + 1;
            dfs(to[x][i], x);
            sz[x] += sz[to[x][i]];
            if (sz[to[x][i]] > sz[son[x]])
                son[x] = to[x][i];
        }
    out[x] = cnt;
}

inline void getans(ll &di, int &de, ll &_di, int &_de) {
    int __di = (_di * 2 % k - di % k) + k;
    if (__di >= k)  __di -= k;
    if (sho[ __di ] != -1 && dis[sho[__di]] + di -_di * 2 == k)
        ans = min(ans, de + dep[sho[__di]] - _de * 2);
}
inline void enter(int & x) {
    int _di = dis[x] % k;
    if (sho[_di] == -1) sho[_di] = x;
    else if (dis[x] < dis[sho[_di]]) sho[_di] = x;
        else if (dis[x] == dis[sho[_di]] && dep[x] < dep[sho[_di]])
            sho[_di] = x;
    sta.pb(_di);
}

void smart(int x, int f, bool top) {
    if (son[x] !=  n) {
        REP(i, sz(to[x])) if (to[x][i] ^ f)
            if (to[x][i] ^ son[x])
                smart(to[x][i], x, 1);
        smart(son[x], x, 0);
    }
    getans(dis[x], dep[x], dis[x], dep[x]);
    enter(x);
    REP(i, sz(to[x])) if (to[x][i] ^ f)
        if (to[x][i] ^ son[x]) {
            int y = to[x][i];
            FORU(i, ent[y], out[y])
                getans(dis[id[i]], dep[id[i]], dis[x], dep[x]);
            FORU(i, ent[y], out[y])
                enter(id[i]);
        }
    if (top) {
        REP(i, sz(sta)) sho[sta[i]] = -1;
        sta.clear();
    }
}

int main() {
    in(n), in(k);
    int x, y, z;
    REP(i, n - 1) {
        in(x), in(y), in(z);
        to[x].pb(y), to[y].pb(x);
        v[x].pb(z), v[y].pb(z);
    }
    dep[0] = 1, dfs(0, -1);
    mmst(sho, -1);
    smart(0, -1, 0);
    printf("%d\n", ans == N ? -1 : ans);
    return 0;
}
```

$MLE$ 了两发, 发现自己的数组范围开大了10倍, 太惨了; 

# 原题信息

[题目链接](http://www.lydsy.com/JudgeOnline/problem.php?id=2599)

<body>
<title>Problem 2599. -- [IOI2011]Race</title><center><h2>2599: [IOI2011]Race</h2><span class="green">Time Limit: </span>70 Sec&nbsp;&nbsp;<span class="green">Memory Limit: </span>128 MB<br><span class="green">Submit: </span>3006&nbsp;&nbsp;<span class="green">Solved: </span>877<br>[<a href="submitpage.php?id=2599">Submit</a>][<a href="problemstatus.php?id=2599">Status</a>][<a href="bbs.php?id=2599">Discuss</a>]</center><h2>Description</h2><div class="content"><p>给一棵树,每条边有权.求一条简单路径,权值和等于K,且边的数量最小.N &lt;= 200000, K &lt;= 1000000</p>
<p></p></div><h2>Input</h2><div class="content"><p>第一行 两个整数 n, k<br>
第二..n行 每行三个整数 表示一条无向边的两端和权值 (注意点的编号从0开始)</p></div><h2>Output</h2><div class="content"><p>一个整数 表示最小边数量 如果不存在这样的路径 输出-1</p></div><h2>Sample Input</h2>
			<div class="content"><span class="sampledata">4 3<br>
0 1 1<br>
1 2 2<br>
1 3 4<br>
</span></div><h2>Sample Output</h2>
			<div class="content"><span class="sampledata">2<br>
</span></div>

</body>