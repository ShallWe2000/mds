---
title: codeforces-选做
date: 2017-03-14 07:47:52
tags:
  - codeforces
  - 线段树
  - STL
  - FFT
  - 倍增
  - 平衡树
  - 贪心
  - 构造
  - DP
categories:
  - 题目集锦
---

> 选做一些codeforces 上的~~无脑码农题~~。

<!--more-->

# 【CF 739C】Alyona and towers
## [题目链接](http://codeforces.com/contest/739/problem/C)
## 题目大意

两个操作：
1. 区间加法；
2. 查询最长上凸单峰子串长； 

## 解题报告

$$a_l < a_{l + 1} < a_{l + 2} < ... < a_k > a_{k + 1} > a_{k + 2} > ... > a_r$$

差分之后， 这个式子等价于一串连续的正数拼一串连续的负数； 

可以直接用**线段树**维护三个信息： 左侧开始最长合法串， 右侧开始最长合法串， 不限制位置的最长合法串；

合并的时候讨论节点处两点的正负； 

区间加，差分后变成单点加； 

## 代码

```c++
#include <bits/stdc++.h> 
using namespace std; 
#define foru(i, s, t) for (i = (int)s; i <= (int) t; ++i)   
#define ford(i, s, t) for (i = (int)s; i ^ (int) t; --i) 
#define pb push_back
#define mp make_pair
#define F first
#define S second
#define check cout << 1 << endl
#define dou (", ") 
#define mao (": ")
typedef pair <int, int> pii; 
typedef vector <int> vi; 
typedef long long ll; 
typedef double ff; 
const int N = 310000; 
const int inf = 1000000000; 
int lan[N << 2], ran[N << 2], ans[N << 2], n, m; ll a[N];
template<typename T> 
inline void in(T &x ) { 
    char ch = getchar(); int f= 1; 
    for (; ch < '0' || ch > '9'; ch = getchar())    
        if ( ch == '-') f =  -1; 
    for (x = 0; ch >= '0'  && ch <= '9'; ch = getchar()) 
        x = x * 10 + ch - 48; 
    x *= f; 
}
inline void newnode(int x, ll va) {
    lan[x] = ran[x] = ans[x] = (va != 0); 
}
inline void push_up(int x, int l, int r) { 
    int mid = l + r >> 1; 
    ans[x] = max(ans[x << 1], ans[x << 1 | 1]); 
    lan[x] = lan[x << 1] , ran[x] = ran[x << 1 | 1]; 
    if (!(a[mid] < 0 && a[mid + 1] > 0) && a[mid] && a[mid + 1])  {
        ans[x] = max(ans[x], ran[x << 1] + lan[x << 1 | 1]);     
        if (lan[x] == (mid - l + 1) )  lan[x] += lan[x << 1 | 1]; 
        if (ran[x] == (r - mid)) ran[x] += ran[x << 1]; 
    }
}
    
void build(int x, int l, int r) { 
    if (l > r) return; 
    if (l == r) newnode(x, a[l]); 
    else { 
        int mid = l + r >> 1;
        build(x << 1, l, mid); 
        build(x << 1 | 1, mid + 1, r) ; 
        push_up(x, l, r); 
    }
}
void add(int x, int l, int r, int go, int delta) {  
    if (l > r) return; 
    if (l == r) a[l] += delta, newnode(x, a[l]); 
    else { 
        int mid = (l + r) >> 1; 
        if (go <= mid) add(x << 1, l, mid, go, delta); 
        if (go > mid)  add(x << 1 | 1, mid + 1, r, go, delta); 
        push_up(x, l, r); 
    }
}
int main() {
    in(n); int i, j, x, y, d; 
    foru(i, 1, n) in(a[i]);    
    foru(i, 1, n) a[i] = a[i+1] - a[i]; 
    --n, build(1, 1, n), in(m);
    foru(i, 1, m) { 
        in(x), in(y), in(d); 
        if (x > 1) add(1, 1, n, x-1, d);
        if (y <= n) add(1, 1, n, y, -d);   
        printf("%d\n", ans[1] + 1); 
    }
    return 0; 
}   
```

# 【CF 739E】Gosha is hunting
## [题目链接](http://codeforces.com/contest/739/problem/E)
## 题目大意

$n$个小精灵，$a$个宝贝球，$b$的大师球，每个小精灵可以用宝贝球或大师球或两个球捕捉，成功概率为$p_i$,$u_i$,$p_i + u_i - p_i * u_i$ 求期望最大捕捉数; （摘自[zky学长的blog](http://kzoacn.xyz/index.php/archives/238))

## 解题报告

首先可以写一个$O(n^3)$的DP,  令$f[n][a][b]$表示前$n$个小精灵，用了$a$个宝贝球，$b$个大师球， 期望最大的捕捉数， 容易得到：
$$f[n][a][b]=max \{ f[n-1][a-1][b]+p_n, f[n-1][a][b-1]+u_n,f[n-1][a-1][b-1]+p_n+u_n-p_n*u_n \}$$

这个DP是可以AC的。。。

从zky学长那里学了一个很优秀的做法，和去年IOI的最后一题做法相似：

对于一个可以求最值的凸函数$f(x)$,  构造$F(x)=f(x)-kx$, 那么$F(x)$为凸函数且可 求最值， 通过二分$k$使得$F(x)$的最值点左右平移， 当最值点位于$x$， 可以得到$f(x)=F(x)+kx$, 从而实现求$f(x)$每个自变量对应的函数值。 

对于这个题， 因为某些不可说的原因，$f[i][j](k)$是一个凸函数， 利用上述方法可以通过$O(\log^{-1}eps \times n^2)$得到答案。

## 代码

```c++
#include <bits/stdc++.h> 
using namespace std; 
#define foru(i, s, t) for (i = (int)s; i <= (int) t; ++i)   
#define ford(i, s, t) for (i = (int)s; i ^ (int) t; --i) 
#define pb push_back
#define mp make_pair
#define F first
#define S second
typedef pair <int, int> pii; 
typedef vector <int> vi; 
typedef long long ll; 
typedef double ff; 
const int N = 2010; 
const int inf = 1000000000; 
const ff eps = 1e-9;
int n, a, b; 
ff p[N], u[N], mix[N], f[N][N], cnt[N][N]; 
int main() {
    scanf("%d%d%d", &n, &a, &b); int i, j; 
    foru(i, 1, n) scanf("%lf", &p[i]); 
    foru(i, 1, n) scanf("%lf", &u[i]), 
    mix[i] = p[i] + u[i] - p[i] * u[i]; 
    ff l = -1e3, r = 1e3, k; 
    while (l + eps < r) {
        k = (l + r) * 0.5; 
        foru(i, 0, n) foru(j, 0, i) 
            f[i][j] = -1e10, cnt[i][j] = 0;     
        f[0][0] = 0; 
        foru(i, 1, n) foru(j, 0, i) {
            if (j < i && f[i-1][j] > f[i][j]) f[i][j] = f[i-1][j], cnt[i][j] = cnt[i-1][j]; 
            if (j && f[i-1][j-1] + p[i] > f[i][j]) f[i][j] = f[i-1][j-1] + p[i], cnt[i][j] = cnt[i-1][j-1]; 
            if (j < i && f[i-1][j] + u[i] + k> f[i][j]) f[i][j] = f[i-1][j] + u[i] + k, cnt[i][j] = cnt[i-1][j] + 1; 
            if (j && f[i-1][j-1] + mix[i] + k> f[i][j]) f[i][j] = f[i-1][j-1] + mix[i] + k, cnt[i][j] = cnt[i-1][j-1] + 1; 
        }
        if (cnt[n][a] < b) l = k; else r = k; 
    }
    printf("%.10lf\n", f[n][a] - (l + r) * 0.5 * b); 
    return 0; 
}
```

# 【CF 736D】Permutations
## [题目链接](http://codeforces.com/contest/736/problem/D)
## 题目大意
要得到一个长度为$n$的排列，给出$m$个信息$(x,y)$，表示$x$可以放置在$y$这个位置，保证使用这$m$个信息，方案数为奇数， 要得到对每个信息而言，如果不使用，方案数是奇数还是偶数；
## 解题报告

![评论](1.jpg)

首先求二分图的完美匹配数， 实际上是联通矩阵的积和式的值， 在$\pmod 2$的意义下， 与行列式等价。 

所以问题等价于求一个行列式去掉$(i,j)$处的$1$，剩下的行列式的值。

考虑将行列式沿第$i$行展开，那么删掉$1$造成的影响就是$j$位置的代数余子式的值， 而又有伴随矩阵$adj(A)=C^{T}(A)$, 其中$C(A)$代数余子矩阵， 且$adj(A)$可以通过$adj(A)=A^{-1}*det(A)$快速求解， 所以利用bitset压位， 预处理复杂度是$O(n^2\log n)$的。  

## 代码

```c++
#include <bits/stdc++.h> 
using namespace std; 
#define FOR(i, a, b) for (int i = (a), _b = (b); i <= b; ++i) 
#define FORD(i, a, b) for (int i = (a), _b = (b); i >= b; --i) 
#define REP(i, b) for (int i = 0, _b = (b); i < b; ++i) 
#define pb push_back
#define fi first
#define se second
#define mp make_pair
#define SZ(x) ((int)(x).size())

typedef vector<int> vii; 
typedef pair<int, int> pii; 
typedef long long ll; 
typedef double ff; 

const int N = 2010; 

int n, m, l[N * N], r[N * N]; 
bitset<N << 1> a[N]; 
int main() { 
	ios :: sync_with_stdio(false); 
	cin >> n >> m; int x, y; 
	FOR(i, 1, m) {
		cin >> x >> y; 
		a[x].set(y), l[i] = x, r[i] = y; 
	}
	FOR(i, 1, n) a[i].set(i + n); 
	FOR(i, 1, n) { 
		int j = i; 
		while (!a[j].test(i)) ++j; 
		swap(a[i], a[j]); 
		FOR(k, 1, n) if ( k != i && a[k].test(i)) 
			a[k] ^= a[i]; 
	}
	FOR(i, 1, m)
		puts(a[ r[i] ].test(n + l[i]) ? "NO" : "YES"); 
	return 0; 
}
```

#  【CF 741C】Arpa’s overnight party and Mehrdad’s silent entering
## [题目链接](http://codeforces.com/contest/741/problem/C)
## 题目大意
每相邻的三个位置不能拥有相同类型的食物, 情侣不能拥有相同类型的食物;

构造一个合法的食物分配方法或者输出无解;
## 解题报告
显然, 如果没有相邻三个的限制, 这个题一定有解; 

如果让第$k* 2$和第$k* 2+1$个类型不同, 那么相邻的三个就一定不全相同, 根据二分图染色的性质, 结合图的特殊性, 可以得到, 这个图不存在奇环, 所以一定可以构造出二分图; 

那么就可以直接二分图染色; 

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
#define REP(i, b) for (int i = a, nn = int(b); i <= nn; ++i)

typedef pair<int, int> pii;
typedef vector<int> vi; 

const int N = 200100;

vi nxt[N];
vector< pii > par;
int n, ans[N];
void dfs(int x, int a) {
	ans[x] = a;
	FORU(i, 0, sz(nxt[x])-1)
		if (ans[nxt[x][i]] == -1)
			dfs(nxt[x][i], a ^ 1) ;
}
int main() {
 	ios :: sync_with_stdio(false);
	int x, y;
	cin >> n;
	FORU(i, 1, n) {
 		cin >> x >> y;
		nxt[x].pb(y), nxt[y].pb(x);
		par.pb( mp(x, y) );
	}
	FORU(i, 1, n) {
		nxt[i * 2 - 1].pb(i * 2);
		nxt[i * 2].pb(i * 2 - 1);
	}
	mmst(ans, -1);
	FORU(i, 1, n*2)
		if (ans[i] == -1)
			dfs(i, 0);
	FORU(i, 0, sz(par) - 1)
		cout << ans[par[i].fi] + 1 << ' ' << ans[par[i].se] + 1 << endl;
	return 0;
}
```

# 【CF 741D】Arpa’s letter-marked tree and Mehrdad’s Dokhtar-kosh paths
## 题目大意
一棵树, 每个节点有一个字符, 求出每个点的子树中, 最长的可以通过改变顺序构成回文串的链;
## 解题报告
翻译一下, 将每个字符当做一个二进制位上的一, 对链的要求实际上等价于异或和为0(长度为偶的回文串), 或者异或和的某个数位为1(长度为奇的回文串);<br>
一个链 $(u,v)$ , 异或和为$v[u] \ xor \ v[v] \ xor \ a[lca]$, 其中$v[x]$ 表示节点$x$到根节点的异或和, $a[x]$ 表示 $x$ 节点自己的二进制数值;<br>
所以可以用一个数组将每个二进制数所表示的 $v[x]$ 对应的深度最大的 $x$ 的深度记录下来, 在 $lca$ 处通过枚举为一的数位, 用 $O(logn)$ 查询答案;<br>
关键是, 怎么用数组记录子树中的信息?<br>
显然需要使用树上的启发式合并, 具体一点, 就是"dsu on the tree", 简单一点, 就是树链剖分;<br>
这个可以看 [yveh的一个学习笔记](http://blog.csdn.net/qaq__qaq/article/details/53455462)<br>
代码很好写, 可以算是裸题;<br>

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
#define REP(i, b) for (int i = a, nn = int(b); i <= nn; ++i)
#define sz(x) (int)(x).size()
#define mmst(a, x) memset(a, x, sizeof(a))

typedef vector<int> vi;
typedef pair<int, int> pii;

const int N = 500010;

int n, toTwo[200], sons[N], in[N], out[N], id[N];
int heav[N], bag[1 << 23], cnt = 0, pref[N], dep[N], ans[N];
vi to[N], alpha[N], tmp;
inline void up(int &ans, int x) {
	if (x > ans ) ans = x;
}
void dfs(int x) {
	sons[x] = 1, id[++cnt] = x, in[x] = cnt, heav[x]= 0;  
	FORU(i, 0, sz(to[x]) - 1) {
		pref[to[x][i]] = pref[x] ^ alpha[x][i];
		dep[to[x][i]] = dep[x] + 1, dfs(to[x][i]);
		if (sons[ to[x][i] ] > sons[ heav[x] ])
			heav[x] = to[x][i];
		sons[x] += sons[ to[x][i] ];
	}
	out[x] = cnt;
}
void got_ans(int x, int type) {
// type: 1-clear, 0-not
//  tmp: a stack recording vertices in the heavy subtrees;
	ans[x] = 0;
	FORU(i, 0, sz(to[x]) - 1) 	
		if (to[x][i] != heav[x])
			got_ans(to[x][i], 1), up(ans[x], ans[to[x][i]]);
	if (heav[x]) 	
		got_ans(heav[x], 0), up(ans[x], ans[ heav[x] ]);
//  information of subtrees updates ans of x;
	int y, z;
	FORU(i, 0, sz(to[x]) - 1)
		if (to[x][i] != heav[x]) {
			y = to[x][i]; vi _tmp;
			FORU(j, in[y], out[y]) {
				z = id[ j ], tmp.pb(z), _tmp.pb(z);
				if (bag[ pref[z] ] != -1)
					up(ans[x], bag[ pref[z] ] + dep[z] - dep[x] * 2);
				FORU(k, 'a', 'v')
					if (bag[ pref[z] ^ toTwo[k] ] != -1)
						up(ans[x], bag[ pref[z] ^ toTwo[k] ] + dep[z] - dep[x] * 2); 				
//  update ans by forming a line with two
			}
			FORU(j, 0, sz(_tmp) - 1) {
				z = _tmp[j];
				if (dep[z] > bag[ pref[z] ])
					bag[pref[z]] = dep[z];
//  update BAGS
			}				
		}
	if (bag[ pref[x] ] != -1)
		up(ans[x], bag[pref[x]] - dep[x]);
	FORU(k, 'a', 'v')
		if (bag[ pref[x] ^ toTwo[k] ] != -1)
			up(ans[x], bag[ pref[x] ^ toTwo[k] ] - dep[x]);
	if (dep[x] > bag[ pref[x] ]) bag[pref[x]] = dep[x];
	tmp.pb(x);
// use x to update
	if (type == 1) {
 		FORU(i, 0, sz(tmp) - 1) {
			y = tmp[i];
			if (bag[ pref[y] ] != -1)
				bag[ pref[y] ] = -1;
		}
		tmp.clear();
// meet the light path and clear
	}
}

int main() {
	ios :: sync_with_stdio(false);
	cin >> n; int pa; string c;
	FORU(i, 'a', 'v')
		toTwo[i] = 1 << i - 'a';
	FORU(i, 2, n) {
		cin >> pa >> c, to[pa].pb(i);
		alpha[pa].pb(toTwo[c[0]]);
	}
	mmst(bag, -1), dfs(1), got_ans(1, 1);
	FORU(i, 1, n)
		cout << ans[i] << ((i == n)? '\n' : ' ')<< endl;
	return 0;
}
```

# 【CF 750E】New Year and Old Subsequence
## [题目链接]( "http://codeforces.com/contest/750/problem/E")
## 题目大意
给出一个字符串, 询问一个区间$[l, r]$中, 使子序列不存在"2016"但存在"2017", 最少需要删除几个字符; 
## 解题报告
这竟然是一道E题, 而且现场只有100个人A; <br> 
一个非常有效的性质, 考虑不存在区间查询这种事情, 就求一个给出的字符串, 需要删除几个字符;<br>
那这个是一个很简单的动态规划, 令 $f[i][j]$ 表示到第 $i$ 个字符, 匹配到第 $j$ 个字符, 最少删几个, 转移显然;<br> 
发现, 某个位置的字符决定这个位置对转移的影响, 也就是经过一个位置, 存在固定不变的转移;<br> 
用一个$5 * 5$的矩阵可以表示每个位置的转移, 然后查询的时候, 只需要做区间矩阵乘法, 用线段树可以解决; <br> 
实际上, 题目中的矩阵可以理解为维护任意两点间的最短路, 即$f[i][j]$表示匹配$i$个字符到匹配$j$个字符, 最少需要删几个字符; <br> 
用了刚学习的zkw线段树, 可能哪里写得不好, 跑得竟然很慢, 但确实很短;<br> 
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
#define sz(x) (int)(x).size()
#define mmst(a, x) memset(a, x, sizeof(a))

typedef vector<int> vi;
typedef pair<int, int> pii;

const int N = 200100; 
inline void minify(int &a, int b) { 
    if (b < a) a = b; 
} 

struct mat { 
    int a[5][5]; 
    mat() {mmst(a, 127/3);} 
    mat operator *(mat b) { 
        mat ans; 
        REP(i, 5) REP(j, 5) REP(k, 5) 
            minify(ans.a[i][k], a[i][j] + b.a[j][k]); 
        return ans; 
    } 
    void init(char x) { 
        REP(i, 5) a[i][i] = 0; 
        if (x == '6') a[3][3] = a[4][4] = 1; 
        string s = "2017"; 
        REP(i, 4) if (x == s[i])
            a[i][i + 1] = 0, a[i][i] = 1; 
    }
} ma[N * 4];
int n, q, m; 
string s;
int main() { 
    ios :: sync_with_stdio(false); 
    cin >> n >> q >> s; 
    for (m = 1; m <= n; m <<= 1); 
    FORU(i, 1, n) ma[m + i].init(s[i - 1]); 
    FORD(i, m-1, 1) ma[i] = ma[i << 1] * ma[i << 1 | 1]; 
    REP(i, q) { 
        int a, b; 
        cin >> a >> b; 
        a += m, b += m; 
        mat lef = ma[a], rig = ma[b]; 
        while (a < b - 1) { 
            if (!(a & 1)) lef = lef * ma[a + 1]; 
            if (b & 1) rig = ma[b - 1] * rig; 
            a >>= 1, b >>= 1; 
        } 
        lef = lef * rig; 
        if (lef.a[0][4] > N) cout << -1 << endl; 
        else
            cout << lef.a[0][4] << endl; 
    } 
    return 0; 
} 
```

# 【CF 717F】Heroes of Making Magic III
## [题目链接]( "http://codeforces.com/contest/717/problem/F")
## 题目大意
主人公可以从 $i$ 走到 $i-1$ 和 $i+1$, 给出 $s$ , $s_i$ 表示点 $i$ 需要经过的次数, 操作: 
1. `1 a b k` 区间加,  $s_i += k,(i \in [a, b])$; 
2. `2 a b` 区间查询, 主人公从 $a$ 出发, 在 $b$ 停止, 不能出区间, 能否找到符合 $s_i(i \in [a,b])$ 的移动方案;

## 解题报告

稳稳地读错了题。。
NB地理解为随便从哪个点出发到随便哪个点停止, 我英死早!<br> 
实际上是很简单的题...<br> 
因为需要从一端到另一端, 所以先走一遍, 即$s_i -= 1, (i \in [l, r])$. <br> 
然后可以增加若干回头路, 也就是相邻的两个位置 $s -= 1$ . <br> 
如果不在区间上考虑, 把这个判断是否有解的子问题看做一个完整的问题的话, 显然有一个贪心的思路, 就是从前向后, 如果第 $i$ 个位置有剩余 $left_i$ , 则  $left_{i+1} -= left_{i}$  , 中途如果有 $left_i < 0$ 那么无解, 如果$ left_n != 0$ , 也就是没有消完, 那么无解;<br> 
向区间上考虑, 需要维护区间$left$的最小值和右端点的值;<br> 
再考虑合并, $a$ 和 $b$ 两个区间的合并, $a$ 的右端点会对 $b$ 的奇数位置偶数位置产生不同的影响, 所以需要把区间的奇数位置, 偶数位置分别记最小值, 并且需要记录区间的大小, 从而更新右端点的值;<br> 
解决区间修改的问题, 通过$left$的计算方式, 区间修改只会影响区间中的偶位置(从0标号), 需要标记和下传; <br> 

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

const int inf = 1000000007;
const int N = 800050; 

struct node { 
    ll siz, evmn, odmn, r, tag; 
    node() {siz = r = tag = 0, evmn = odmn = inf;} 
    node(ll siz, ll ev, ll od, ll r, ll tag)
        : siz(siz), evmn(ev), odmn(od), r(r), tag(tag){} 
    void operator += (ll k) { 
        if (siz & 1) r += k; 
        tag += k, evmn +=k;
    } 
    void print() { 
        cout << siz << ' ' << evmn << ' ' << odmn << ' ' << r <<  ' ' << tag << endl; 
    } 

} t[N], ans; 
node operator +(node a, node b) { 
    node c(a.siz + b.siz, 0, 0, 0, 0); 
    if (a.siz & 1) { 
        c.evmn = min(a.evmn, b.odmn + a.r); 
        c.odmn = min(a.odmn, b.evmn - a.r); 
        c.r = b.r + ((b.siz & 1)?-1:1) * a.r; 
    } else { 
        c.evmn = min(a.evmn, b.evmn - a.r); 
        c.odmn = min(a.odmn, b.odmn + a.r); 
        c.r = b.r + ((b.siz & 1)?-1:1) * a.r; 
    } 
    return c; 
} 
int n, q; 
inline void push_up(int x) { 
    t[x] = t[x << 1] + t[x << 1 | 1]; 
} 
inline void push_down(int x) { 
    if (t[x].tag) { 
        t[x << 1] += t[x].tag; 
        t[x << 1 | 1] += t[x].tag; 
        t[x].tag = 0; 
    } 
} 

void build(int x, int l, int r) { 
    if (l == r) { 
        t[x].siz = 1, cin >> t[x].evmn; 
        --t[x].evmn, t[x].r = t[x].evmn; 
        t[x].tag = 0; 
    } else { 
        int mid = (l + r) >> 1; 
        build(x << 1, l, mid); 
        build(x << 1 | 1, mid + 1, r); 
        push_up(x); 
    } 
} 
void add(int x, int l, int r, int L, int R, int v) { 
    if (L <= l && r <= R) { 
        t[x] += v; 
    } else { 
        int mid = (l + r) >> 1; push_down(x);
        if (L <= mid) add(x << 1, l, mid, L, R, v); 
        if (R > mid) add(x << 1 | 1, mid + 1, r, L, R, v); 
        push_up(x); 
    } 
} 
void query(int x, int l, int r, int L, int R) { 
    if (L <= l && r <= R) { 
        ans = ans + t[x]; 
    } else { 
        int mid = (l + r) >> 1; push_down(x); 
        if (L <= mid) query(x << 1, l, mid, L, R); 
        if (R > mid ) query(x << 1 | 1, mid + 1, r, L, R); 
    } 
}
int main() {
    ios :: sync_with_stdio(false); 
    cin >> n;
    build(1, 1, n); 
    int type, a, b, k; 
    cin >> q; 
    while (q--) { 
        cin >> type; 
        if (type == 1) { 
            cin >> a >> b >> k; 
            add(1, 1, n, a + 1, b + 1, k); 
        } else { 
            cin >> a >> b ; 
            ans = node(0, inf, inf, 0, 0); 
            query(1, 1, n, a + 1, b + 1); 
            if (ans.evmn >= 0 && ans.odmn >=0 && ans.r == 0) 
                puts("1"); 
            else 
                puts("0"); 
        } 
    } 
    return 0; 
} 
```
# 【CF 720D】Slalom
## [题目链接](http://codeforces.com/contest/720/problem/D)
## 题目大意
从 $(1, 1)$ 点到 $(m, n)$ 点的不同路径数, 定义两路径不同当且仅当存在一个矩形障碍物分别在这个路径的左右两侧; 
## 解题报告
这个题很像TA哥模拟题中的栅栏(之后补出题解); <br> 
利用差分, 令 $f[i][j]$ 表示 $(i, j)$ 位置比 $(i, j-1)$ 位置多出的方案数; <br> 
容易得到, $f[i][j]$ 会继承到 $f[i+1][j]$ , 如果$f[i+1][j]$没有障碍物;<br> 
考虑障碍物的作用: ①将当前行列位置的差分信息清空, ②使上一行的部分差分信息汇聚到障碍物下端右侧位置;<br> 
可以来个图感受一下:
![one](2.png)
白色的部分是障碍物, 那么红色位置的差分信息会汇聚到绿色位置; 
![two](3.png)
同样的, 这个图中, ①位置会分别对其上方的位置和绿色位置贡献一个差分+1; <br> 
总结一下, 对于当前行, 没有障碍物的位置, 继承上一行的差分, 障碍物最靠下位置差分信息清空, 右端汇聚前一行可以到达的差分信息;<br> 
继承信息, 不需要进行操作, 汇聚信息, 需要区间求和和单点修改, 清空信息, 需要区间清零;<br> 
所以线段树可以搞定; <br> 
查上一行可以到达的左端位置, 使用平衡树(set)很轻松; <br> 

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
#define sz(x) (int)(x).size()
#define mmst(a, x) memset(a, x, sizeof(a))

typedef long long ll;

typedef vector<int> vi;
typedef pair<int, int> pii;

const int p = 1000000007;
const int N = 1000100; 

vector <pii> in[N], out[N]; 
multiset<int> block; 
int sm[N << 2], n, m, k; 
bool zero[N << 2]; 
inline int stop(int &x) { 
    multiset<int> :: iterator i; 
    i = block.upper_bound(x); 
    if (i != block.begin()) {
        -- i; return *i + 1; 
    } else return 1; 
}  
inline void notify(int x) { 
    zero[x] = 1, sm[x] = 0; 
} 
inline void push_up(int x) { 
    sm[x] = (sm[x << 1] + sm[x << 1 | 1]) % p; 
}
inline void push_down(int x) { 
    if (zero[x]) notify(x << 1), notify(x << 1 | 1), zero[x] = 0; 
} 

void build(int x, int l, int r) { 
    if (l == r) sm[x] == 0; 
    else { 
        int mid = (l + r) >> 1; 
        build(x << 1, l, mid); 
        build(x << 1 | 1, mid + 1, r); 
        push_up(x); 
    } 
} 
void change(int x, int l, int r, int w, int v) { 
    if (l == r) sm[x] = v; 
    else { 
        int mid = (l + r) >> 1; 
        push_down(x); 
        if (w <= mid) change(x << 1, l, mid, w, v); 
        else change(x << 1 | 1, mid + 1, r, w, v); 
        push_up(x);
    }
}
void clear(int x, int l, int r, int L, int R) { 
    if (L <= l && r <= R) notify(x); 
    else { 
        int mid = (l + r) >> 1; 
        push_down(x); 
        if (L <= mid) clear(x << 1, l, mid, L, R); 
        if (R > mid) clear(x << 1 | 1, mid +1, r, L, R); 
        push_up(x); 
    } 
}
int query(int x, int l, int r, int L, int R) { 
    if (L <= l && r <= R) return sm[x]; 
    else { 
        int mid = (l + r) >> 1, ans = 0; 
        push_down(x); 
        if (L <= mid) (ans += query(x << 1, l, mid , L, R)) %= p; 
        if (R > mid) (ans += query(x << 1 | 1, mid + 1, r, L, R)) %= p; 
        return ans; 
    } 
} 

int main() { 
    ios :: sync_with_stdio(false); 
    cin >> n >> m >> k; 
    int x, y, _x, _y; 
    REP(i, k) { 
        cin >> x >> y >> _x >> _y; 
        in[y].pb( mp(x, _x) ); 
        out[_y + 1].pb( mp(x, _x) ); 
    } 
    build(1, 1, n); 
    change(1, 1, n, 1, 1); 
    REP(i, sz(in[1])) {
        block.insert(in[1][i].fi); 
        block.insert(in[1][i].se); 
    } 
    FORU(i, 2, m) { 
        vi delta; 
        REP(j, sz(in[i])) { 
            int r = in[i][j].se + 1; 
            int l = stop(r); 
            if (l <= r && r <= n)
                delta.pb(query(1, 1, n, l, r)); 
            else delta.pb(0);
        } 
        REP(j, sz(in[i])) { 
            int r = in[i][j].se + 1; 
            if (r <= n) 
                change(1, 1, n, r, delta[j]); 
        }
        REP(j, sz(in[i])) {
            clear(1, 1, n, in[i][j].fi, in[i][j].se); 
            block.insert(in[i][j].fi); 
            block.insert(in[i][j].se); 
        } 
        REP(j, sz(out[i])) { 
            block.erase(block.find(out[i][j].fi)); 
            block.erase(block.find(out[i][j].se)); 
        } 
    }
    int l = stop(n); 
    int ans = query(1, 1, n, l, n); 
    return cout << ans << endl, 0; 
} 
```

# 【CF 750G】New Year and Binary Tree Paths
## [原题链接](http://codeforces.com/contest/750/problem/G)
## 题目大意

无限大的一棵完全二叉树上, 节点编号和等于 $k$ 的路径条数;

## 解题报告

考虑路径的 $LCA$ 为 $x$ , 左右两条路径长度(不包括 $x$ )为 $len_a$
, $len_b$ , 路径编号和 $S$ 是多少; <br>
$$S = x * (2^{len_a+1} + 2^{len_b+1} - 3) + 2^{len_b} - 1 + bonus_a + bonus_b$$

一条自上而下的路径, 如果一个位置增加 $x$ , 从这个位置开始的长度为 $n$ , 那么增加 $x$ 的收益是 $x * (2^n - 1)$ ;<br>
所以 $S$ 的计算公式中的 $bonus_a$ 和 $bonus_b$ 代表左右两条长度分别为 $len_a - 1$ 和 $len_b - 1$ 的路径通过向右下方移动得到的收益;<br>
计算 $bonus$ 的时候, 因为$2^n -1$的形式不好统计方案, 而对每一个为一的二进制位获得的收益加一, 再对总收益/2, 正好可以得到向左移动为0, 向右移动为1对应的二进制数;<br>
所以得到 $bonus_a = a * 2 - ones_a$ ,其中 $a$ 表示左右移动得到的二进制数,  $bonus_b$ 同理; <br>
一个有效的性质是, <br>
$$bonus_a + bonus_b +2^{len_b} - 1 < 2^{len_a} - 2 + 2^{len_b  + 1} - 3 < 2^{len_a+1} + 2^{len_b+1} - 3$$

所以, 如果确定 $len_a$ 和 $len_b$ , 就可以直接计算 $LCA$ ;<br>
然后可以枚举 $ones_a$ 和 $ones_b$ , 将问题转化成, 两个二进制数, 分别有 $len_a - 1$ 和 $len_b - 1$ 位, $1$ 的个数分别为 $ones_a$ , $ones_b$ , 和为 $(bonus_a + bonus_b - ones_a - ones_b)/2$ 的方案数, 数位 $dp$ ; <br>

时间复杂度 $O(log^6S)$ ;

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
#define Tgetw(a, w) ((a >> w) & 1)

typedef long long ll;
typedef vector<int> vi;
typedef pair<int, int> pii;

const int N = 150;

ll s, ans;
inline ll dp(int ot, int la, int lb, ll sm) {
    static ll f[N][2], las[N][2];
    REP(a, ot + 1) REP(re, 2) f[a][re] = 0;
    f[0][0] = 1;
    for (int i = 0; (1LL << i) <= sm; ++i) {
        REP(a, ot + 1) REP(re, 2) {
            las[a][re] = f[a][re];
            f[a][re] = 0;
        }
        REP(ca, 2) if (i < la || ca == 0)
            REP(cb, 2) if (i < lb || cb == 0) {
                REP(re, 2) {
                    int dg = ca + cb + re;
                    if ((dg & 1) != Tgetw(sm, i))
                        continue;
                    for (int al = 0; al + ca + cb <= ot; ++al)
                        f[al + ca + cb][dg / 2] += las[al][re];
                }
        }
    }
    return f[ot][0];
}
int main() {
    cin >> s, ans = 0;
    REP(a, 59) REP(b, 59) {
        ll mul = (1LL << a + 1) + (1LL << b + 1) - 3;
        ll lca = s / mul ; if (!lca) continue;
        ll res = s - lca * mul - ( (1LL << b) - 1);
        REP(ones, a + b + 1) {
            ll aandb = res + ones;
            if (aandb < 0 || aandb % 2 == 1) continue;
            aandb /= 2;
            ans += dp(ones, max(a-1, 0), max(b-1, 0), aandb);
        }
    }
    cout << ans << endl;
    return 0;
}
```
# 【CF 702F】T-shirt
## [原题链接](http://codeforces.com/contest/702/problem/F)
## 题目大意
每个人会先买质量高的东西, 质量相同会先买便宜的东西, 给出 $n$ 个物品各自的质量和价格, 给出 $k$ 个人每个人带的钱, 问每个人能卖几个东西;
## 解题报告

这个裸数据结构题场上没人A, 你逗我?<br>
首先, 针对每个人, 计算他可以买哪些东西真是搞不了, 咋想都别扭. <br>
古人有句话说得好... 正难则反易.<br>
所以考虑每个东西被哪个人买好了..<br>
把人按照剩下的钱排个平衡树, 然后, 物品按照双关键字排序后, 从前向后扫, 对大于当前物品价格的平衡树上的区间, 做答案 $+1$ , 剩余钱数 $-price_i$ 的操作; <br>
有一个问题, 就是做完区间减法后, 平衡树的性质可能被破坏, 那么需要把操作后, 剩余钱数 $<price_i$ 的节点重新插入平衡树的合适位置;<br>
这样不会 $T$ 成狗 ? <br>
好像不会, 证明下:
1. 一个节点被重新加入平衡树的条件是 $price_i < money < 2price_i$;
2. 这个节点进行完重构操作后, $money' < price_i$ ;
3. $delta_{money} > money * \frac{1}{2}$ , 也就是每个节点最多进行 $log \ money$ 次操作;
4. 时间复杂度是 $O(nlog^2n)$ ;
<br>

写的是 `FHQ Treap`, 仍然经历了痛苦的打残->调试的过程; <br>

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
#define RAND (rand() << 15 | rand())
typedef long long ll;
typedef vector<int> vi;
typedef pair<int, int> pii;

const int N = 200100;

int n, c[N], q[N], lin[N], k, b[N], usr[N];
int ans[N];
struct node {
	int l, r, at, bt, sz, a, b, id, key;
	node(){}
	node(int b, int id, int key)
		: b(b), id(id), key(key){
			l = r = at = bt = a = 0, sz = 1;
		}
} t[N];
int root;
inline void in(int &x) {
	char ch = getchar(); int f= 1;
	for (; ch < '0' || ch > '9'; ch = getchar())
		if (ch == '-') f = -1;
	for (x = 0; ch >= '0' && ch <= '9'; ch = getchar())
		x = x * 10 + ch - 48;
	x *= f ;
}

inline bool cmp1(int x, int y) {
	return q[x] == q[y] ? c[x] < c[y] : q[x] > q[y];
}
inline bool cmp2(int x, int y) {
	return b[x] < b[y];
}

inline void add_a(int x, int v) {
	t[x].a += v, t[x].at += v;
}
inline void add_b(int x, int v) {
	t[x].b += v, t[x].bt += v;
}
inline void push_up(int x) {
	t[x].sz = 1;
	if (t[x].l) t[x].sz += t[t[x].l].sz;
	if (t[x].r) t[x].sz += t[t[x].r].sz;
}
inline void push_down(int x) {
	if (t[x].at) {
		if (t[x].l) add_a(t[x].l, t[x].at);
		if (t[x].r) add_a(t[x].r, t[x].at);
		t[x].at  = 0;
	}
	if (t[x].bt) {
		if (t[x].l) add_b(t[x].l, t[x].bt);
		if (t[x].r) add_b(t[x].r, t[x].bt);
		t[x].bt = 0;
	}
}
inline int build(int *lis) {
	int x ;
	static int stack[N], top, last;
	FORU(i, 1, k) {
		x = lis[i], last = 0;
		t[i] = node(b[x], x, RAND);
		while (top && t[stack[top]].key > t[i].key) {
			last = stack[top], push_up(stack[top]);
			stack[top--] = 0;
		}
		if (top) t[stack[top]].r = i;
		stack[++top] = i, t[i].l = last;
	}
	while (top) push_up(stack[top--]);
	return stack[1];
}
pii split(int x, int k) {
	if (x == 0) return mp(0, 0);
    push_down(x); pii y;
	if (k == 0) return mp(0, x);
	int lsize = t[t[x].l].sz;
	if (lsize >= k) {
		y = split(t[x].l, k);
		t[x].l = y.se, push_up(x);
		y.se = x;
	} else {
		y = split(t[x].r, k - lsize  - 1);
		t[x].r = y.fi, push_up(x);
		y.fi = x;
	}
	return y;
}
int merge(int a, int b) {
	if (a * b == 0) return a + b;
	push_down(a), push_down(b);
	if (t[a].key < t[b].key) {
		t[a].r = merge(t[a].r, b);
		push_up(a); return a;
	} else {
		t[b].l = merge(a, t[b].l);
		push_up(b); return b;
	}
}
inline int lower(int rt, int v) {
	int x = rt, tmp = 0;
	while (x) {
        push_down(x);
		if (t[x].b < v) {
			tmp += t[t[x].l].sz + 1;
			x = t[x].r;
		} else x = t[x].l;
	}
	return tmp;
}
inline void insert(int x, int &rt) {
	int k = lower(rt, t[x].b);
	pii droot = split(rt, k);
	droot.se = merge(x, droot.se);
	rt = merge(droot.fi, droot.se);
}
void del(int x, int &to) {
	push_down(x);
	if (t[x].l) del(t[x].l, to);
	if (t[x].r) del(t[x].r, to);
	t[x].l = t[x].r = 0;
    insert(x, to);
}
inline void buy(int price) {
	int k = lower(root, price);
	pii droot = split(root, k);
	add_b(droot.se, -price);
	add_a(droot.se, 1);
	int _k = lower(droot.se, price);
	pii dright = split(droot.se, _k) ;
	del(dright.fi, droot.fi);
	root = merge(droot.fi, dright.se);
}
void query(int x) {
	ans[ t[x].id ] = t[x].a;
	push_down(x);
	if (t[x].l) query(t[x].l);
	if (t[x].r) query(t[x].r);
}
inline void itera(int x) {
    push_down(x);
    if (t[x].l) itera(t[x].l);
    cout << x << " , " << t[x].id << " : " << t[x].a << " , " << t[x].b << endl;
    if (t[x].r) itera(t[x].r);
}

int main() {
	srand(time(0) + 217);
	in(n);
	FORU(i, 1, n) in(c[i]), in(q[i]), lin[i] = i;
	in(k);
	FORU(i, 1, k) in(b[i]), usr[i] = i;
	sort(lin + 1, lin + 1 + n, cmp1);
	sort(usr + 1, usr + 1 + k, cmp2);
	root = build(usr);
//    itera(root);
	FORU(i, 1, n) buy(c[lin[i]]);// itera(root);
	query(root);
	FORU(i, 1, k)
		printf("%d%c", ans[i], (i == k)?'\n':' ');
	return 0;
}
```

# 【CF 678D】Lena and Queries

## [原题链接](http://codeforces.com/contest/678/problem/F)
## 题目大意
三个操作：
1. 加入一个有序数对$(a, b)$; 
2. 查询对于给出的$x$, $ax+b$的最大值;
3. 取消$i$号插入操作;

##解题报告

开始想的是时间分治, 离线后, 每个插入相当于对一个区间有效;<br>
每次将完全包含当前区间的插入点做成一个凸包, 对区间内的询问点进行查询;<br> 
算一下这样的时间复杂度, 每个区间可以表示成 $logn$ 的分治区间, 也就是被插入 $logn$ 次凸包; 询问需要在凸包上三分, 每个询问需要三分 $logn$ 次, 令三分的复杂度是 $T$ , 那么时间复杂度为 $O(nTlogn)$ ; <br>
对于三分复杂度这个问题, 我用的是三等分的做法, 没有去借鉴网上的 $O(3log_{3}n)$ 的说法, 鉴于每次问题规模缩为原问题的 $\frac{2}{3}$, 每次需要运算两个点的值, 所以我把复杂度写成 $O(2log_{\frac{3}{2}}n)$ , 这样看就是比二分慢不少, 貌似挺科学?<br>
所以这样做的复杂度就是 $O(nlognlog_{\frac{3}{2}}n)$ ; <br>
要写的时候, 突然想到向量集那个题, 受到一些启发, 干脆把时间建成线段树, 然后把点放到线段树上去, 在线段树上查询;<br> 
使用标记永久化, 每个点会加入 $logn$ 次, 复杂度是不变的; <br> 
这样空间翻倍了, 但写得更板子...<br> 
A了之后看了眼题解, woc, 这个分块有点科学呀, 而且复杂度根本不差!!!<br>
所以还是应该继续思考...<br>

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
#define sz(x) (int)(x).size()

typedef long long ll;
typedef vector<int> vi;
typedef pair<int, int> pii;

const int N = 301000;
inline void in(int & x) {
    char ch = getchar();int f = 1;
    for (; ch < '0' || ch > '9'; ch = getchar())
    if (ch == '-') f = -1;
    for (x = 0; ch >= '0' && ch <= '9'; ch = getchar())
    x = x * 10 + ch - 48;
    x *= f;
}
ll cross(pll a, pll b, pll c) {
    ll x = a.fi - b.fi, y = a.se - b.se;
    ll _x = b.fi - c.fi, _y = b.se - c.se;
    return x * _y - y * _x;
}
ll val(pll a, ll x) {
    return x * a.fi + a.se;
}

struct convex {
    vector< pll > po;
    void add(pll x) {po.pb(x);}
    void build() {
        vector< pll > tp; tp = po;
        po.clear();
        sort(tp.begin(), tp.end());
        REP(i, sz(tp)) {
            while (sz(po) >= 2 && cross(tp[i], po.back(), po[sz(po)-2]) <= 0)
                po.pop_back();
            po.pb(tp[i]);
        }

    }
    ll query(ll x) {
        int l = 0, r = sz(po) -1;
        ll ans = LLONG_MIN;
        while (l + 2 < r) {
            int _l = l + (r - l + 1) / 3;
            int _r = r - (r - l + 1) / 3;
            if (val(po[_l], x) > val(po[_r], x))
                r = _r;
            else l = _l;
        }
        FORU(i, l, r) ans = max(ans, val(po[i], x));
        return ans;
    }
} c[N << 2];


int t[N], l[N], r[N], q[N], n;
pll poi[N];
ll ans;
void add(int x, int l, int r, int _l, int _r, pll p) {
    if (_l <= l && r <= _r) {
        // cout << l << ' ' << r << ' '<< _l << ' '<< _r << endl;
        c[x].add(p); return;
    }
    if (l == r) return;
    int mid = (l + r) >> 1;
    if (_l <= mid) add(x << 1, l, mid, _l, _r, p);
    if (_r > mid) add(x << 1 | 1, mid + 1, r, _l, _r, p);
}
void build(int x, int l, int r) {
     c[x].build();
     if (l == r) return;
     int mid = l + r >> 1;
     build(x << 1, l, mid);
     build(x << 1 | 1, mid + 1, r);
 }
ll query(int x, int l, int r, int p, ll _x) {
    ll ans = LLONG_MIN;
    ans = max(ans, c[x].query(_x));
    if (l == r) return ans;
    int mid = (l + r) >> 1;
    if (p <= mid)
        ans = max(ans, query(x << 1, l, mid, p, _x));
    else
        ans = max(ans, query(x << 1 | 1, mid + 1, r, p, _x));
    return ans;
}
int main() {
    in(n); int x, y;
    FORU(i, 1, n) {
        in(t[i]);
        if (t[i] == 1) {
            in(x), in(y), poi[i] = mp(x, y);
            l[i] = i, r[i] = n;
        } else if (t[i] == 2)
                in(x), r[x] = i;
            else in(q[i]);
    }
    // cout << endl;
    FORU(i, 1, n) if (t[i] == 1) {
//        cout << i << ": " <<  l[i] << ' ' << r[i] << endl;
        add(1, 1, n, l[i], r[i], poi[i]);
    }
    build(1, 1, n);
    FORU(i, 1, n) if (t[i] == 3) {
        ans = query(1, 1, n, i, q[i]);
        if (ans == LLONG_MIN) puts("EMPTY SET");
        else printf("%lld\n", ans);
    }
    return 0;
}
```

# 【CF 700D】 Huffman Coding on Segment
## [原题链接:](http://codeforces.com/contest/700/problem/D)
## 题目大意
区间查询huffman编码
## 解题报告
学习huffman编码/树的美好童年时光QWQ.<img src="4.jpg" alt="简直一派胡言" width="80" height="80">

利用huffman编码的定义, 每次最后一个0/1位相当于把两个数合并, 类似于果子合并中的合并; <br>
显然, 出现次数小的应该先进行合并, 这样可以减小总的编码长度, 答案是所有数的huffman编码长度和, 也就是编码过程中huffman树的每个点的深度和;<br> 
所以问题实际上是对一个区间内, 每个数的出现次数为权值构建huffman树...(<del>怎么神犇根本不需要想就知道呀...</del>)<br>
但是一个直接每个区间构建huffman树是 $O(n^2logn)$ 的...<img src="5.jpg" width="80" height="70" >

考虑优化, 貌似出现次数相同的果子(...)可以一起处理, 也就是如果出现次数是 $times$ 的果子一共有 $k$ 个, 那么如果当前最小出现次数果子的出现次数就是 $times$ , 就可以直接合并出 $\lfloor k/2 \rfloor$ 个 $times * 2$, 如果还剩下一个果子, 就向后面寻找配对的果子...<br> 
这个复杂度就很科学了! 因为不同的出现次数只有 $\sqrt{n}$ 种, 而类似桶排, 可以对 $times <= O(sqrt{n})$ 做这样的处理, 使得剩下的果子出现次数 $times > \sqrt{n}$ ... 而这样的果子最多 $n / \sqrt{n}$ 个... <img src="6.gif">

现在的问题是, 怎样统计区间内每个数出现的次数和出现相应次数的有几个数, 连我这个蒟蒻都会的莫队算法可以搞呀..<br> 
妙! 复杂度是 $O(n \sqrt{n} + n \sqrt{n} + n \sqrt{n}log \sqrt{n})$  <br>

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

const int N = 100100;

int car[N], tim[N], n, q, a[N];
int id[N], B, blo[N], l[N], r[N] ;
ll ans[N];
vi large;
inline bool cmp(int a, int b) {
    return blo[l[a]] == blo[l[b]] ? r[a] < r[b]
        : blo[l[a]] < blo[l[b]];
}
inline void colorify(int w, int v) {
    tim[car[w]] --, car[w] += v, tim[car[w]] ++;
}
inline void got_ans(ll &ans) {
    ans = 0;
    static int tis[N];
    FORU(i, 1, B * 2) tis[i] = tim[i];
    FORU(i, 1, B) {
        ans += tis[i] / 2 * 2 * i;
        tis[i * 2] += tis[i] / 2;
        tis[i] = tis[i] % 2;
        if (tis[i]) FORU(j, i + 1, B) {
            if (tis[j]) {
                ans += i + j;
                tis[i] --, tis[j]--;
                tis[i + j] ++;
                break;
            }
        }
    }
    priority_queue<int, vi, greater<int> > q;
    FORU(i, 1, B * 2) REP(j, tis[i]) q.push(i);
    REP(i, sz(large))
        if (car[large[i]] > B * 2)
            q.push(car[large[i]]);
    int x, y;
    while (q.size() >= 2) {
        x = q.top(), q.pop();
        y = q.top(), q.pop();
        ans += x + y;
        q.push(x + y);
    }
}
int main() {
    ios :: sync_with_stdio(false);
    cin >> n , B = pow(2 * n, 0.5);
    FORU(i, 1, n) cin >> a[i], ++car[a[i]];
    FORU(i, 1, 100000) if (car[i] > B * 2) large.pb(i);
    mmst(car, 0);
    int tmp = 0;
    FORU(i, 1, n) if (i % B == 0)
        blo[i] = tmp ++;
    else blo[i] = tmp;
    cin >> q;
    FORU(i, 1, q) id[i] = i;
    FORU(i, 1, q) cin >> l[i] >> r[i];
    sort(id + 1, id + 1 + q, cmp);
    int _l = 1, _r = 0, __l, __r, x;
    FORU(i, 1, q) {
        x = id[i], __l = l[x], __r = r[x];
        while (_r < __r)
            ++_r, colorify(a[_r], 1);
        while (_l > __l)
            --_l, colorify(a[_l], 1);
        while (_r > __r)
            colorify(a[_r], -1), --_r;
        while (_l < __l)
            colorify(a[_l], -1), ++_l;
        got_ans(ans[x]);
    }
    FORU(i, 1, q) cout << ans[i] << endl;
    return 0;
}
```

#【CF 755G】PolandBall and Many Other Balls
## [题目链接](http://codeforces.com/contest/755/problem/G)
## 题目大意

一排石子, 相邻的一块\两块可以圈成一个单位, 问$n$个石子, 圈出$k$个单位的方案数. 

## 解题报告

很容易想到**倍增**, 分为两种情况考虑 : 

1. 当前石子数$x$, 变为$2x$ , 对于选出$k$个单位的方案, 显然可以分成左边选$j$个, 右边选$k-j$个两部分. 但是有一种特殊的情况, 就是就是将第$x$块和第$x+1$块圈成一块, 那么方案数还需要加上在左边$(x-1)$块中选$j$, 右边$(x-1)$块中选$k-j$这种情况. 
2. 当前石子数为$x$, 变为$x+1$, 这个是简单的dp, $f[x][k] = f[x-1][k]+f[x-1][k-1]+f[x-2][k-1]$ 使用当前的$x$和$x-1$的信息尽可以转移. 

上面的两种情况有一个问题, 就是当$x->2x$ , 没有统计出$2x-1$的信息. 

考虑$2x-1 = x+x-1$ , 也就是左边$x$和右边$x-1$或者左边$x-1$右边$x$ . 但是这两个情况又有重叠, 其实就是第一种情况左边的$x$将最后一个位置单独圈成一个单位, 和第二种情况右边的$x$将第一个位置单独圈成一个单位后, 得到的答案是重复的, 所以需要减去重复的答案.  

整理一下, 可以得到完整的转移: 
$$
\left.
\begin{array}{l}
&f[x][k]=f[x-1][k]+f[x-1][k-1]+f[x-2][k-1]\\
&f[2x] (k) = (f[x] \otimes f[x])(k) + (f[x-1] \otimes f[x-1])(k-1)\\
&f[2x-1] (k) = 2(f[x] \otimes f[x-1])(k) - (f[x-1] \otimes f[x-1])(k-1)
\end{array}
\right.
$$

NTT, 大力优化一波卷积. 

## 代码

```c++
#include <bits/stdc++.h> 

using namespace std;

#define FORU(i, a, b) for (int i = int(a), nn = int(b); i <= nn; ++i) 
#define FORD(i, a, b) for (int i = int(a), nn = int(b); i >= nn; --i) 
#define REP(i, b) for (int i = 0, nn = int(b); i < b; ++i) 
#define DEBUG(x) cout << (#x) << " : " << x << endl

typedef long long ll; 
typedef double ff; 

const int K = 1 << 16; 
const int p = 998244353; 

int n, k, _n, re[K]; 
ll f[K], g[K]; 

inline ll _fast(ll x, int k) {
	ll ans = 1; 
	for (; k; k >>= 1) { 
		if (k & 1) ans = ans * x % p; 
		x = x * x % p ; 
	} 
	return ans; 
} 

inline int ad(int x) {
 	if (x >= p) x -= p ; 
 	if (x < 0) x += p; 
 	return x; 
}
  
inline void fft(ll *a, int f = 1) {
 	REP(i, _n) if (i < re[i]) 
 		swap(a[i], a[ re[i] ]); 
 	ll G = (f == 1) ? 3 : _fast(3, p-2); 
 	for (int m=1; m < _n; m <<= 1) { 
 		ll wn = _fast(G, (p-1) / m / 2); 
 		for (int i=0; i < _n; i += m << 1) { 
 			ll w = 1; 
 			for (int j=0; j < m; ++j) {
 				ll x = a[i+j], y = a[i+j+m] * w % p; 
 				a[i+j] = ad(x + y), a[i+j+m] = ad(x - y); 
 				w = w * wn % p;
 			} 
 		} 
 	} 
 	if (f == -1) { 
	 	ll po = _fast(_n, p-2); 
	 	REP(i, _n) a[i] = a[i] * po % p; 
	}
 	
} 

inline void one() {
 	FORD(i, k, 1) {
 	 	g[i] = f[i]; 
 	 	f[i] = ad(ad(f[i] + f[i-1]) + g[i-1]); 
 	}
 	g[0] = f[0]; 
}  

inline void two() {
 	fft(f), fft(g); 
 	ll w = 1, wn = _fast(3, (p-1)/_n), _f, _g; 
 	REP(i, _n) { 
 		_f = (f[i] * f[i] + g[i] * g[i] % p * w) % p; 
 		_g = (f[i] * g[i] * 2 - (g[i] * g[i]) % p * (w + 1)) % p; 
 		f[i] = _f, g[i] = ad(_g), w = w * wn % p; 
 	} 
 	fft(f, -1), fft(g, -1); 
 	FORU(i, k + 1, _n - 1) f[i] = g[i] = 0; 
} 

int main() {
 	
 	ios :: sync_with_stdio(false); 
	cin >> n >> k;  
	for (_n = 1; _n < k * 2; _n <<= 1); 

	for (int i = 0, j = 0; i < _n; ++i) {
 		re[i] = j; 
 		for (int _k = _n>>1; (j^=_k) < _k; _k >>= 1); 

 	} 
 	
	f[0] = 1; 	
	for (int i = 1<<29, s=0; i; i >>= 1) {
	 	if (s) two() , s <<= 1; 
	 	if (n & i) one(), ++ s; 
	} 
	FORU(i, 1, k) 
		cout << f[i] <<( (i==k)?'\n':' '); 
	return 0; 
}	
```
