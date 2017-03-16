---
title: codechef-选做
date: 2017-03-14 10:59:26
tags:
  - 平衡树
  - 矩阵/行列式
  - bitset
  - DP
  - 网络流
  - 树状数组
  - 分块
  - 树链剖分
  - 并查集
  - 数论
categories:
  - 题目集锦
  
---

> 选做CC hard 题中最简单，通过人数最多的做...

<!--more-->

# 【CC MGCHGYM】Misha and Gym
## 题目大意

$[1,n]$ 每个位置有一个重量$w_i$ , 支持三个操作： 

1. 修改一个位置的重量
2. 翻转一个区间
3. 查询一个区间内的重量能否拼出重量$w$ . 

保证出现的重量不超过10种， 3操作只有$1000$个。 


## 解题报告

很trivial的一个题呀...

因为出现的重量很少， 所以每个区间维护每个重量出现的次数，  因为存在翻转操作， 所以需要大力的写一个平衡树, 手选FHQTreap.

现在就差操作三了， 因为次数只有1000， 所以试着大力的背包一下， 咦 -_- 复杂度是$O(wn)$的， 好像不行呀...

优化背包的办法好多呀， 因为维护0\1信息， 所以强上bitset呀， 复杂度就是$O(wn/64)$了 . 有理想的交了好几发.都T了...

所以继续优化背包， 因为种类很少， 所以每一种用二进制优化一下， 复杂度就是$O(10w*log(n)/64)$了， 好像还是不科学， 但是竟然A了 。 WLGC。。。

## 代码

```c++
#include <bits/stdc++.h> 
#pragma GCC optimize("O3")
using namespace std;
 -
		} else { 
			go[b][0] = merge(a, go[b][0]); 
			push_up(b); 
			return b; 
		} 
	} 
		
			
}
 
inline void in(int &x) { 
	for (; *cp < '0' || *cp > '9'; cp++) ; 
	for (x = 0; *cp >= '0' && *cp <= '9'; cp++) 
		x = x * 10 + *cp - 48; 
} 
 
inline void enter(int x) { 
	if (id[x]==-1) {
		id[x] = tot ++; 
		v.push_back(x); 
	} 
	
} 
 
inline void chang() { 
	using namespace FHQ; 
	int l, w; in(l), in(w); 	
	droot le = split(root, l-1); 
	droot ri = split(le.se, 1); 
	int x = ri.fi; 
	enter(w); 
	va[x] = id[w], info[x] = data(va[x]);
	ra[x] = RAND;
	le.se = merge(x, ri.se); 
	root = merge(le.fi, le.se); 
} 
 
inline void revers() {
	using namespace FHQ;  
	int l, r; in(l), in(r); 
	droot le = split(root, l - 1); 
	droot ri = split(le.se, r - l + 1); 
	_rev(ri.fi) ; 
	le.se = merge(ri.fi, ri.se); 
	root = merge(le.fi, le.se); 
} 
 
inline void query() { 
	using namespace FHQ; 
	int l, r, w; in(l), in(r), in(w); 
	droot le = split(root, l - 1); 
	droot ri = split(le.se, r-l+1); 
	int x= ri.fi; bool check=0; 
	dp.reset(); dp.set(0); 
	for (int i = 0; i < tot; ++i) {
		int lef = info[x].sts[i]; 
		for (int j = 1; j <= lef; j <<= 1) {
			dp |= dp << (j * v[i]); 
			lef -= j, check |= dp.test(w); 
			if (check) break; 
		} 
		if (lef) dp |= dp << (lef * v[i]); 
		check |= dp.test(w); 
		if (check) break; 
	} 
	
		
//	cout << ", " << w << ": "; 
	puts((check) ? "Yes" : "No"); 
	le.se = merge(ri.fi, ri.se); 
	root = merge(le.fi, le.se); 
} 
 
//--- line to be a line ---- 
 
 
int main() { 
	fread(cp, 1, 2000000, stdin); 
	in(n), in(q);
	memset(id, -1, sizeof(id));  
	FORU(i, 1, n)  
		in(_w[i]), enter(_w[i]); 
	srand( time(0) + 217);   
	FHQ :: build() ; 
	
	int type; 
	while (q --) {
		in(type); 
		switch (type) { 
		case 1 : chang(); break; 
		case 2 : revers();break; 
		case 3 : query(); break; 
		} 
	} 
	return 0; 
} 
```
# 【CC BWGAME】Black-white Board Game

## 题目大意

一个棋盘， 第$i$行的$[l_i, r_i]$被染成黑色， 两个人Alex和Fedor， 每个人选择排列$P_i$, 要求$[i,P_i]$为黑色， 且对Alex, $P_i$的逆序对个数为偶数， 对Fedor， 逆序对个数为奇数， 问谁能选的排列多.

## 解题报告

逆序对个数的奇偶？ 可以联想到行列式的一个求法 $\sum (-1)^{\iota(P_i)} \prod_{j=1}^{n} a[j][P_i[j]]$ ; 

所以就是把棋盘看成一个O\1矩阵， 对其行列式求值。

简单的方法是高斯消元， 但是$O(n^3)$显然过不了， 考虑矩阵中$1$的分布的特殊性，就是每一行， $1$的分布是一个连续的区间。

可以将$[l_i,r_i]$放入第$l_i$棵左偏树， 根据$r_i$建小根堆， 枚举到第$i$列， 拿出第$i$棵左偏树中$r_i$最小的区间， 使用这个区间消元这棵树中剩余的区间， 并加入$r_i+1$这个树 ， 并根据是否交换两行改变行列式的值。

这样消元， 对角线一定全是$1$(或者有$0$, 也就是平手的情况), 只需要判断行列式的值的正负就好。 

## 代码

不知道为什么CC这道题不能提交， 所以只能当过样例选手了。

```c++
#include <bits/stdc++.h> 

using namespace std;

#define FORU(i, a, b) for (int i = int(a), nn = int(b); i <= nn; ++i) 
#define FORD(i, a, b) for (int i = int(a), nn = int(b); i >= nn; --i) 
#define REP(i, b) for (int i = 0, nn = int(b); i < b; ++i) 

typedef long long ll; 
typedef double ff; 

using namespace std; 

const int N = 200010; 
int test, n, rt[N], ans, tot, a[N]; 

struct LTREE {
	int de, ls, rs, v, id; 
	LTREE() {de=ls=rs=v=0;}
	LTREE(int v, int id) :v(v),id(id) {
		de = ls = rs = 0; 
	} 
} no[N]; 

inline void in(int &x) {
 	char ch = getchar(); 
 	for (;ch < '0' || ch > '9'; ch=getchar()); 
 	for (x=0; ch>='0' && ch<='9'; ch=getchar())
 		x = x * 10 + ch - 48; 
 } 
 
inline int birth(int v, int x) {
 	++ tot; no[tot] = LTREE(v, x); 
 	return tot; 
} 
 
int merge(int a, int b) {
 	if (!(a * b)) return a + b; 
 	if (no[a].v > no[b].v) swap(a, b); 
 	no[a].rs = merge(no[a].rs, b); 
 	if (no[no[a].rs].de > no[no[a].ls].de) 
 		swap(no[a].ls, no[a].rs); 
 	no[a].de = no[no[a].rs].de + 1; 
 	return a; 
} 

inline void ins(int &rt, int v, int x) {
 	rt = merge(rt, birth(v, x)); 
}

inline int top(int &rt) { 
	if (!rt) return 0; 
	int tmp = rt; 
	rt = merge(no[tmp].ls, no[tmp].rs); 
	return tmp; 
} 

int main() {
//	freopen("A.in", "r", stdin); 
	in(test);
	while (test --) {
	 	in(n); int l, r; ans = 1; 
	 	for (int i = 1; i<=n; ++i) 
	 		in(l), in(r), ins(rt[l], r, i), a[i] = i; 
	 		
	 	for (int i = 1; i<=n; ++i) 
	 		if (!rt[i]) { ans = 0; break;} 
	 		else { 
	 			int x = top(rt[i]), y = top(rt[i]); 
//	 			cout << i << ", " << no[x].v << ": " << no[x].id << endl;
	 			if (no[x].v == no[y].v) { 
	 				ans = 0; break; 
	 			} else rt[i] = merge(rt[i], y); 
	 			if (no[x].id != a[i]) 
	 				ans *= -1, a[no[x].id]=a[i];
	 			rt[no[x].v+1] = merge(rt[no[x].v+1], rt[i]); 
	 		} 
	 	(ans == 0) ? puts("Draw")
	 	: (ans == 1) ? puts("Alex")
	 		: puts("Fedor"); 
	 }
	 return 0; 
}
```
#【CC CBAL】Chef and Balanced Strings
## 题目大意

区间内， 每个字符出现偶数次的子区间的个数 。

## 解题报告

这个破玩意显然只能分块乱搞呀 ...

首先， 字符集很小， 可以把每个字符出现的前缀奇偶性压成一个int, 离散后搞成一个状态信息.

大力预处理， 从第$i$个块开始，向后到第$j$个位置的三类答案$fans_{1|2|3}$ , 从第$i$个块开始， 向前到第$j$个位置的三类答案$bans_{1|2|3}$ , 这个可以通过扫的同时记录每个奇偶状态信息的出现次数， 出现位置$i$的和， 出现次数$i^2$的和 得到。 

然后， 对于每个查询， 设$l$所在的块为$Li$, 从$Lm$开始， $r$所在的块为$Ri$, 到$Rm$结束， 可以先得到左右的零散位置和中间若干整块的答案加上整块中的答案， 就差左右零散位置之间的答案， 这部分暴力扫一遍就好. 

时空复杂度$O(n\sqrt{n})$  .

## 代码

写得真·TM爽 。

```c++
#include <bits/stdc++.h> 

using namespace std;

#define FORU(i, a, b) for (int i = int(a), nn = int(b); i <= nn; ++i) 
#define FORD(i, a, b) for (int i = int(a), nn = int(b); i >= nn; --i) 
#define REP(i, b) for (int i = 0, nn = int(b); i < b; ++i) 

typedef long long ll; 
typedef double ff; 


const int N = 100010; 
const int SN = 250; 

char s[N];

// forward ans
ll fans0[SN][N], fans1[SN][N], fans2[SN][N];

// backward ans
ll bans0[SN][N], bans1[SN][N], bans2[SN][N];

int n, Q, S, idx[1<<26], values[N];

ll sums0[N], sums1[N], sums2[N];

int cnt;

int get_idx(int x) {
    if (!~idx[x]) idx[x] = cnt++;
    return idx[x];
}

ll ans0, ans1, ans2;

void initalize() {
    S = max(n/200, int(n/sqrt(Q)));
    int x;
    idx[x = 0] = -1;
    for (int i = 0; i < n; i++) {
        idx[x ^= 1 << s[i] - 'a'] = -1;
    }
    
    cnt = 0;
    values[0] = get_idx(x = 0);
    for (int i = 0; i < n; i++) {
        values[i+1] = get_idx(x^=(1<<s[i]-'a'));
    }

    for (int i = 0; i*S <= n; i++) {
        
        for (int j = i*S; j <= n; j++) {
            int v = values[j];
            sums0[v] = 0;
            sums1[v] = 0;
            sums2[v] = 0;
        }
        ans0 = ans1 = ans2 = 0;
        for (int j = i*S; j <= n; j++) {
            int v = values[j];
            fans0[i][j] = ans0 += sums0[v];
            fans1[i][j] = ans1 += sums0[v]*j - sums1[v];
            fans2[i][j] = ans2 += sums0[v]*j*j - sums1[v]*2*j + sums2[v];
            sums0[v]++;
            sums1[v] += j;
            sums2[v] += j*(ll)j;
        }

        int end = min(n, i*S + S - 1);
        for (int j = end; j >= 0; j--) {
            int v = values[j];
            sums0[v] = 0;
            sums1[v] = 0;
            sums2[v] = 0;
        }
        ans0 = ans1 = ans2 = 0;
        for (int j = end; j >= 0; j--) {
            int v = values[j];
            bans0[i][j] = ans0 += sums0[v];
            bans1[i][j] = ans1 += sums1[v] - sums0[v]*j;
            bans2[i][j] = ans2 += sums0[v]*j*j - sums1[v]*2*j + sums2[v];
            sums0[v]++;
            sums1[v] += j;
            sums2[v] += j*(ll)j;
        }
    }
}




ll solve0(int L, int R) {
    L--;
    int Li = (L+1)/S, Ri = (R-1)/S;
    if (Ri - Li <= 1) {
        // normal
        for (int j = L; j <= R; j++) {
            int v = values[j];
            sums0[v] = 0;
        }
        ans0 = 0;
        for (int j = L; j <= R; j++) {
            int v = values[j];
            ans0 += sums0[v];
            sums0[v]++;
        }
    } else {
        int Lm = Li*S+S-1, Rm = Ri*S;
        ans0 = fans0[Li+1][R] + bans0[Ri-1][L] - fans0[Li+1][Rm-1];
        for (int j = L; j <= Lm; j++) {
            int v = values[j];
            sums0[v] = 0;
        }
        for (int j = Rm; j <= R; j++) {
            int v = values[j];
            sums0[v] = 0;
        }
        for (int j = L; j <= Lm; j++) {
            int v = values[j];
            sums0[v]++;
        }
        for (int j = Rm; j <= R; j++) {
            int v = values[j];
            ans0 += sums0[v];
        }
    }
    return ans0;
}

ll solve1(int L, int R) {
    L--;
    int Li = (L+1)/S, Ri = (R-1)/S;
    if (Ri - Li <= 1) {
        // normal
        for (int j = L; j <= R; j++) {
            int v = values[j];
            sums0[v] = 0;
            sums1[v] = 0;
        }
        ans1 = 0;
        for (int j = L; j <= R; j++) {
            int v = values[j];
            ans1 += sums0[v]*j - sums1[v];
            sums0[v]++;
            sums1[v] += j;
        }
    } else {
        int Lm = Li*S+S-1, Rm = Ri*S;
        ans1 = fans1[Li+1][R] + bans1[Ri-1][L] - fans1[Li+1][Rm-1];
        for (int j = L; j <= Lm; j++) {
            int v = values[j];
            sums0[v] = 0;
            sums1[v] = 0;
        }
        for (int j = Rm; j <= R; j++) {
            int v = values[j];
            sums0[v] = 0;
            sums1[v] = 0;
        }
        for (int j = L; j <= Lm; j++) {
            int v = values[j];
            sums0[v]++;
            sums1[v] += j;
        }
        for (int j = Rm; j <= R; j++) {
            int v = values[j];
            ans1 += sums0[v]*j - sums1[v];
        }
    }
    return ans1;
}

ll solve2(int L, int R) {
    L--;
    int Li = (L+1)/S, Ri = (R-1)/S;
    if (Ri - Li <= 1) {
        // normal
        for (int j = L; j <= R; j++) {
            int v = values[j];
            sums0[v] = 0;
            sums1[v] = 0;
            sums2[v] = 0;
        }
        ans2 = 0;
        for (int j = L; j <= R; j++) {
            int v = values[j];
            ans2 += sums0[v]*j*j - sums1[v]*2*j + sums2[v];
            sums0[v]++;
            sums1[v] += j;
            sums2[v] += j*(ll)j;
        }
    } else {
        int Lm = Li*S+S-1, Rm = Ri*S;
        ans2 = fans2[Li+1][R] + bans2[Ri-1][L] - fans2[Li+1][Rm-1];
        for (int j = L; j <= Lm; j++) {
            int v = values[j];
            sums0[v] = 0;
            sums1[v] = 0;
            sums2[v] = 0;
        }
        for (int j = Rm; j <= R; j++) {
            int v = values[j];
            sums0[v] = 0;
            sums1[v] = 0;
            sums2[v] = 0;
        }
        for (int j = L; j <= Lm; j++) {
            int v = values[j];
            sums0[v]++;
            sums1[v] += j;
            sums2[v] += j*(ll)j;
        }
        for (int j = Rm; j <= R; j++) {
            int v = values[j];
            ans2 += sums0[v]*j*j - sums1[v]*2*j + sums2[v];
        }
    }
    return ans2;
}


void decode(){
    ll A = 0, B = 0, ans=0 ;
    while( Q -- ) {
        int X, Y, type;
        scanf("%d%d%d", &X, &Y, &type);
        int L = ( X + A ) % n + 1; 
        int R = ( Y + B ) % n + 1; 
        if (L > R) swap(L, R);
        if (type == 0) ans = solve0(L, R); 
        if (type == 1) ans = solve1(L, R); 
        if (type == 2) ans = solve2(L, R); 
        
        printf("%lld\n", ans);
        A = B, B = ans;   
    } 
}


int main() {
    int test;
    scanf("%d", &test);
    while(test--) {
        scanf("%s%d", s, &Q);
        n = strlen(s);
        initalize();
        decode();
    }
}
```

# 【CC STREETTA】The Street
## 题目大意

两个序列$a$ , $b$ ， 三个操作：
1. 对序列$a$的区间加等差数列；
2. 对序列$b$的区间进行等差数列取$max$ ;
3. 查询$a_i+b_i$ . 

## 解题报告

区间加等差数列是很简单的， 因为$ax+b$中的$a,b$分别作为一个标记， 很容易合并（我是直接把两个标记永久化了，听说常数很好？）

区间等差数列求$max$?是李超线段树的裸题， 其实就是加入一个线段， 线段树的每个节点保留优势$a,b$(就是对$mid$位置占优的$a,b$)，其他的下传， 同样标记永久化。

## 代码

win 切 linux 的时候忘同步git了，所以暂时找不到了...

# 【CC LYRC】Music & Lyrics
## 题目大意

对于每个询问串， 查询匹配串中一共出现了多少次。

## 解题报告

这个是裸题， 对于询问串形成的trie建立一个AC自动机， 然后每个匹配串在AC自动机上匹配， 对于每个前缀匹配的位置打一个$+1$的标记。

在fail树上将标记进行子树和，每个询问串出现的次数就是对应节点的子树和。

时间复杂度$O(52*length)$


## 代码

```c++
#include <bits/stdc++.h> 
 
using namespace std;
#define int ll 
#define REP(i,a,b) for(int i=int(a),nn=int(b);i<=nn;++i)
#define VEP(i,a,b) for(int i=int(a),nn=int(b);i>=nn;--i) 
#define rep(i,b) for(int i=0,nn=int(b);i<b;++i)
 
typedef long long ll;
typedef double ff; 
 
char s[50010]; 
int w, n, cnt=1, son[3000000][52], las[510], fal[3000000]; 
ll v[3000000];
 
inline int insert(char *s) { 
	int x = 1, tp; 
	rep(i, strlen(s)) { 
		if (s[i] >= 'a') tp=s[i]-'a'; 
		else tp = s[i]-'A'+26; 
		if (!son[x][tp]) son[x][tp]=++cnt; 
		x = son[x][tp]; 
	} 
	return x; 
} 
 
inline void AC_auto() {
	queue<int> q; q.push(1); 
	int x, y; 
	while (!q.empty()) { 
		x = q.front(); q.pop(); 
		rep(i, 52) { 
			int fa = fal[x]; 
			while (fa && !son[fa][i]) fa = fal[fa]; 
			if (y=son[x][i]) 
				fal[y] = (son[fa][i]?son[fa][i]:1), q.push(y);
			else son[x][i] = (son[fa][i]?son[fa][i]:1); 
		} 
	}
} 
	
inline void mark(char *s) { 
	int x = 1, tp, y; 
	rep(i, strlen(s)) 
		if (s[i] == '-') x = 1; 
		else { 
			if (s[i] >= 'a') y=s[i]-'a'; 
			else y = s[i]-'A'+26; 
			while (x&&!son[x][y]) x=fal[x]; 
			x = son[x][y]; 
			if (!x) x = 1;
			v[x] ++; 
		} 
}
		
main() { 
	scanf("%d", &w); 
	rep(i, w) scanf("%s", s), las[i]=insert(s); 
	AC_auto(); 
	scanf("%d", &n); 
	rep(i, n) scanf("%s", s), mark(s); 
	VEP(i, cnt, 2) v[fal[i]] += v[i]; 
	rep(i, w) printf("%lld\n", v[las[i]]); 
	return 0; 
}
```

# 【CC QTREE6】Query on a tree VI
## 题目大意

开始一棵树全都是黑色的， 每次两个操作：
1. 将一个点变色；
2. 询问和一个点颜色相同的联通块大小；

## 解题报告

这个题是Qtree系列比较简单的?

对黑色和白色分别维护子树中与他单色联通的点的数量， 这个使用树链剖分+树状数组就好了， 每次修改相当于修改到根路径， 区间修改单点查询。

还需要得到与某个点单色联通的深度最前的祖先， 这个需要使用树链剖分+维护区间黑色数（树状数组就足够），

每次如果一条重链单色联通就像上跳， 否则二分单色联通的最靠上的位置。 

时间复杂度？ $O(n \log^2n)$

## 代码

```c++
#include <bits/stdc++.h> 
 
using namespace std; 
 
#define REP(i,a,b) for(int i=int(a),nn=int(b);i<=nn;++i)
#define VEP(i,a,b) for(int i=int(a),nn=int(b);i>=nn;--i)
#define rep(i,b) for(int i=0,nn=int(b);i<nn;++i) 
 
const int N = 100100; 
 
struct edge {
	int nxt, to; 
	edge(int nxt=0, int to=0):nxt(nxt), to(to){}
} e[N << 1]; 
int s0[N], s1[N], oo[N], n, m, head[N], tot, a[N]; 
int f[N], de[N], dn[N], tp[N], sn[N], bd[N], sz[N], cnt; 
 
inline void in(int &x) { 
	char ch = getchar(); int f=1; 
	for(;ch<'0'||ch>'9';ch=getchar())
		if (ch=='-') f = -1; 
	for(x=0; ch>='0'&&ch<='9';ch=getchar())
		x = x*10 + ch - 48; 
	x *= f; 
} 
 
inline void link(int x, int y) {
	e[++tot]=edge(head[x], y), head[x]=tot; 
} 
 
void dfs(int x, int y=0) { 
	de[x] = de[f[x]]+1, sz[x]=1, sn[x]=0;  
	for (int i=head[x];i;i=e[i].nxt) 
		if (y=e[i].to, y!=f[x]) {
			f[y]=x, dfs(y), sz[x] += sz[y]; 
			if (sz[y]>sz[sn[x]]) sn[x]=y; 
		} 
}
	
void odr(int x, int top) { 
	tp[x]=top, dn[x]=++cnt, bd[cnt]=x; int y; 
	if (sn[x]) {
		odr(sn[x], top); 
		for (int i=head[x];i;i=e[i].nxt) 
			if (y=e[i].to, y!=sn[x]&&y!=f[x])
				odr(y, y); 
	} 
}
inline void ad0(int x, int v) {
	for (; x <= n; x += x&-x) s0[x]+=v; 
}
inline void ad1(int x, int v) { 
	for (; x <= n; x += x&-x) s1[x]+=v; 
} 
inline void add(int x, int v) { 
	for (; x <= n; x += x&-x) oo[x]+=v;
}
inline int ak1(int x, int sm=0) {
	for (; x; x-=x&-x) sm+=s1[x]; return sm; 
} 
inline int ak0(int x, int sm=0) {
	for (; x; x-=x&-x) sm+=s0[x]; return sm; 
}
inline int ask(int x, int sm=0) {
	for (; x; x-=x&-x) sm+=oo[x]; return sm; 
} 
 
inline void chng(int x, int y, int v, int t) { 
//	cout << x << ' ' << y <<  ' ' << v << ' ' << t << endl;
	if (de[x] < de[y]) return; 
	for (; tp[x]!=tp[y]; x=f[tp[x]])
		(t) ? (ad1(dn[tp[x]],v), ad1(dn[x]+1,-v))
		  : (ad0(dn[tp[x]],v), ad0(dn[x]+1,-v));
	t ? (ad1(dn[y], v), ad1(dn[x]+1,-v))
	  : (ad0(dn[y], v), ad0(dn[x]+1,-v)); 
} 
 
inline int dd(int l, int r, int t) {
	int ans=r, md;
	while (l<=r) { md = (l+r) >> 1; 
		if (t) {
			if (ask(r) - ask(md-1) == r-md+1) 
			 	ans = md, r = md-1; 
			else l = md + 1; 
		} else {
			if (ask(r) - ask(md-1)) l = md+1; 
			else ans =md, r = md-1;  
		}
	}
	return bd[ans]; 
}
			
inline int up(int x, int t) { 
	while (tp[x] != 1) if (t) 
		if (dn[x]-dn[tp[x]]+1==ask(dn[x])-ask(dn[tp[x]]-1))
			if (a[f[tp[x]]]) x=f[tp[x]]; else return tp[x]; 
		else return dd(dn[tp[x]], dn[x], t); 
	else 
		if (ask(dn[x])-ask(dn[tp[x]]-1)) return dd(dn[tp[x]],dn[x],t); 
		else if (a[f[tp[x]]]) return tp[x]; else x=f[tp[x]]; 
	return dd(1, dn[x], t); 
}
int main() { 
//	freopen("A.in", "r", stdin); 
//	freopen("A.out", "w", stdout); 
	in(n); int x, y; 
	rep(i, n-1) in(x), in(y),link(x,y),link(y,x); 
	dfs(1); odr(1, 1); 
	REP(i, 1, n) {
		ad1(dn[i], sz[i]), ad1(dn[i]+1,-sz[i]);
		add(dn[i], 1), a[i] = 1; 
	} 
	ad0(1, 1), in(m); int t; 
	while (m--) { 
		in(t), in(x); 
		if (t) {
			a[x] ? (chng(f[x], max(1, f[up(x,a[x])]), -ak1(dn[x]), 1),
			add(dn[x], -1), chng(f[x], max(1, f[up(x,a[x]=0)]),ak0(dn[x]), 0))
		:	(chng(f[x], max(1, f[up(x,a[x])]),-ak0(dn[x]), 0),
			add(dn[x], 1), chng(f[x], max(1, f[up(x,a[x]=1)]),ak1(dn[x]),1)); 
		} else {
			printf("%d\n", (a[x]?ak1(dn[up(x,a[x])]):ak0(dn[up(x,a[x])])));
		}
	}
	puts(""); 
	return 0; 
} 
```
# 【CC RIN】Course Selection
## 题目大意

每个学期学一个课程会有一个收益， 每个课程有若干个前置课程， 问获得的最大收益。

## 解题报告

这个是鸟哥(@faebdc)论文里的题目。

把每个课程在不同的学期拆点， 然后相邻的两个学期连边， $S$连到第一个学期， 最后一个学期连到$T$, 流量限制是与最大收益的差值； 

如何对选课的顺序进行限制? 通过连inf边强制前置课程在前选择。

## 代码

```c++
#include <bits/stdc++.h> 
 
using namespace std;
 
#define FORU(i, a, b) for (int i = int(a), nn = int(b); i <= nn; ++i) 
#define FORD(i, a, b) for (int i = int(a), nn = int(b); i >= nn; --i) 
#define REP(i, b) for (int i = 0, nn = int(b); i < b; ++i) 
 
typedef long long ll; 
typedef double ff; 
 
const int N = 110; 
const int inf=1000000000;
 
struct edge { 
	int nxt, to, f; 
	edge(int nxt=0,int to=0,int f=0) 
		:nxt(nxt), to(to), f(f) {}
} e[100000]; 
 
int head[100000], X[N][N], n, m, k, S, T, tot=1; 
 
template <typename T> 
 
inline void cmax(T &x, T a) {
	if (a > x) x = a; 
} 
 
inline int id(int x, int y) { 
	return (x-1) * m + y; 
} 
 
inline void add(int x, int y, int f) {
	e[++tot]=edge(head[x], y, f), head[x]=tot; 
	e[++tot]=edge(head[y], x, 0), head[y]=tot; 
} 
 
int de[100000]; 
queue<int> q; 
 
inline bool bfs() {
	q.push(S); memset(de, 0, sizeof(de)); 
	de[S] = 1; 
	int x, y; 
	while (!q.empty()) { 
		x = q.front(), q.pop(); 
		for (int i=head[x]; i; i=e[i].nxt) 
			if (y = e[i].to, e[i].f&&!de[y])
				de[y] = de[x] + 1, q.push(y); 
	} 
	return de[T]; 
}
 
int dfs(int x, int mx) { 
	if (x == T) return mx; 
	int hv = 0, y, f;
	for (int i=head[x];i;i=e[i].nxt)
		if (y=e[i].to, e[i].f&&de[y]==de[x]+1) { 
			f = dfs(y, min(mx, e[i].f)); 
			e[i].f -= f, e[i^1].f += f, hv += f; 
			mx -= f; if (!mx) return hv; 
		} 
	de[x] = -1; return hv; 
} 
 
int MINCOST() { 
	int rec = 0; 
	while (bfs()) rec += dfs(S, inf); 
	return rec; 
} 
 
int main() { 
	scanf("%d%d%d", &n, &m, &k); 
	S = 0, T = n*m + 1; int ans = 0;  
	
	for (int i = 1; i <= n; ++i) {
		int mx = 0; 
		for (int j=1; j<=m; ++j) 			
			scanf("%d", &X[i][j]), cmax(mx, X[i][j]); 
		for (int j=1; j<m; ++j) 
			add(id(i, j), id(i, j+1), (X[i][j]==-1)?inf:mx-X[i][j]); 
		add(id(i, m), T, (X[i][m]==-1)?inf:mx-X[i][m]); 
		add(S, id(i, 1), inf), ans += mx; 
	} 
	
	for (int i=1; i<=k; ++i) { 
		int A, B; scanf("%d%d", &A, &B); 
		for (int j=1; j<=m-1; ++j) 
			add(id(A, j), id(B, j+1), inf); 
		add(id(A, m), T, inf); 
	} 
	
	ans -= MINCOST(); 
	
	return printf("%.2lf\n", (ff)ans/(ff)n), 0; 
} 
```

# 【CC GNUM】Game of Numbers
## [题目链接](https://www.codechef.com/problems/GNUM)

## 解题报告

（这个人怎么懒到不想写题目大意呀）

把$gcd(i,j)!=1$的点分为$b_j > a_i$ 和 $b_j < a_i$ 两类， 建成连$S$ \ $T$ 的两排点， 把出现过的质因数建成中间的一排点。

$(i,j)$ 连 $gcd(i,j)$ 包含的质数的点， 流量是$1$, 直接最大流。

## 代码

```c++
#include <bits/stdc++.h> 
 
using namespace std;
 
#define FORU(i, a, b) for (int i = int(a), nn = int(b); i <= nn; ++i) 
#define FORD(i, a, b) for (int i = int(a), nn = int(b); i >= nn; --i) 
#define REP(i, b) for (int i = 0, nn = int(b); i < b; ++i) 
 
typedef long long ll; 
typedef double ff; 
 
const int inf = 1000000000; 
const int N = 410; 
 
struct edge { 
	int nxt, to, f; 
	edge(int nxt=0, int to=0, int f=0) 
		:nxt(nxt), to(to), f(f) {}
} e[N*N*10]; 
 
int test, S, T, n, totA, totB, totp, tot;
int a[N], b[N], c[N], head[N*N*10];
 
int A[N*N], Ap[N*N][12];
int B[N*N], Bp[N*N][12];
 
map <int,int> HA, HB, Hp;
 
inline void add(int x, int y, int f) {
	e[++tot]=edge(head[x], y, f), head[x]=tot; 
	e[++tot]=edge(head[y], x, 0), head[y]=tot; 
} 
 
int de[100000]; 
queue<int> q; 
 
inline bool bfs() {
	q.push(S); memset(de, 0, sizeof(de)); 
	de[S] = 1; 
	int x, y; 
	while (!q.empty()) { 
		x = q.front(), q.pop(); 
		for (int i=head[x]; i; i=e[i].nxt) 
			if (y = e[i].to, e[i].f&&!de[y])
				de[y] = de[x] + 1, q.push(y); 
	} 
	return de[T]; 
}
 
int dfs(int x, int mx) { 
	if (x == T) return mx; 
	int hv = 0, y, f;
	for (int i=head[x];i;i=e[i].nxt)
		if (y=e[i].to, e[i].f&&de[y]==de[x]+1) { 
			f = dfs(y, min(mx, e[i].f)); 
			e[i].f -= f, e[i^1].f += f, hv += f; 
			mx -= f; if (!mx) return hv; 
		} 
	de[x] = -1; return hv; 
} 
 
inline int dinic() { 
	int rc = 0; 
	while (bfs()) rc += dfs(S, inf);
	return rc; 
} 
 
int main() {
	scanf("%d", &test);
	
	while (test --) { 
	
		HA.clear(), HB.clear(), Hp.clear();
		totA = totB = totp = 0;
		
		scanf("%d",&n);
		
		int u, v, k, i, j; 
		for(i=1; i<=n; ++i) scanf("%d",&a[i]);
		for(i=1; i<=n; ++i) scanf("%d",&b[i]);
		
		for(i=1; i<=n ; ++i) for(j=1 ; j<=n ; ++j)
			if(a[i] != b[j]) {
				u = __gcd(a[i],b[j]);
				v = 1, c[0] = 0;
				for(k=2; k*k <= u; ++k)
					if(u%k == 0) {
						v *= k, c[++c[0]] = k;
						while (u % k==0) u/=k; 
					}
				if(u!=1) v *= u, c[++c[0]]=u;
				
				if(a[i] > b[j]) {
				
					if(!HA[v]) {
						HA[v] = ++totA;
						A[totA]=1, Ap[totA][0] = c[0];
						for(k=1; k<=c[0]; ++k) {
							if (!Hp[c[k]]) Hp[c[k]] = ++totp;
							Ap[totA][k] = Hp[c[k]];
						}
					} else ++A[HA[v]];
					
				} else {
				
					if(!HB[v]) {
						HB[v] = ++totB;
						B[totB] = 1, Bp[totB][0] = c[0];
						for(k=1; k <= c[0]; ++k) {
							if(!Hp[c[k]]) Hp[c[k]] = ++totp;
							Bp[totB][k] = Hp[c[k]];
						}
					} else ++B[HB[v]];
					
				}
			}
		S=0, T=totA + totp + totB + 1;
		
		memset(e, 0, sizeof(e)); 
		memset(head, 0, sizeof(head)); 
		tot = 1; 
		
		for(i=1; i <= totA; ++i) {
			add(S, i, A[i]);
			for(j=1; j <= Ap[i][0]; ++j)
				add(i, totA+Ap[i][j], inf);
		}
		
		for(i = 1; i <= totB; ++i) {
			for(j=1; j <= Bp[i][0]; ++j)
				add(totA+Bp[i][j], totA+totp+i, inf);
			add(totA+totp+i, T, B[i]);
		}
		
		printf("%d\n",dinic());
	}
	return 0; 
}  
```
# 【CC QTREE2】Counting on a Tree
## 题目大意

每次修改一个点的权值， 求路径$gcd$为$1$的路径条数。

## 解题报告

设$f(i)$表示路径$gcd$是$i$的倍数的路径条数， 答案是$\sum_i \mu(i)f(i) $ .

如果没有修改， 那么每条边会被作为$2^t$个数的倍数， 其中$t$是边权的质因子个数（因为如果幂指数不是$1$，$\mu$就是$0$), 通过并查集得到$f(i)$是十分容易的。

因为修改的次数很少， 所以可以把从来没有修改的边先加入并查集中， 然后对于每一个询问， 把有修改的边修改好加入并查集， 统计好答案后还原。

因为需要还原， 所以并查集需要用启发式合并。

复杂度$O(2^t(n+Q^2)\log n)$

## 代码

```c++
#include<cstdio>
#include<stack>
#include<vector>
#include<algorithm>
using namespace std;
 
#define fi first
#define se second
#define mp make_pair
#define pb push_back
 
typedef long long ll; 
typedef pair<int, int> pii; 
 
const int N = 100005, P = 1000005;
 
int n, Q, f[N], size[N], a[N], c[N];
ll res, ans[N];
int mu[P], p[P], vis[P], t;
bool lock[N];
 
struct edge { int u, v, w; } s[N];
 
vector<int> del;
vector< pii > q[P];
 
inline void init() {
	static int prime[P], tot = 0;
	mu[1] = 1; ++t;
	for(int i = 2; i < P; ++i) {
		if(vis[i] != t) mu[prime[++tot] = p[i] = i] = -1;
		for(int j = 1; j <= tot && prime[j] * i < P; ++j) {
			vis[prime[j] * i] = t;
			p[prime[j] * i] = prime[j];
			if(i % prime[j] == 0) break;
			mu[i * prime[j]] = -mu[i];
		}
	}
}
 
inline void divide(pii op, int val) {
	vector<int> ds; ds.pb(1);
	while(val != 1) {
		int x = p[val];
		for(int i = ds.size(); i--; ) ds.pb(x * ds[i]);
		while(val % x == 0) val /= x;
	}
	for(int i = ds.size(); i--; ) q[ds[i]].pb(op);
}
 
inline ll ask(int i) { 
	return 1LL * i * (i - 1) >> 1; 
}
	 
inline int _find(int x) { 
	return f[x] == x ? x : _find(f[x]); 
}
 
stack< pii > S;
 
inline void _union(int u, int v) {
	u = _find(u); v = _find(v);
	if(u == v) return;
	if(size[u] < size[v]) swap(u, v);
	res -= ask(size[u]) + ask(size[v]);
	size[f[v] = u] += size[v];
	res += ask(size[u]);
	S.push(make_pair(u, v));
}
 
inline void _restore() {
	int u = S.top().fi, v = S.top().se; S.pop();
	res -= ask(size[u]);
	size[u] -= size[f[v] = v];
	res += ask(size[u]) + ask(size[v]);
}
 
inline void solve(int p) {
	vector< pii >& op = q[p];
	if(op.empty()) return;
	
	stack< pii >().swap(S);
	res = 0;
	for(auto x : op) {
		int u = s[x.fi].u, v = s[x.fi].v;
		size[f[u] = u] = size[f[v] = v] = 1;
	}
	
	for(int l = 0, r = 0; l < op.size(); l = r) {
		for(; r < op.size() && op[r].se == op[l].se; ++r)
			_union(s[op[r].fi].u, s[op[r].fi].v);
		if(op[l].se < 0) {
			stack< pii >().swap(S);
			++t;
			for(int i = r; i < op.size(); ++i) vis[op[i].se] = t;
			for(int i = 0; i <= Q; ++i) if(vis[i] != t) ans[i] += res * mu[p];
		}
		else {
			ans[op[l].se] += res * mu[p];
			while(!S.empty()) _restore();
		}
	}
}
 
int main() {
	init();
	scanf("%d", &n);
	for(int i = 1; i < n; ++i) 
		scanf("%d%d%d", &s[i].u, &s[i].v, &s[i].w);
 
	scanf("%d", &Q);
	for(int i = 1; i <= Q; ++i) {
		scanf("%d%d", a + i, c + i);
		if(!lock[a[i]]) {
			del.pb(a[i]);
			lock[a[i]] = true;
		}
	}
	
	for(int i = 1; i < n; ++i)
		if(!lock[i]) divide(make_pair(i, -1), s[i].w);
	
	for(int i = 1; i < n; ++i)
		if(lock[i]) divide(make_pair(i, 0), s[i].w);
	
	for(int i = 1; i <= Q; ++i) {
		s[a[i]].w = c[i];
		for(int j : del) divide(make_pair(j, i), s[j].w);
	}
	
	for(int i = 2; i < P; ++i)
		if(mu[i] != 0) solve(i);
	
	for(int i = 0; i <= Q; ++i) 
		printf("%lld\n", ans[i] + ask(n));
		
	return 0;
}
```


