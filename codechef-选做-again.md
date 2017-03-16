---
title: codechef-选做-again
date: 2017-03-14 11:28:39
tags:
  - 分块
  - 倍增
  - DP
  - 数论
  - 矩阵乘法
  - 匹配
  - 线段树
  - 树链剖分
  - 组合数学
  - BSGS
categories:
  - 题目集锦
---

> 继续做CC中难度比较低的题目

<!--more-->

# 【CC TRIPS】Children Trips
## 题目大意

一棵树上每个边的权值为1\2， 每次给出一个小盆友的体力值， 每天可以行进的距离$\leqslant$体力值，问走$(x,y)$这条路径需要几天。

## 解题报告

首先考虑体力值为1 、2 、3...的分别走步数$n$ , 总共需要走的距离是$n^3$级别的，实际总步数是$n^2$级别， 所以每个体力值走过的步数是$\sqrt{n}$级别， 每一步走到哪一个节点通过倍增（二分？）$\log n$可以得到, 复杂度$O(n\sqrt{n}\log n)$

如果一个体力值走的步数$> n$，那么可以对这个体力值进行重建图， 每个点连向上能到达的最远的点， 在呈现做倍增， 这个重构图的过程最多进行$\sqrt{n}$次， 单次复杂度$O(n\log n)$.

重构图之后， 通过倍增就可以$\log n$解决一次询问， 复杂度$O(n\log n)$.

所以总体的复杂度就是$O(n \sqrt{n} \log n)$

## 代码

```c++
#include <bits/stdc++.h> 
 
using namespace std;
 
#define FORU(i,  a,  b) for (int i = int(a),  nn = int(b); i <= nn; ++i) 
#define FORD(i,  a,  b) for (int i = int(a),  nn = int(b); i >= nn; --i) 
#define REP(i,  b) for (int i = 0,  nn = int(b); i < b; ++i) 
 
typedef long long ll; 
typedef double ff; 
 
const int N = 100010; 
 
struct node {
	int u, v, c, id;
	node(int u=0, int v=0, int c=0, int id=0) 
		:u(u), v(v), c(c), id(id) {}
} t[N];
struct edge { 
	int nxt, to, v; 
	edge(int nxt=0,int to=0,int v=0) 
		:nxt(nxt), to(to), v(v) {}
} e[N*2]; 
 
int n, m, d, lca, sum, c;
int head[N], tot; 
int dep[N], deep[N], fa[N][20], _fa[N][20];
int st[N], top, top2;
int ans[N], s[N*2];
 
inline bool cmp(node a, node b){return a.c>b.c;}
 
inline void add(int x, int y, int v) {
 	e[++tot] = edge(head[x], y, v), head[x]=tot;
} 
 
void dfs(int x) {
	if (deep[x] > d) d = deep[x];
	int y; 
	
	for (int i=0; fa[x][i]; ++i) 
		fa[x][i+1] = fa[ fa[x][i] ][i]; 
 
	for (int i = head[x]; i; i = e[i].nxt)
		if (y = e[i].to, y != fa[x][0]) {
			deep[y] = deep[x] + e[i].v;
			dep[y] = dep[x] + 1, fa[y][0]=x; 
			dfs(y);
		}
}
 
void _dfs(int x) {
	int _top = top2, y;
	for ( st[ ++top ] = x; deep[x]-deep[st[top2]] > c; ++top2);
	_fa[x][0] = st[top2];
	for (int i=0; _fa[x][i]; ++i)
		_fa[x][i+1] = _fa[_fa[x][i]][i];
		
	for (int i=head[x]; i; i = e[i].nxt)
		if (y = e[i].to, y != fa[x][0]) 
			_dfs(y);
			
	--top, top2 = _top;
}
 
int LCA(int u, int v) {
	int i, j;
	if(dep[u] < dep[v]) swap(u, v);
	for (i=0, j=dep[u]-dep[v]; j; ++i, j>>=1)	
		if(j & 1) u = fa[u][i];
	if (u == v) return u;
	for (i=17; i>=0; --i) if(fa[u][i] != fa[v][i])
		 u = fa[u][i], v = fa[v][i];
	return fa[u][0];
}
 
int main() {
	scanf("%d", &n); int u, v, i, j, k; 
	for (i=1; i<n; ++i) {
		scanf("%d%d%d", &u, &v, &c);
		add(u, v, c), add(v, u, c); 
	}
	dep[1] = deep[1]=1, dfs(1);
	 
	scanf("%d", &m);
	for (i=1; i<=m; ++i){
		scanf("%d%d%d", &u, &v, &c); 
		t[i] = node(u, v, c, i); 
		
		s[ t[i].c ] += n / t[i].c;
		if(s[ t[i].c ]>n) s[t[i].c] = n;
	}
	
	sort(t+1, t+m+1, cmp);
	for (i=1, c=-1; i <= m; ++i) {
		u = t[i].u, v = t[i].v;
		lca = LCA(u, v), sum=0;
		if(s[ t[i].c ] < n) {
			c = t[i].c;
			for (; deep[u]-deep[lca]>=c; u=k, ++sum) {
				for (j=17, k=u; j>=0; --j)
					if(deep[u] - deep[fa[k][j]] <= c)
						k = fa[k][j];
				if(k == lca)break;
			}
			
			for (; deep[v]-deep[lca]>=c; v=k, ++sum) {
				for (j=17, k=v; j>=0; --j)
					if(deep[v]-deep[fa[k][j]] <= c)
						k = fa[k][j];
				if(k == lca)break;
			}
			
			if(u!=lca || v!=lca)
				sum += 1 + (deep[u]+deep[v]-2*deep[lca]>c);
		} else {
			if (t[i].c != c) c = t[i].c, _dfs(1);
			for (j=17; j>=0; --j)
				if (deep[_fa[u][j]] - deep[lca]>0)
					u=_fa[u][j], sum += 1<<j;
			for (j=17; j>=0; --j)
				if (deep[_fa[v][j]] - deep[lca]>0)
					v=_fa[v][j], sum += 1<<j;
			if(u!=lca || v!=lca)
				sum += 1+(deep[u]+deep[v]-2*deep[lca]>c);
		}
		ans[t[i].id] = sum;
	}
	for (i=1; i<=m; ++i) 
		printf("%d\n", ans[i]);
	return 0; 
} 
```

# 【CC SEAEQ】Sereja and Equality
## 题目大意

子串相等在题目中的定义就是离散后对应位置相同的意思， 求所有长度为$n$的排列中， 相等的逆序对个数不超过$E$的子串对数。

## 解题报告

首先求出长度为$i$的逆序对个数不超过$j$的排列个数， 设$f[i][j]$表示长度为$i$逆序对个数为$j$的方案数， $$f[i][j]=\sum_{k=0}^{i-1} f[i-1][j-k]$$ $$g[i][j]=\sum_{k=0}^{j}f[i][k]$$
直接转移就可以得到， 复杂度$O(n^2)$

对于一个询问$n,E$ , 枚举相等子串的长度$l$ , 对答案的贡献是$$((n-l)!C_n^l(n-l+1))^2 \times g[l][E]$$

复杂度$O(n^2+nq)$

## 代码

```c++
#include <bits/stdc++.h> 
 
using namespace std;
 
#define FORU(i, a, b) for (int i = int(a), nn = int(b); i <= nn; ++i) 
#define FORD(i, a, b) for (int i = int(a), nn = int(b); i >= nn; --i) 
#define REP(i, b) for (int i = 0, nn = int(b); i < b; ++i) 
#define int ll 
typedef long long ll; 
typedef double ff; 
 
const int N = 510; 
const int p = 1000000007; 
 
int f[N][N*N/2], test, C[N][N], fac[N];
 
main() { 
	f[0][0] = 1, f[1][0] = 1; 
	
	int n = 500, E, pE, dl; 
	for (int i = 2; i <= n; ++i) { 
		E = i*(i-1) / 2, pE = (i-1)*(i-2)/2; 
		for (int j = 0; j <= E; ++j) {
			f[i][j] = f[i-1][min(j, pE)]; 
			dl = j - i; 
			if (dl >= 0) f[i][j] -= f[i-1][dl]; 
			if (f[i][j] >= p) f[i][j] -= p; 
			if (f[i][j] < 0) f[i][j] += p; 
		} 
		for (int j = 1; j <= E; ++j) {
			f[i][j] += f[i][j-1]; 
			if (f[i][j] >= p) f[i][j] -= p ;
		} 
	}
 
	
	C[0][0] = 1; 
	for (int i=1; i<=n; ++i) {
		C[i][0] = 1; 
		for (int j=1; j<=i; ++j) {
			C[i][j] = C[i-1][j] + C[i-1][j-1]; 
			if (C[i][j] >= p) C[i][j] -= p; 
		}
	}

	fac[0] = 1; 
	for (int i=1; i<=n; ++i) 
		fac[i] = 1LL * fac[i-1] * i % p; 
	
	scanf("%d", &test);  
	 
	while (test --) { 
		scanf("%d%d", &n, &E); 
		int ans = 0; 
		for (int l = 1; l <= n; ++l) {
			ll x = C[n][l]  * fac[n-l] % p  ; 
			ans += x * x % p * (n-l+1) % p * f[l][min(l*(l-1)/2, E)] % p;
			if (ans >= p) ans -= p; 
		} 
		printf("%d\n", ans); 
	} 
	
	return 0; 
} 
```

# 【CC PARSIN】Sine Partition Function
## 题目大意

$$f(n, m, x) = \sum_k sin(k_1x)sin(k_2x)···sin(k_mx)$$

其中$k_1+k_2+...+k_m=n$ .

求$f(n, m, x)$


## 解题报告

考虑$n+1$ , 有两个情况， 一种在末尾添加一个$sin(x)$， 一种是将$k_m+1$。

后一种情况出现和和差角问题， 利用高中（我并没有学过）的三角函数知识， 需要利用$cos(k_mx)$进行转移。

令$g(n,m,x)=\sum_k sin(k_1x)sin(k_2x)···cos(k_mx)$

可以得到：
$$f(n,m,x)=f(n-1,m,x)cosx+(f(n-1,m-1,x)+g(n-1,m-1,x))sinx$$
$$g(n,m,x)=(f(n-1,m-1,x)-f(n-1,m,x))sinx+g(n-1,m,x)cosx$$

因为$m$很小， 所以把两个数组都压到矩阵里， 大力快速幂+矩阵乘法就搞了。

## 代码

```c++
#include <bits/stdc++.h> 
 
using namespace std;
 
#define FORU(i, a, b) for (int i = int(a), nn = int(b); i <= nn; ++i) 
#define FORD(i, a, b) for (int i = int(a), nn = int(b); i >= nn; --i) 
#define REP(i, b) for (int i = 0, nn = int(b); i < b; ++i) 
 
typedef long long ll; 
typedef long double ff; 
 
const int M = 31; 
 
struct matrix { 
	int x, y; 
	ff a[M+M][M+M]; 
	matrix() {
		memset(a, 0, sizeof(a)); 
		x = y = 0; 
	} 
	matrix(int x, int y,int t=0) 
	:x(x), y(y) { 
		memset(a, 0, sizeof(a)); 
		if (t == 1) 
			REP(i, x) a[i][i] = 1; 
	} 
	matrix operator *(const matrix b)const{ 
		matrix c(x, b.y) ; 
		REP(i, c.x) REP(j, c.y) REP(k, y) 
			c.a[i][j] += a[i][k]*b.a[k][j]; 
		return c; 
	}
	void print() {
		cout << endl << endl;
		
	 	cout << x << ' ' << y << endl; 
	 	REP(i, x) REP(j, y) cout << a[i][j] << ((j+1==y)?'\n':' ');
	 	
	 	cout << endl; 
	 } 
	  
} trans, ini; 
 
int test, n, m; ff X;
 
matrix fast(matrix x, int k) {
 	matrix ans(x.x, x.y, 1); 
 	for (; k; k>>=1, x=x*x) 
 		if (k & 1) ans = ans*x; 
 	return ans; 
} 
 
int main() {
	ios :: sync_with_stdio(false); 
	
	cin >> test; 
	while (test --) { 
		cin >> m >> n >> X; 
		
		ini = trans = matrix(); 
		ini.x = (m+1)*2, ini.y = 1; 
		ini.a[0][0]=ini.a[m+1][0]=1; 
		
		trans.x=trans.y=(m+1)*2; 
		ff _s = sin(X), _c = cos(X); 
		REP(i, trans.x) 
			if (i == 0 || i == m+1) 
				continue; 
			else if (i <= m) { 
					trans.a[i][i-1]=_s; 
					trans.a[i][i]=_c; 
					trans.a[i][m+1+i]=_s; 
				} else { 
					trans.a[i][i-(m+1)-1]=_c; 
					trans.a[i][i-(m+1)]=-_s; 
					trans.a[i][i]=_c; 
				} 
		
		trans = fast(trans, n); 
		
//		trans.print(); 
		trans = trans * ini; 
		
		cout << trans.a[m][0] << endl; 
	} 
	return 0;
} 
```

# 【CC LECOINS】Little Elephant and Colored Coins
## 题目大意

有$n$种硬币， 每种硬币无限使用， 有一个价值$v$和颜色$c$， 每次询问能否凑出价值$S$， 最多用多少颜色。 

## 解题报告

对价值最小的硬币做模意义， 因为价值最大是$200000$。

设$f[i][j]$表示使用$i$种颜色，价值模意义下为$j$的最小能凑出的价值， 把同一种颜色的硬币放在一起， 对$f[i][j]$从$f[i-1][k]$进行更新， 然后利用取模成环的性质，对$f[i][j]$从$f[i][k]$更新， 复杂度是$O(n^2*V)$ 

（要把和价值最小硬币颜色相同的硬币单独拿出来， 更新$f[i][j]$得到$g[i][j]$）

询问的时候， 枚举$i$, 看$f[i][S%V]$和$S$的大小关系， 在判断$g[i][S%V]$和$S$的大小关系， 更新答案

时间复杂度$O(n^2V+Q*n)$

## 代码

```c++
#include <bits/stdc++.h> 
 
using namespace std;
 
#define FORU(i, a, b) for (int i = int(a), nn = int(b); i <= nn; ++i) 
#define FORD(i, a, b) for (int i = int(a), nn = int(b); i >= nn; --i) 
#define REP(i, b) for (int i = 0, nn = int(b); i < b; ++i) 
#define pb push_back
#define se second 
#define fi first
 
typedef long long ll; 
typedef double ff; 
typedef vector<int> vi; 
typedef pair<int, int> pii; 
 
const ll INF = 1000000000000000000LL;
const int N = 200005;
 
int n, m, cnt; ll mV, mC, V[N]; 
pii A[64];
ll R[64][N],  _R[64][N];
bool U[N];
 
inline void update(vector< pii > T) {
	FORD (i, cnt, 0) REP (it, T.size())
		for (int j = mV-1; j >= 0; j--)
			if (R[i][j] < INF) R[i+1][(j+T[it].se)%mV] = 
					min(R[i+1][(j+T[it].se)%mV], R[i][j]+T[it].se);
					
	REP(i, cnt+2) REP(j, T.size()) {
		ll v = T[j].se, vm = v%mV;
		if (vm == 0) continue;
		int d = __gcd((ll)vm, mV);
		REP(x, d) {
			int z = x, y = x;
			while (1) {
				if (R[i][y] < R[i][z]) z = y;
				y = (y + vm) % mV;
				if (y == x) break;
			}
			y = z;
			while (1) {
				int w = (y + vm) % mV;
				R[i][w] = min(R[i][w], R[i][y]+v);
				y = w; if (y == z) break;
			}
		}
	}
}
 
int main() {
	REP(i, 64) REP(j, N) R[i][j] = INF;
	R[0][0] = 0;
	cin >> n;
	mV = INF, mC = -1;
	REP(i, n) {
		cin >> A[i].se >> A[i].fi;
		if (A[i].se < mV) 
			mV = A[i].se, mC = A[i].fi;
	}
	
	sort(A, A+n);
	int i = 0;
	vector< pii > _T;
	while (i < n) {
		if (A[i].first == mC) {
			if (A[i].se != mV) _T.pb(A[i]);
			i ++; continue;
		}
		int j = 0;
		vector< pii > T;
		while (i+j<n && A[i+j].fi==A[i].fi) 
			T.pb(A[i+j]), j ++;
			
		update(T), i += j, cnt ++;
	}
	
	REP(i, 64) REP(j, mV) _R[i][j] = R[i][j];
	update(_T), cnt ++;
 
	cin >> m;
	
	REP(i, m) { 
		ll s; int res = -1;
		scanf("%lld", &s);
		REP(j, n+1) {
			ll d = _R[j][s % mV];
			if (d >= INF) continue;
			if (d > s) continue;
			if (d == s) res = max(res, j);
			else res = max(res, j + 1);
		}
		
		if (!_T.empty()) {
			REP(j, n+1)  {
				ll d = R[j][s % mV];
				if (d >= INF) continue;
				if (d > s) continue;
				res = max(res, j);
			}
		}
		printf("%d\n", res);
	}
	return 0;
} 
```
#【CC ANUDTQ】 Dynamic Trees and Queries
## 题目大意

四种操作： 
1. 向一个点连一个权值为$v$的叶子；
2. 向一个点的子树加一个权值$v$;
3. 删除一棵子树；
4. 询问一棵子树的权值和；

## 解题报告

本来增加连边的问题要先试试LCT, 但设计子树加和子树查询要维护若干信息， 比较麻烦。

如果用SPLAY维护出栈入栈序就十分的方便， 具体的， 就是每个点有一个入栈点和一个出栈点， 出栈点是没有权值的size的。

剩下就是模板了。。

## 代码

```c++
#include <bits/stdc++.h> 
 
#define REP(i,a,b) for(int i=int(a),nn=int(b);i<=nn;i++)
#define pb push_back
#define mp make_pair
#define fi first
#define se second
#define DEBUG(x) cout << (#x) << ": " << x <<endl;
#define int ll
using namespace std;
 
typedef long long ll;
typedef double f;
typedef pair<int,int> pii;
typedef vector<int> vi;
 
const int N=400020;
 
struct edge{
	int to, nxt; 
	edge(int nxt=0,int to=0):nxt(nxt),to(to){}
} e[N];
int head[N], n, cnt, tot, rec, rt; 
int sn[N][2], fa[N], sz[N], s[N], key[N], sm[N], ps[N], dn[N], ot[N], Eu[N]; 
ll ans;
 
inline int in() { 
	char ch = getchar(); int f=1,x =0; 
	for (;ch<'0'||ch>'9';ch=getchar())
		if (ch=='-') f = -1; 
	for (;ch>='0'&&ch<='9';ch=getchar())
		x=x*10+ch-48; 
	return x*f; 
} 
 
void add(int x,int y){
	e[++tot]=edge(head[x], y), head[x]=tot;
	e[++tot]=edge(head[y], x), head[y]=tot;
}
 
inline int birth(int _v, int _s) {
	++cnt, sn[cnt][0]=sn[cnt][1]=fa[cnt]=0, sz[cnt]=s[cnt]=_s; 
	key[cnt] = sm[cnt] = _v; return cnt;  
}
 
inline void ad(int x, int _v) {
	if (!x) return;  
	ps[x] += _v, sm[x] += _v * sz[x]; 
	key[x] += s[x] * _v; 
} 
	
inline void up(int x) { 
	sz[x] = s[x], sm[x] = key[x]; 
	if (sn[x][0]) sz[x] += sz[sn[x][0]], sm[x] += sm[sn[x][0]]; 
	if (sn[x][1]) sz[x] += sz[sn[x][1]], sm[x] += sm[sn[x][1]]; 
} 
 
inline void down(int x) { 
	if (!x) return; 
	if (ps[x]) ad(sn[x][0], ps[x]), ad(sn[x][1], ps[x]), ps[x]=0; 
} 
 
inline void rotate(int x) { 
	int y = fa[x], z = fa[y], d = (sn[y][1]==x); 
	fa[x] = z; if (z) sn[z][sn[z][1]==y] = x; 
	if (sn[x][d^1]) fa[sn[x][d^1]] = y; sn[y][d] = sn[x][d^1]; 
	fa[y] = x, sn[x][d^1] = y, up(y); 
} 
 
inline void splay(int x, int aim=0){
	static int stk[N], top=0, tmp; tmp = x;  
	while (tmp != aim) stk[++top]=tmp, tmp=fa[tmp]; 
	while (top) down(stk[top]), --top; 
	for (int y=fa[x]; y!=aim; rotate(x), y=fa[x]) {
		if (fa[y] == aim) continue; 
		if ((sn[y][0]==x)^(sn[fa[y]][0]==y)) rotate(x); 
		else rotate(y); 
	}
	up(x); if (!aim) rt = x; 
} 
 
void dfs(int x,int fa){
	Eu[ ++rec ] = dn[x];
	for (int i=head[x], y; i; i=e[i].nxt)
		if (y = e[i].to, y != fa) dfs(y, x);
	Eu[ ++rec ] = ot[x];
}
 
int build(int l, int r){
	if (l > r) return 0; 
	int mid = (l+r) >> 1, x = Eu[mid];
	sn[x][0] = build(l, mid-1);
	if (sn[x][0]) fa[sn[x][0]] = x;
	sn[x][1] = build(mid+1, r); 
	if (sn[x][1]) fa[sn[x][1]] = x; 
	up(x); return x;
}
int pre(int x){
	splay(x);
	if (!sn[x][0]) return 0;
	for (x = sn[x][0], down(x); sn[x][1];) 
		x = sn[x][1], down(sn[x][1]);
	return x;
}
 
main() {
//	freopen("A.in", "r", stdin); 
	n = in();
	REP(i, 1, n) {
		dn[i] = birth(in(), 1);
		ot[i] = birth(0, 0);
	}
	REP(i, 1, n-1) add(in()+1, in()+1);
	dfs(1, 0), rt = build(1, rec);
	
	int t, x, tmp; 
	for (int m = in(); m--; ) {
		t = in(), x = in()+1+ans;
		if (t == 1){
			dn[++n] = birth(in(), 1), ot[n] = birth(0, 0);
			sn[dn[n]][1] = ot[n], fa[ot[n]] = dn[n]; 
			splay(ot[x]), tmp = sn[ot[x]][0];
			sn[ot[x]][0] = dn[n], fa[dn[n]] = ot[x];
			sn[dn[n]][0] = tmp;
			if (tmp) fa[tmp] = dn[n];
			up(dn[n]), up(ot[x]); 
//			DEBUG(sz[ot[x]]); 
		}  	
		if (t == 2) {
			if (!pre(dn[x])) 
				splay(ot[x]), ad(sn[rt][0], in()), up(rt); 
			else {
				splay(pre(dn[x])), splay(ot[x], rt);
				ad(sn[sn[rt][1]][0], in());
				up(sn[rt][1]), up(rt); 
			}
		}
		if (t == 3) {
			splay(pre(dn[x])), splay(ot[x], rt);
			sn[sn[rt][1]][0] = 0;
			up(sn[rt][1]), up(rt); 
		}
		if (t == 4) {
			if (!pre(dn[x]))
				splay(ot[x]), printf("%lld\n",ans=sm[sn[rt][0]]);
			else {
				splay(pre(dn[x])), splay(ot[x], rt);
				printf("%lld\n", ans=sm[sn[sn[rt][1]][0]]);
			}
		}
	}
	return 0; 
}  
```

# 【CC MATCH】Expected Maximum Matching
## 题目大意

一个二分图， 点$i$和点$j$之间有边相连的概率为$p[i][j]$, 求期望最大匹配数。 

## 解题报告

有一个叫做Hall定理的东西， 内容大概是， 二分图存在完美匹配的充分必要条件是一类点的任意一个大小为$x$的点集都与另一类点中$>=x$个点相连。 

可以注意到， 左边的点数量少到出奇， 这使得合法的状态非常的少， 使用$2^{2^n}$表示子集是否满足Hall定理，因为子集满足Hall定理具有包含关系，也就是一个集合满足， 他的子集都需要满足， 所以通过打表， 发现状态数十分的有限。 

然后预处理转移， 对于一个状态$x$, 预处理出与一个右侧的点联通情况为$S$时，得到的Hall定理状态为$nxt$ ，再预处理出与左侧每个点联通情况为$S$的概率。 

之后直接dp, $f[i][j]$表示右侧前$i$个点， Hall定理状态为$j$的概率， 直接计算期望就好了。 


## 代码

```c++
#include <bits/stdc++.h> 

using namespace std; 

#define REP(i,a,b) for(int i=int(a),nn=int(b);i<=nn;++i)
#define VEP(i,a,b) for(int i=int(a),nn=int(b);i>=nn;--i) 
#define rep(i,b) for(int i=0,nn=int(b);i<nn;++i)
#define B(x) (1u<<(x))
typedef long long ll; 
typedef double ff; 
typedef unsigned int ui; 

const int N = 5010; 
const int M = 105; 

int trans[N][40],ful[6],as[N],n,m,w,tot; 
ff f[M][N],g[M][40],p[M][6]; ui q[N]; map<ui,int> hsh; 
inline int find(ui S) {
	if (hsh.count(S) == 0) q[hsh[S] = tot++]=S; 
	return hsh[S]; 
}

int main() { 
	scanf("%d%d", &n, &m); w = 1<<n; 
	rep(i, n) rep(j, m) scanf("%lf", &p[j][i]); 
	rep(i,w) {int cnt=0; rep(j,n) if(i&B(j)) ++cnt; ful[cnt]|=B(i);}
	int k = 0; 
	for (find(1); k<tot; ++k) { 
		ui x = q[k]; 
		static ui nxt[6];
		rep(i, n) nxt[i] = 0; 
		rep(i, w) if (x & B(i)) rep(j, n) nxt[j]|=B(i|B(j)); 
		rep(i, w) { ui S = x ; 
			rep(j, n) if (i&B(j)) S|=nxt[j]; 
			trans[k][i] = find(S); 
		} 
		VEP(i, n, 0) if (x&ful[i]) {as[k]=i; break;}
	}
	rep(i, m) rep(x, w) { g[i][x]=1;rep(j,n) 
		if(x&B(j))g[i][x]*=p[i][j];else g[i][x]*=(1-p[i][j]);} 
	
	
	f[0][0]=1; rep(i, m) rep(k, tot) rep(x, w) 
		f[i+1][trans[k][x]] += f[i][k] * g[i][x]; 

	ff ans = 0; 
	rep(i, tot) ans += f[m][i] * (ff)as[i]; 
	
	printf("%.10lf\n", ans); 
	return 0; 
} 	
```

# 【CC COT5】Count on a Treap
## 题目大意

维护一颗Treap, 支持:
1. 插入一个key和weight分别为k, w的点； 
2. 删除一个key是k的点；
3. 返回key是ku和kv的两个点在树上的距离； 

## 解题报告

主要就是挖掘Treap的性质， 也就是Treap上的点$x$是点$y$的祖先， 当且仅当key在$x$，$y$之间的数重量都小于$w_x$; 

treap上两个点的LCA就是key在两点之间的重量最大数， 这个使用区间最大值的 位置就可以确定；

那么关键就是确定一个点在treap中的深度， 这个利用点$x$是点$y$的祖先的条件， 分别寻找点左右两侧的祖先数量， 求和就可以得到深度。

点$x$左、右祖先的个数， 实际是维护一个向左、右单调递增的序列，这个通过记录区间max和区间向左向右单增序列长度， 可以单次$\log n$(分入左右子树)进行合并， 从而$O(n\log^2 n)$修改+查询； 

## 代码

```c++
#include <bits/stdc++.h> 
 
using namespace std;
 
#define REP(i,a,b) for(int i=int(a),nn=int(b);i<=nn;++i)
#define VEP(i,a,b) for(int i=int(a),nn=int(b);i>=nn;--i)
#define rep(i,b) for(int i=0,nn=int(b);i<nn;++i)
 
typedef unsigned int ui; 
 
const int N=200010;
 
struct ask  {
	int tp; ui x, y;
	void read() {
		scanf("%d", &tp);
		if (tp==1) scanf("%u",&x); 
		else scanf("%u%u", &x, &y);
	}
} as[N];
 
ui x[N],y[N];
 
int n, m, now, ans;
 
int find(ui *a, ui x) { 
	int l=1, r=m+1, ans=0;
	while (l < r) {
		int mid = l+r>>1;
		if (a[mid] >= x)
			ans = mid, r = mid;
		else l = mid+1;
	}
	return ans;
}
 
int l[N<<2], r[N<<2], we[N<<2], mx[N<<2];
 
int findl(int x, int v) {
	if (v==0) return l[x];
	if (v > mx[x]) return 0;
	if (l[x] == 1) return 1;
	if (v <= mx[x<<1]) 
		return l[x]-l[x<<1]+findl(x<<1, v);
	return findl(x<<1|1, v);
}
 
int findr(int x, int v) {
	if (v==0) return r[x];
	if (v > mx[x]) return 0;
	if (r[x] == 1) return 1;
	if (v <= mx[x<<1|1]) 
		return r[x]-r[x<<1|1]+findr(x<<1|1,v);
	return findr(x<<1,v);
}
 
void push_up(int x) {
	if (mx[x<<1]>mx[x<<1|1]) {
		mx[x]=mx[x<<1]; we[x]=we[x<<1];
	} else  {
		mx[x]=mx[x<<1|1]; we[x]=we[x<<1|1];
	}
	l[x]=l[x<<1]+findl(x<<1|1, mx[x<<1]);
	r[x]=r[x<<1|1]+findr(x<<1,mx[x<<1|1]);
}
 
void insert(int x,int ll,int rr,int _l,int _r) {
	if (ll == rr) {
		l[x]=1, r[x]=1, we[x]=ll, mx[x]=_r;
	} else {
		int mid=ll+rr>>1;
		if (_l<=mid) insert(x<<1,ll,mid,_l,_r); 
		else insert(x<<1|1,mid+1,rr,_l,_r);
		push_up(x);
	}
}
 
void del(int x,int ll,int rr,int p) {
	if (ll==rr) {
		l[x]=0, r[x]=0, we[x]=0, mx[x]=0;
	} else {
		int mid=ll+rr>>1;
		if (mid>=p) del(x<<1,ll,mid,p); 
		else del(x<<1|1,mid+1,rr,p);
		push_up(x);
	}
}
 
int findl(int x,int ll,int rr,int p) {
	if (ll==rr) {
		now=mx[x]; return 1; 
	} else{
		int mid=ll+rr>>1;
		if (mid>=p) {
			int ans=findl(x<<1,ll,mid,p);
			ans+=findl(x<<1|1,now);
			now=max(now,mx[x<<1|1]); return ans;
		} else return findl(x<<1|1,mid+1,rr,p);
	}
}
 
int findr(int x,int ll,int rr,int p) {
	if (ll==rr) {
		now=mx[x]; return 1;
	} else {
		int mid=ll+rr>>1;
		if (mid>=p) return findr(x<<1,ll,mid,p);
		else  {
			int ans=findr(x<<1|1,mid+1,rr,p);
			ans+=findr(x<<1,now);
			now=max(now,mx[x<<1]); return ans;
		}
	}
}
 
int findd(int x) {
	return findl(1,1,m,x)+findr(1,1,m,x);
}
 
void findmax(int x,int ll,int rr,int _l,int _r) {
	if (ll>_r||rr<_l) return;
	if (ll>=_l&&rr<=_r) {
		if (mx[x]>now) 
			now=mx[x], ans=we[x];
		return;
	}
	int mid=ll+rr>>1;
	findmax(x<<1,ll,mid,_l,_r);
	findmax(x<<1|1,mid+1,rr,_l,_r);
}
 
int main() {
	scanf("%d",&n);
	for (int i=1;i<=n;i++) as[i].read();
	m=0;
	for (int i=1;i<=n;i++) if (as[i].tp==0) {
		x[++m]=as[i].x; y[m]=as[i].y;
	}
	
	sort(x+1,x+m+1);
	sort(y+1,y+m+1);
	
	for (int i=1;i<=n;i++) {
		if (as[i].tp==0) {
			int l=find(x,as[i].x),r=find(y,as[i].y);
			insert(1,1,m,l,r);
		} 
		if (as[i].tp==1) {
			int l=find(x,as[i].x); del(1,1,m,l);
		} 
		if (as[i].tp==2) {
			int l=find(x,as[i].x),r=find(x,as[i].y); now=0; ans=0;
			if (l>r) swap(l,r);
			findmax(1,1,m,l,r);
			printf("%d\n",findd(l)+findd(r)-findd(ans)*2);
		}
	}
	return 0;
}  
```

# 【CC FN】Fibonacci Number
## 题目大意

在$\mod \ p$意义下, 斐波那契数列第$n$项$f_n=C$的最小的$n$.

## 解题报告

题目中说$p \mod 10=1/9$, 这个性质就保证了$\sqrt{5}$在$\mod p$意义下存在对应整数（$5^{(p-1)/2}=1(\mod p)$）；

现在令$x=\frac{\sqrt{5}+1}{2}$, 则$x^n-(-x)^{-n}=\sqrt{5}C$

从而， 对$n$分奇偶情况讨论， 得到： 

1. $n$为偶数， 且$(x^n)^2-\sqrt{5}C \times (x^n)-1=0$;
2. $n$为奇数， 且$(x^n)^2-\sqrt{5}C \times (x^n)+1=0$;

利用求根公式+二次剩余+BSGS可以得到最小的$n$ ;

关于二次剩余， 求$x^2=n(\mod p)$ , 先随机找到一个$w=a^2-n$, 且$w^{(p-1)/2}=-1(\mod p)$, 那么$x=(a+ \sqrt{w})^{(p+1)/2}$就是一个合法的$x$;

因为$x^2 = (a+\sqrt{w})^{p+1}=a^2-w=n$


## 代码

```c++
#include<cstdio>
#include<cstring>
#include<algorithm>
#include<cmath>
 
using namespace std;
 
typedef long long LL;
 
const int Hm = 4705219;
 
unsigned int seed;
int C,p,w,r2,s5,L,A;
 
struct Hmap {
 
	int info[Hm][2],nxt[Hm],head[Hm],use[Hm],cnt,tot; 
	void clear() {
		for(int i = 1;i <= cnt;i ++) head[use[i]] = 0;
		cnt = 0, tot = 0;
	}
 
	void push(int mi,int va) {
		int p = va % Hm;
		for(int i = head[p];i;i = nxt[i])
			if (info[i][0] == va) return;
		info[++ tot][0] = va,info[tot][1] = mi;
		if (!head[p]) use[++ cnt] = p;
		nxt[tot] = head[p],head[p] = tot;
	}
 
	int find(int va) {
		int p = va % Hm;
		for(int i = head[p];i;i = nxt[i])
			if (info[i][0] == va) return info[i][1];
		return -1;
	}
} hsh[2];
 
struct Z {
	LL a,b; Z(void){}
	Z(LL a,LL b) : a(a),b(b){}
};
 
Z operator *(const Z &a,const Z &b) {
	return Z((a.a * b.a % p + a.b * b.b % p * w % p) % p,(a.a * b.b % p + a.b * b.a % p) % p);
}
 
void updt(int &a,int b) {
	if (b == -1) return;
	if (a == -1) a = b; else
		a = min(a,b);
}
 
unsigned RAND() {
	return (seed = (seed * 31 + 998244353));
}
 
int fast(int a,int b) {
	LL as = 1;
	for(;b;b >>= 1) {
		if (b & 1) as = as * a % p;
		a = a * 1ll * a % p;
	}
	return as;
}
 
Z fast(Z a,int b) {
	Z as = Z(1,0);
	for(;b;b >>= 1) {
		if (b & 1) as = as * a;
		a = a * a;
	}
	return as;
}
 
int lerend(int n) {
	if (fast(n,(p - 1) / 2) == 1) return 1;
	return -1;
}
 
int _sqrt(int n) {
	if (!n) return 0;
	if (lerend(n) == -1) return -1;
	int a;
	while (1) {
		a = RAND() % p;
		w = (a * 1ll * a - n + p) % p;
		if (lerend(w) == -1) break;
	}
	Z cur = Z(a,1);
	cur = fast(cur,(p + 1) / 2);
	return cur.a;
}
 
int _solve(int sig,int tar) {
	int mi = -1,least = fast(fast(A,L),p - 2);
	for(int i = 0,tmp = tar;i <= p / L;i ++,tmp = 1ll * tmp * least % p) {
		int fr = i * L,v = hsh[(sig - (fr & 1) + 2) % 2].find(tmp);
		if (v == -1) continue;
		updt(mi,v + fr);
	}
	return mi;
}
 
void work() {
	scanf("%d%d", &C, &p);
	r2 = fast(2,p - 2),s5 = _sqrt(5);
	L = sqrt(p);
	A = (1 + s5) % p * 1ll * r2 % p;
	hsh[0].clear();
	hsh[1].clear();
	for(int i = 0,tmp = 1;i < L;i ++,tmp = tmp * 1ll * A % p)
		hsh[i & 1].push(i,tmp);
	C = C * 1ll * s5 % p;
	int Ans = -1;
	for(int od = 0,sig = 1;od < 2;od ++,sig *= -1) {
		int delta = (C * 1ll * C % p + 4 * sig % p + p) % p;
		delta = _sqrt(delta);
		if (delta == -1) continue;
		updt(Ans,_solve(od,(C + delta) % p * 1ll * r2 % p));
		updt(Ans,_solve(od,(C - delta + p) % p * 1ll * r2 % p));
	}
	printf("%d\n", Ans);
}
 
int main() {
	seed = 17;
	int T;
	scanf("%d", &T);
	for(;T;T --) work();
	return 0;
} 
```	

# 【CC SEALCM】Sereja and LCM
## 题目大意

对于$k \in [l,r]$, 其中满足$k|LCM(a_1,a_2...a_n)$且$max(a_i)<=m$的数列个数。 

## 解题报告

这是一个千载难逢的水题呀！

因为$m$和$k$的范围都很小， 所以质因子的个数最多只有$5$个； 

用一个$2^5$表示$k$的每个质因子是否满足的状态， 然后压到矩阵里， 大力转移就可以了。

## 代码

```c++
#include <bits/stdc++.h> 
 
using namespace std; 
 
#define REP(i,a,b) for(int i=int(a),nn=int(b);i<=nn;++i)
#define VEP(i,a,b) for(int i=int(a),nn=int(b);i>=nn;--i)
#define rep(i,b) for(int i=0,nn=int(b);i<nn;++i) 
#define B(x) (1<<(x))
 
typedef long long ll; 
 
const int p = 1000000007;
 
int n, test, m, l, r, t, a[5]; 
 
inline void divide(int d) {
 	t = 0; memset(a, 0, sizeof(a)); 
 	for (int i=2; i*i<=d; ++i) 
 		if (d % i==0) { 
 			a[t++]=1; while(d%i==0)a[t-1]*=i,d/=i;  
 		} 
  	if (d > 1) a[t++]=d, d=1; 
} 
  
struct matrix { 
	int x, y; ll a[20][20]; 
	matrix(){x=y=0; memset(a,0,sizeof(a));}
	matrix(int x,int y,int t=0):x(x),y(y) {
		memset(a, 0, sizeof(a)); 
		if(t) rep(i, x) a[i][i]=1; 
	} 
	matrix operator * (const matrix b) const { 
		matrix c(x, b.y) ;
		rep(i, c.x) rep(j, c.y) rep(k, y)
			c.a[i][j]=(c.a[i][j]+a[i][k]*b.a[k][j]%p)%p; 
		return c; 
	}
	void print() { 
		cout << x << ", " << y << endl; 
		rep(i, x) rep(j, y) cout<<a[i][j]<<((j+1==y)?'\n':' '); 
	} 
	
} ;
 
matrix _fast(matrix x, int k) {
	matrix as(x.x, x.y, 1); 
	for (;k; k>>=1, x=x*x)if(k&1)as=as*x;
	return as; 
} 
 
int main() { 
	scanf("%d", &test); 
	while (test --) { 
		scanf("%d%d%d%d", &n, &m, &l, &r); 
		int as = 0; 
		REP(d, l, r) { 
			divide(d); int w=B(t); matrix trs(w, w);  
			rep(i, w) REP(x,1,m) { int S=i; 
				rep(j, t) if (x%a[j]==0) S|=B(j); 
				trs.a[S][i]++;
			}
//			trs.print(); 
			trs=_fast(trs, n); 
			matrix ini(w, 1); ini.a[0][0]=1; 
			ini = trs*ini; 
			as += ini.a[w-1][0]; if(as>=p)as-=p; 
		} 
		printf("%d\n", as); 
	} 
	return 0; 
} 
```

# 【CC QTREE】Queries on tree again!
## 题目大意

环套树上两点的最短路径的最大连续字段和；

## 解题报告

强行树->环套树 。。。

首先把环拆开， 变成一棵树+一条边， 对于每个操作$x,y$进行特判， 判断是否走零散的那一条边。

然后用线段树，树链剖分后，对于每一个区间维护最大连续子段和，最小连续子段和，前缀最大子段和， 后缀最大子段和，前缀最小子段和， 后缀最小子段和， 区间变号标记。

修改的话非常的normal，就是直接树链剖分+线段树区间修改，查询的话比较麻烦， 因为需要完成从深到浅+由浅入深两个过程， 对应线段树的两种区间查询， 写的时候就比较冗长。

## 代码

```c++
#include <bits/stdc++.h> 
 
using namespace std; 
 
#define REP(i,a,b) for(int i=int(a),nn=int(b);i<=nn;++i) 
#define VEP(i,a,b) for(int i=int(a),nn=int(b);i>=nn;--i) 
#define rep(i,b) for (int i=0,nn=int(b);i<nn;++i) 
 
typedef long long ll; 
 
const int N=100001; 
const int S=400001;
struct edge { int nxt,to,c; 
	edge(int nxt=0,int to=0,int c=0):nxt(nxt),to(to),c(c){}
} e[N<<1]; 
int n,hed[N],tot,va[N],sn[N],sz[N],de[N],tp[N],f[N],dn[N],bd[N],cnt,u,v,w; 
int mnl[S],mnr[S],mn[S],mxl[S],mxr[S],mx[S],sm[S],rv[S],as,asr,m; 
 
inline void add(int x,int y,int c) {e[++tot]=edge(hed[x],y,c),hed[x]=tot;}
inline void in(int &x) { char ch=getchar(); int f=1; 
	for (;ch<'0'||ch>'9';ch=getchar()) if(ch=='-')f=-1; 
	for (x=0;ch>='0'&&ch<='9';ch=getchar())x=x*10+ch-48; x*=f; 
} 
void dfs(int x) { 
	de[x]=de[f[x]]+1,sz[x]=1,sn[x]=0; int y;
	for (int i=hed[x];i;i=e[i].nxt) if(y=e[i].to,!sz[y]) 
		va[y]=e[i].c,f[y]=x,dfs(y),sz[x]+=sz[y],sn[x]=((sz[y]>sz[sn[x]])?y:sn[x]); 
}
void odr(int x,int top) {
	tp[x]=top,dn[x]=++cnt,bd[cnt]=x; 
	if (sn[x]) { odr(sn[x],top); int y;
		for (int i=hed[x];i;i=e[i].nxt) 
			if(y=e[i].to,de[y]==de[x]+1&&y!=sn[x]) odr(y,y); 
	}
}
inline int _lca(int x, int y) { 
	int fx=tp[x], fy=tp[y]; 
	while (fx!=fy) { if (de[fx]<de[fy]) swap(fx,fy),swap(x,y); 
		x=f[fx], fx=tp[x];
	}  if (de[x]>de[y]) swap(x,y); return x; 
} 
inline int _dis(int x, int y) { return de[x]+de[y]-de[_lca(x,y)]*2;}
 
inline void _rev(int x) { 
	sm[x]=-sm[x],rv[x]^=1,swap(mn[x],mx[x]),swap(mnl[x],mxl[x]),swap(mnr[x],mxr[x]); 
	mn[x]=-mn[x],mx[x]=-mx[x],mnl[x]=-mnl[x],mnr[x]=-mnr[x],mxl[x]=-mxl[x],mxr[x]=-mxr[x]; 
}
inline void up(int x) { 
	sm[x]=sm[x<<1]+sm[x<<1|1], mx[x]=max(mxr[x<<1]+mxl[x<<1|1],max(mx[x<<1],mx[x<<1|1])); 
	mxl[x]=max(mxl[x<<1],sm[x<<1]+mxl[x<<1|1]); 
	mxr[x]=max(mxr[x<<1|1],sm[x<<1|1]+mxr[x<<1]); 
	mn[x]=min(mnr[x<<1]+mnl[x<<1|1],min(mn[x<<1],mn[x<<1|1])); 
	mnl[x]=min(mnl[x<<1],sm[x<<1]+mnl[x<<1|1]); 
	mnr[x]=min(mnr[x<<1|1],sm[x<<1|1]+mnr[x<<1]); 
}
inline void down(int x) { if (rv[x]) _rev(x<<1),_rev(x<<1|1),rv[x]=0;} 
void build(int x,int l,int r) {
	if (l==r) {  mn[x]=mnl[x]=mnr[x]=mx[x]=mxl[x]=mxr[x]=sm[x]=va[bd[l]];} 
	else { int mid=(l+r)>>1; build(x<<1,l,mid),build(x<<1|1,mid+1,r); up(x);} 
} 
void _rivers(int x,int l,int r,int _l,int _r) { 
	if (_l<=l&&r<=_r) _rev(x);
	else { int mid=(l+r)>>1;  down(x); 
		if (_l<=mid)_rivers(x<<1,l,mid,_l,_r); 
		if (_r>mid)_rivers(x<<1|1,mid+1,r,_l,_r); up(x); 
	} 
}
void _qry_r(int x,int l,int r,int _l,int _r) { 
	if (_l<=l&&r<=_r) { if(mx[x]>as)as=mx[x]; if(asr+mxr[x]>as)as=asr+mxr[x];
		asr+=sm[x]; if (mxl[x]>asr) asr=mxl[x]; 
	} else { int mid=(l+r)>>1; down(x); 
		if (_r>mid) _qry_r(x<<1|1,mid+1,r,_l,_r); 
		if (_l<=mid) _qry_r(x<<1,l,mid,_l,_r); 
	} 
} 
void _qry_l(int x,int l,int r,int _l,int _r) { 
	if (_l<=l&&r<=_r) { if(mx[x]>as)as=mx[x]; if(asr+mxl[x]>as)as=asr+mxl[x]; 
		asr+=sm[x]; if (mxr[x]>asr) asr=mxr[x]; 
	} else { int mid=(l+r)>>1; down(x); 
		if (_l<=mid) _qry_l(x<<1,l,mid,_l,_r); 
		if (_r>mid) _qry_l(x<<1|1,mid+1,r,_l,_r); 
	}
}
inline void reverse(int x,int y) { 
	int fx=tp[x], fy=tp[y]; while (fx!=fy) { 
		if (de[fx]<de[fy]) swap(fx,fy),swap(x,y); 
		_rivers(1,1,n,dn[fx],dn[x]),x=f[fx],fx=tp[x]; 
	} 
	if (de[x]>de[y]) swap(x,y); 
	if (x!=y) _rivers(1,1,n,dn[x]+1,dn[y]); 
} 
inline void query(int x, int y) { int lca=_lca(x,y); 
	while (tp[x]!=tp[lca]) _qry_r(1,1,n,dn[tp[x]],dn[x]),x=f[tp[x]]; 
	if (x!=lca) _qry_r(1,1,n,dn[lca]+1,dn[x]); 
	static int _l[N],_r[N],top; top=0; 
	while (tp[y]!=tp[lca]) _l[++top]=dn[tp[y]],_r[top]=dn[y],y=f[tp[y]]; 
	if (y!=lca) _l[++top]=dn[lca]+1,_r[top]=dn[y]; 
	while (top) _qry_l(1,1,n,_l[top],_r[top]),--top;
} 	
	
int main() { 
	in(n); int x,y,c;char type; rep(i,n) in(x),in(y),in(c),add(x,y,c),add(y,x,c); 
	dfs(1), odr(1,1); REP(x,1,n) for(int i=hed[x];i;i=e[i].nxt)
		if (y=e[i].to, f[x]!=y&&f[y]!=x) u=x,v=y,w=e[i].c; 
//	cout << u<<"->"<<v<<": "<<w<<endl; 
//	REP(i,1,n) cout << f[i] <<"->" <<i<<": "<<tp[i]<<", "<<va[i]<<endl;
	for(build(1,1,n),in(m); m; --m) { 
		type=getchar(); while(type!='?'&&type!='f')type=getchar(); 
		in(x),in(y);  if (type=='f') { 
			if (_dis(x,u)+_dis(y,v)>_dis(x,v)+_dis(y,u))swap(x,y); 
			if (_dis(x,y)<_dis(x,u)+_dis(y,v)+1) reverse(x,y); 
			else reverse(x,u), w=-w, reverse(v,y); 
		} else { 
			as=0, asr=0; if (_dis(x,u)+_dis(y,v)>_dis(x,v)+_dis(y,u)) swap(x,y); 
			if (_dis(x,y)<_dis(x,u)+_dis(y,v)+1) query(x,y); 
			else { query(x,u);if(w>as)as=w;if(asr+w>as)as=asr+w; 
				asr+=w; if(w>asr) asr=w; query(v,y);
			} 
			printf("%d\n", as); 
		} 
	}
	return 0; 
} 		 
```

# 【CC PERMUTE】Just Some Permutations 3
## 题目大意

长度为$n$的，相邻两个数的和不超过的$m$的合法排列个数。 

## 解题报告

考虑从大到小插入$n$个数， 每次插入一个数$n$， 需要用$\leqslant m-n$的数将他包裹起来， 也就是形如$XnY$, 或者$|nX$,或者$Xn|$的形式。 

现在考虑有几个可以选择来包裹$n$的数，开始显然是$m-n$个， 如果使用形如$|nX$、$Xn|$的形式进行包裹， 那么$X$将不能再被使用， 成为左右边界$|$的一部分，如果使用形如$XnY$的形式包裹， 那么$XnY$可以看成$X'$, 一个新的用来包裹的“数”； 

可以发现， 随着插入$n$,$n-1$...$\frac{m+1}{2}+1$ , 用来包裹的数的个数是不变的。

1. 如果$m$为奇数， 那么$(m+1)/2$插入时， 剩余的数的个数是$m-n+1$， 方案数是$(m-n+1)!$,总的方案数是$(2k+k(k-1))^{n-\frac{m+1}{2}}*(m-n+1)!$
2. 如果$m$为偶数， 那么$(m+1)/2$插入时， 剩余的数的个数是$m-n$, 方案数是$(2k+k(k-1))^{n-\frac{m+1}{2}}*(m-n)!$

其中$k=m-n$, 表示可以用来包裹的数的个数。 

## 代码

```c++
#include <bits/stdc++.h>
 
using namespace std; 
 
#define REP(i,a,b) for(int i=int(a),nn=int(b);i<=nn;++i) 
#define VEP(i,a,b) for(int i=int(a),nn=int(b);i>=nn;--i) 
#define rep(i,b) for(int i=0,nn=int(b);i<nn;++i) 
 
typedef long long ll; 
 
const int p =1000000007;
 
int n,m,T,k,fac[1000001]; 
 
inline ll fast(ll x,int k) {ll as=1;for(;k;k>>=1,x=x*x%p)if(k&1)as=as*x%p;return as;}
 
int main() { 
	scanf("%d", &T); 
	fac[0]=1; REP(i,1,1000000)fac[i]=1ll*fac[i-1]*i%p;
	while (T--) { 
		scanf("%d%d", &n, &m); k=m-n; 
		if (m & 1) printf("%d\n", fast(1ll*k*(k+1)%p,n-(m+1)/2)*fac[k+1]%p); 
		else printf("%d\n", fast(1ll*k*(k+1)%p,n-m/2)*fac[k]%p); 
	}
	return 0; 
} 
```