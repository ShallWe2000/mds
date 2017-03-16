---
title: codechef-选做-3
date: 2017-03-14 13:48:10
tags:
  - LCT
  - 平衡树
  - 网络流
  - meet in the middle
  - DP
  - 后缀自动机
  - 博弈
  - 数位DP
categories:
  - 题目集锦
---

> 仍然在做CC中简单一点的题目

<!--more-->

# 【CC GERALD07】Chef and Graph Queries

## 题目大意

保留区间$[l,r]$的边，问图的联通块个数。 

## 解题报告

这个题目是一个动态树裸题，结合扫描线+离线。

使用LCT维护一个最大时间生成树，如果当前边相连的两个点是不连通的，那么直接加入当前边， 否则每次查询一条路径的最早加入的边， 删除后再加入当前边。

用一个树状数组维护单点修改区间和，边加入删除的时候在树状数组上也加入删除。  

然后联通块数=点数-边数；

## 代码

```c++
#include <bits/stdc++.h> 
 
using namespace std; 
 
#define REP(i,a,b) for(int i=int(a),nn=int(b);i<=nn;++i)
#define VEP(i,a,b) for(int i=int(a),nn=int(b);i>=nn;--i) 
#define rep(i,b) for(int i=0,nn=int(b);i<nn;++i) 
#define mp make_pair
#define fi first
#define se second
 
typedef long long ll; 
typedef pair<int,int> pii; 
 
const int N=200010;
const int INF=1000000000;
 
struct query { 
	int l,r,id; query(int l=0,int r=0,int id=0):l(l),r(r),id(id){}
} qry[N]; 
	
int n,T,m,q,bt[N],as[N]; 
pii pth[N];
 
template <typename T> 
inline void cmin(T &x, T a) { if (a<x) x=a; }
 
template <typename T>
inline void in(T &x) { char ch=getchar();
	int f=1; for(;ch<'0'||ch>'9';ch=getchar())if(ch=='-')f=-1;
	for (x=0; ch>='0'&&ch<='9';ch=getchar()) x=x*10+ch-48; x*=f; 
} 
 
namespace LCT { 
	int sn[N*2][2],f[N*2],mn[N*2],v[N*2],rv[N*2]; 
	
	inline void birth(int x) { 
		sn[x][0]=sn[x][1]=f[x]=rv[x]=0,mn[x]=v[x]=((x>n)?x:INF); 
	}
	inline int _d(int x) { return sn[f[x]][1]==x;}
	inline bool isroot(int x) {return !f[x]||(sn[f[x]][0]!=x&&sn[f[x]][1]!=x);}
	inline void up(int x) { 
		mn[x]=v[x]; 
		if (sn[x][0]) cmin(mn[x],mn[sn[x][0]]); 
		if (sn[x][1]) cmin(mn[x],mn[sn[x][1]]); 
	}
	inline void _rev(int x) {
		if (!x) return; swap(sn[x][0],sn[x][1]),rv[x]^=1;
	}
	inline void down(int x) { 
		if (rv[x]) _rev(sn[x][0]),_rev(sn[x][1]),rv[x]=0; 
	}
	inline void rotate(int x) {
	 	int y=f[x],d=_d(x),z=f[y]; 
	 	f[x]=z; if (!isroot(y)) sn[z][_d(y)]=x; 
	 	if (sn[x][d^1])f[sn[x][d^1]]=y; sn[y][d]=sn[x][d^1]; 
	 	sn[x][d^1]=y, f[y]=x, up(y); 
	} 
	inline void splay(int x) {
	 	static int stk[N],tp,nw; tp=0,nw=x; 
		while (!isroot(nw)) stk[++tp]=nw,nw=f[nw]; 
		stk[++tp]=nw; while(tp) down(stk[tp--]);
		for (int y=f[x];!isroot(x);rotate(x),y=f[x])
			if (isroot(y)) continue; 
			else rotate((_d(y)==_d(x)?y:x)); 
		up(x); 
	}
	inline void access(int x) {
	 	for (int r=0;x;r=x,x=f[x])splay(x),sn[x][1]=r; 
	}
	inline void mkroot(int x) { 
		access(x),splay(x),_rev(x); 
	}
	inline int root(int x) { 
		access(x),splay(x);while(sn[x][0])x=sn[x][0];return x;
	}
	inline void link(int x,int y) {
		mkroot(x),f[x]=y,access(x);
	}
	inline void cutf(int x) { 
		access(x),splay(x),f[sn[x][0]]=0,sn[x][0]=0,up(x);
	}
	inline void cut(int x,int y) { 
		mkroot(x), cutf(y); 
	} 
	inline int _min(int x,int y) { 
		mkroot(x),access(y),splay(y);return mn[y];
	}
	inline bool cnct(int x,int y) {
		return x!=y&&root(x)==root(y); 
	}
}
using namespace LCT; 
 
inline bool cmp(query a, query b) {return a.r<b.r;}
inline void add(int x,int v) {for(;x<=m;x+=x&(-x))bt[x]+=v;}
inline int qury(int x,int sm=0) {for(;x;x-=x&(-x))sm+=bt[x];return sm;}
 
int main() {
	for (in(T);T;--T) { 
		in(n),in(m),in(q); int x,y,z;
		REP(i,1,m) in(x),in(y),pth[i]=mp(x,y); 
		REP(i,1,q) in(x),in(y),qry[i]=query(x,y,i); 
		sort(qry+1, qry+1+q, cmp); int j=0; 
		REP(i,0,n+m) birth(i);REP(i,0,m)bt[i]=0; 
		REP(i,1,m) { x=pth[i].fi,y=pth[i].se;
			
			if (cnct(x,y)) { 
				z=_min(x,y); 
				cut(pth[z-n].fi,z),cut(pth[z-n].se,z), add(z-n, -1); 
			}
			if (x!=y) link(x,i+n),link(y,i+n),add(i,1); 
			if (qry[j+1].r==i) { int as2=qury(i),as1; 
				while (j<q&&qry[j+1].r==i) 
					as1=qury(qry[++j].l-1),as[qry[j].id]=as2-as1;
			}
		}
		REP(i,1,q) printf("%d\n",n-as[i]); 
	} 
	return 0; 
} 
```

# 【CC DIFTRIP】Different Trips

## 题目大意

每个点的权值是度数，问从每个点到根的所有字符串中，本质不同的子串个数。

## 解题报告

就用每个点的度数离散一下， 做了一个后缀自动机， 字符集的大小开的是$\sqrt{n}$的，因为本质不同的度数最多有$\sqrt{n}$个，然后果断跑得好慢（<del>比map还慢</del>)

貌似标算是树上的后缀数组。

广义后缀自动机是从这个地方学的：
[传送门]( "http://dwjshift.logdown.com/posts/304570")
## 代码

```c++
#include <bits/stdc++.h> 
 
using namespace std; 
 
#define REP(i,a,b) for(int i=int(a),nn=int(b);i<=nn;++i)
#define VEP(i,a,b) for(int i=int(a),nn=int(b);i>=nn;--i) 
#define rep(i,b) for(int i=0,nn=int(b);i<nn;++i) 
#define pb push_back
#define fi first
#define se second
 
typedef long long ll; 
typedef vector<int> vi; 
 
const int N=100010;
int ls[N],rt,ct,sn[N*2][320],pa[N*2],le[N*2],n,hh[N],tt;
vi nxt[N]; 
 
inline int birth(int len,int *go,int par=0) {
	++ct, le[ct]=len, pa[ct]=par; 
	if (go!=NULL) REP(i,1,tt) sn[ct][i]=go[i]; 
	return ct;
}
 
inline int ins(int x,int ls) {
	if (!hh[x])hh[x]=++tt; x=hh[x];
	int p=ls;  
	if (sn[p][x]&&le[sn[p][x]]==le[p]+1) 
		return sn[p][x];
	if (sn[p][x]) {
		int q=sn[p][x],nq=birth(le[p]+1,sn[q],pa[q]); 
		pa[q]=nq; while(p&&sn[p][x]==q)sn[p][x]=nq,p=pa[p]; 
		return nq; 
	}
	int np=birth(le[p]+1,NULL); 
	while (p&&!sn[p][x]) sn[p][x]=np,p=pa[p];
	if (!p) pa[np]=rt; else { 
		int q=sn[p][x]; 
		if (le[q]==le[p]+1) pa[np]=q; 
		else { 	
			int nq=birth(le[p]+1,sn[q],pa[q]); 
			pa[q]=pa[np]=nq;
			while(p&&sn[p][x]==q)sn[p][x]=nq,p=pa[p];
		}
	}
	return np;
}
	
void dfs(int x, int f) { 
	int sz=(f!=0);
	for (auto y:nxt[x])if(y!=f)++sz;
	ls[x]=ins(sz,ls[f]); 
	for (auto y:nxt[x])if(y!=f)dfs(y,x);
}
 
int main() { 
	scanf("%d", &n); int x,y;
	rep(i,n-1) scanf("%d%d",&x,&y),nxt[x].pb(y),nxt[y].pb(x);
	ls[0]=rt=ct=1, dfs(1,0); ll as=0; 
	REP(i,2,ct) as+=le[i]-le[pa[i]];
	return printf("%lld\n", as),0; 
} 
```

# 【CC EVILBOOK】Evil Book 

## 题目大意

每一个物品可以通过付出$C_i$的代价得到$M_i$的收益， 可以通过失去$X$点收益使某一个物品的$C_i$和$M_i$同时$\div 3$ , 问获得$666$点收益，最少付出多少的代价。 

## 解题报告

我竟然是一个连暴力都打不好的蒟蒻。。。

考虑最暴力的求解方法， 枚举每个物品选择的之前进行了多少次的$\div 3$操作，因为物品的收益和代价都是$\leqslant 10^7$的， 所以复杂度是$(\log_3 10^7)^n$, （<del>这就有点太暴力了</del>)

优化这个暴力？每个物品的$M$如果不小于$666*3$, 那么显然是很多余的（因为我们只需要$666$), 就可以必要的使用一次$\div 3$操作；

如果付出的$X$已经超过了可以获得的$M$？那就不可能选了对吧。题目中$X$有一个下界，是$10$, 可以发现， 如果对$666$进行$2$次$\div 3$操作，那么$2X>666 \div 3^2$,也就是$M$消到$666$后,最多进行$2$次$\div 3$操作。 

那么一个物品可能进行的$\div 3$操作数只有四种连续的情况，现在复杂度是$4^n$了，这就可以过了？ 再顺手加一个最优化剪枝， 最对进行的$\div 3$操作次数强加一个递增的顺序，就跑得比较快了。

## 代码

```c++
#include <bits/stdc++.h> 
 
using namespace std;
 
#define REP(i,a,b) for(int i=int(a),nn=int(b);i<=nn;++i) 
#define VEP(i,a,b) for(int i=int(a),nn=int(b);i>=nn;--i) 
#define rep(i,b) for(int i=0,nn=int(b);i<nn;++i) 
 
typedef double ff; 
const ff eps=1e-9; 
const int N=21; 
int c[N],m[N],n,X,pd[N];
ff sm,ans;
void dfs(int us,int pre,ff iu,ff hv,ff co){
//	cout<<hv<<' '<<co<<endl; 
	if (hv>666-eps){ ans=min(ans,co); return;}
	if (hv+sm/iu<666-eps) return;
	if (co>ans) return;
	if (us*X>hv) return;
	REP(i,pre,n) if (pd[i]==0&&us*X<m[i]/iu){
		pd[i]=1; sm-=m[i];
		dfs(us,i+1,iu,hv-us*X+m[i]/iu,co+c[i]/iu);
		pd[i]=0; sm+=m[i];
	}
	dfs(us+1,1,iu*3,hv,co);
}
void solve(){
	scanf("%d%d",&n,&X); sm=0;
	REP(i,1,n) scanf("%d%d",&c[i],&m[i]), sm+=m[i];
	if (sm<666) { puts("impossible"); return;}
	ans=1e18; dfs(0,1,1,0,0);
	printf("%.0lf\n",ans);
}
int main(){
	int t; scanf("%d",&t);
	for (;t;t--) solve();
	return 0;
}  
```
# 【CC LEBOXES】Little Elephant and Boxes

## 题目大意

$n$个宝箱，每一个宝箱有$p_i/100$的概率转换成$V_i$的dollars, 有$1-\frac{p_i}{100}$的概率转成$1$个diamond. 

$m$个商品，每个需要$C_i$的dollars和$D_i$的diamonds. 

问聪明至极的情况下， 买到的商品数的期望值。 

## 解题报告

首先$n$的范围十分的精妙，$2^n$跑不了， $2^{n/2}$次方没问题， 像极了meet in the middle的数据范围。

虽然dollars的数量级很大，但diamond很小， 考虑让diamond数和买到的商品数做下标， 具体的，$f[i][j]$表示使用$i$个钻石，买$j$个商品，最少需要多少的dollars, 这个可以用O/1背包解决。

对 $2^x$ ， 即前 $x$ 个商品的转换情况进行预处理，得到 $sto_{i,j}(x,y)$ , 表示得到$i$个钻石和不超过$x$个dollars，的概率为 $y$ .

对$2^{n-x}$， 即剩下的商品进行查询， 通过二分+查表， 得到期望值。 

复杂度$O(2^x\log 2^x+n^3+2^{n-x} \log 2^x)$

## 代码

```c++
#include <bits/stdc++.h>
 
using namespace std; 
 
#define REP(i,a,b) for(int i=int(a),nn=int(b);i<=nn;++i)
#define VEP(i,a,b) for(int i=int(a),nn=int(b);i>=nn;--i) 
#define rep(i,b) for(int i=0,nn=int(b);i<nn;++i) 
#define pb push_back
#define lobo lower_bound
#define fi first
#define se second 
#define mp make_pair
#define sz(x) (x.size()) 
typedef long long ll; 
typedef double ff; 
typedef pair<int,ff> pif; 
 
const int N=35; 
int les[N][N],T,n,m,B,v[N],c[N],d[N],p[N]; 
vector<pif> sto[N]; 
ff as; 
 
template <typename T> 
inline void cmin(T &x, T a) {if (a<x) x=a;}
 
void dfs(int x,int dos,int das,ff pr) { 
	if (x>B) { sto[das].pb(mp(dos,pr)); return;} 
	dfs(x+1,dos+v[x],das,pr*((ff)p[x]/100.00)); 
	dfs(x+1,dos,das+1,pr*(1.0-(ff)p[x]/100.00)); 
} 
 
void getans(int x,int dos,int das,ff pr) { 
	if (x>n) { 
		REP(i,0,B) REP(j,1,m) { 
			vector<pif>::iterator
				k=lobo(sto[i].begin(),sto[i].end(),
				mp(les[i+das][j]-dos,(ff)0)); 
			if (k!=sto[i].end()) { 
				ff pl=sto[i][sz(sto[i])-1].se; 
				if (k!=sto[i].begin()) --k,pl-=k->se; 
				as+=(ff)j*pl*pr; 
			}
		}
		return; 
	}
	getans(x+1,dos+v[x],das,pr*((ff)p[x]/100.00)); 
	getans(x+1,dos,das+1,pr*(1.0-(ff)p[x]/100.00)); 
} 
 
inline void prepare() {
//get les[i][j] : use i diamonds to buy j objects, what the least cost is; 
	memset(les,127,sizeof(les)); 
	les[0][0]=0; REP(i,1,m) VEP(j,n,d[i]) VEP(k,m,1) 
		if (les[j-d[i]][k-1]+c[i]<300000010) 
			cmin(les[j][k],les[j-d[i]][k-1]+c[i]); 
//get sto[i][j](fi,se) : use the first B boxes, get i dias and no more than fi dols, what the prob; 
	B=n/3*2; dfs(1,0,0,1); 
	REP(i,0,B) sort(sto[i].begin(),sto[i].end()); 
	REP(i,0,B) REP(j,1,sz(sto[i])-1) sto[i][j].se+=sto[i][j-1].se; 
} 
 
int main() {
	scanf("%d",&T); while (T--) { as=0; 
		rep(i,22) sto[i].clear(); 
		scanf("%d%d",&n,&m);  REP(i,1,n) scanf("%d%d",&v[i],&p[i]); 
		REP(i,1,m) scanf("%d%d",&c[i],&d[i]);  
		prepare(), getans(B+1,0,0,1); 
		printf("%.4lf\n",as); 
	} 
	return 0; 
} 
```

# 【CC GTHRONES】A Game of Thrones
## 题目大意

Bran 和 Tyrion 进行博弈， Bran先手，指定一个数作为当前数$x$， 之后交替进行操作： 选择一个与当前数满足条件的数$y$, 将$y$设为当前数， 删掉一个$x$, 不能操作者负。 

需要满足的条件：
1. 如果$x>y$ , $y|x$, $x/y$是一个质数；
2. 如果$y>x$ , $x|y$, $y/x$是一个质数。

## 解题报告

首先，如果把满足条件的点$x$,$y$连边， 那么原图会成为一个二分图，可以安装质因子个数的奇偶进行分类。 

如果这个二分图存在完美匹配，那么先手必负， 因为后手只需要选取匹配中与当前点相连的点就可以不负。

否则， 先手必胜， 考虑先手选择一个不在某个最大匹配中的点作为初始状态， 后手需要使用存在于最大匹配中的点（否则不满足最大匹配的定义），然后，先手只需要使用最大匹配中对应的点回敬即可不负。 

所以只需要判断是否存在完美匹配就可以判断胜负， 然后找到最小的不在某个最大匹配中的点作为最小代价，可以通过网络流的退流操作进行判断。 
 
## 代码
 
```c++
#include <bits/stdc++.h> 
 
using namespace std; 
 
#define REP(i,a,b) for(int i=int(a),nn=int(b);i<=nn;++i) 
#define VEP(i,a,b) for(int i=int(a),nn=int(b);i>=nn;--i) 
#define rep(i,b) for(int i=0,nn=int(b);i<nn;++i) 
#define pb push_back
 
typedef unsigned long long ll; 
typedef long double ff; 
 
const int N=510; 
const ll inf=1000000000000000ll; 
 
struct edge { 
	int nxt,to; ll f; 
	edge(int nxt=0,int to=0,ll f=0) 
		:nxt(nxt),to(to),f(f){}
} e[N*N*3]; 
 
int hed[N*2],tot=1,n,c[N*2],co[N*2],p[N],od[N]; bool viz[N];
int S,T,s,t,de[N*2],q[N*N],he,ta; vector<int> v[N]; ll u[N*2];  
 
inline int add(int x,int y,ll f) {
	e[++tot]=edge(hed[x],y,f), hed[x]=tot; 
	e[++tot]=edge(hed[y],x,0), hed[y]=tot; return tot^1; 
} 
template <typename T> 
inline void in(T &x) { char ch=getchar(); int f=1; 
	for (;ch<'0'||ch>'9';ch=getchar()) if(ch=='-')f=-1; 
	for (x=0; ch>='0'&&ch<='9';ch=getchar()) x=x*10+ch-48; 
	x*=f; 
} 
void _color(int x,int fm=0) { 
	viz[x]=1;co[x]=co[fm]^1; 
	rep(i,v[x].size()) if (!viz[v[x][i]]) _color(v[x][i],x); 
} 
 
inline bool bfs() { 
	he=ta=0,q[ta++]=s; memset(de,0,sizeof(de)),de[s]=1; int x,y;  
	while (he<ta) { x=q[he++]; for(int i=hed[x];i;i=e[i].nxt) 
		if (y=e[i].to, !de[y]&&e[i].f) de[y]=de[x]+1,q[ta++]=y; 
	} 
	return de[t]; 
} 
ll dfs(int x, ll lm) { ll f,hv=0; int y; if (x==t) return lm;  
	for (int i=hed[x];i;i=e[i].nxt) if (y=e[i].to,de[y]==de[x]+1&&e[i].f) { 
		f=dfs(y,min(e[i].f, lm)); e[i].f-=f,e[i^1].f+=f,hv+=f,lm-=f; 
		if (!lm) return hv; 
	} de[x]=-1; return hv; 
} 
 
inline ll dinic(int _s, int _t) { 
	s=_s, t=_t; ll as=0; while (bfs()) as+=dfs(s,inf); return as; 
} 
 
const int pm[10]={2,3,5,7,11,13,17,19,23,29}; 
inline ll mul(ll a,ll b,ll p) {
	a%=p,b%=p; if(a<0)a+=p; if(b<0)b+=p; 
	ll k=(ll)(((ff)a*b/p)); ll as=a*b-p*k; 
	as%=p; if(as<0) as+=p; return as; 
} 
inline ll fast(ll x,ll k,ll p) { 
	ll as=1; for(;k;k>>=1,x=mul(x,x,p)) if(k&1) as=mul(as,x,p);
	return as;
} 
inline bool miller_rabin(ll x) { 
	if (x==2) return 1; if (x<=1) return 0; if(x%2==0) return 0; 
	ll a,m=x-1,k=0; while (m%2==0) m/=2,++k; 
	for (int i=0; i<10;++i) { if (pm[i]==x) return 1; 
		a=fast(pm[i],m,x); if(a==1) continue; 
		int j; for (j=1;j<=k;++j) { 
			if(a==x-1) break; a=mul(a,a,x);
		} 
		if (j>k) return 0; 
	} 
	return 1; 
} 
 
inline bool chek(ll a,ll b) { 
	if (a<b) swap(a,b); if(a%b) return 0; 
	if (!miller_rabin(a/b)) return 0; return 1;
} 
inline bool cmp(int a,int b) { return u[a]<u[b];} 
 
int main() {
//	freopen("A.in","r",stdin); 
//	freopen("A.out","w",stdout); 
	in(n); REP(i,1,n) in(u[i]),in(c[i]);
	REP(i,1,n) REP(j,i+1,n) if(chek(u[i],u[j])) 
		v[i].pb(j),v[j].pb(i); 
	REP(i,1,n) if (!viz[i]) _color(i); 
 
	S=0, T=n+1; ll s1=0, s2=0;  
	REP(i,1,n) if (co[i]==1) { 
		p[i]=add(S,i,c[i]),s1+=c[i]; 
		rep(j,v[i].size()) add(i,v[i][j],inf); 
	} else p[i]=add(i,T,c[i]),s2+=c[i]; 
	ll flo=dinic(S,T); 
 
	if (flo==max(s1,s2)) {puts("Tyrion");return 0;}
	REP(i,1,n) od[i]=i; sort(od+1,od+1+n,cmp); 
	ll fno,old1,old2; int pa,k;  
	printf("Bran"); REP(i,1,n) { k=od[i],pa=p[k]; 
		if (e[pa].f) {printf(" %lld\n",u[k]);return 0;}
		old1=e[pa].f, old2=e[pa^1].f; 
		e[pa].f=e[pa^1].f=0; 
		if (co[k]==1) fno=dinic(S,k); 
		else fno=dinic(k,T); 
		if (fno) { printf(" %lld\n",u[k]); return 0;} 
		e[pa].f=old1, e[pa^1].f=old2;
	} 
	return 0; 
} 
```

# 【CC CARDCHUF】Card Shuffle
## 题目大意

初始一个$1$到$n$的序列，每次进行几个操作：
1. 拿走顶端的$A$个数。
2. 再拿走顶端的$B$个数。 
3. 把$A$个数放回顶端。
4. 在顶端拿走$C$个数。 
5. 把$B$个数倒序放回顶端。
6. 放回$C$个数。

## 解题报告

天朝数据结构碾压终生。

FHQTreap超容易写的一道裸题。HARD？

## 代码

```c++
#include <bits/stdc++.h> 
 
using namespace std; 
 
#define REP(i,a,b) for(int i=int(a),nn=int(b);i<=nn;++i)
#define VEP(i,a,b) for(int i=int(a),nn=int(b);i>=nn;--i) 
#define rep(i,b) for(int i=0,nn=int(b);i<nn;++i) 
#define pb push_back
#define mp make_pair
#define fi first
#define se second 
 
typedef long long ll; 
typedef double ff; 
typedef pair<int,int> droot; 
 
const int N=100100; 
 
int n,m,sn[N][2],ke[N],v[N],sz[N],ct,rv[N],rt;
 
inline void in(int &x) {
	char ch=getchar(); int f=1; 
	for(;ch<'0'||ch>'9';ch=getchar()) 
		if (ch=='-') f=-1; 
	for (x=0;ch>='0'&&ch<='9';ch=getchar()) 
		x=x*10+ch-48; x*=f; 
} 
 
inline int birth(int x) { 
	++ct; ke[x]=rand(),sn[ct][0]=sn[ct][1]=0,v[ct]=x,sz[ct]=1;return ct;
} 
 
inline void _rev(int x) {
	if (!x) return ; 
	rv[x]^=1; swap(sn[x][0],sn[x][1]); 
} 
 
inline void down(int x) { 
	if (rv[x]) _rev(sn[x][0]), _rev(sn[x][1]), rv[x]=0; 
} 
 
inline void up(int x) { 
	sz[x]=1; if (sn[x][0]) sz[x]+=sz[sn[x][0]]; 
	if (sn[x][1]) sz[x]+=sz[sn[x][1]]; 
} 
 
int build() { 
	static int stk[N],tp,x,las; 
	for (int i=1; i<=n; ++i) { 
		x=birth(i); while (tp&&ke[stk[tp]]>ke[x]) 
			las=stk[tp], --tp, up(las); 
		if (tp) sn[stk[tp]][1]=x;
		sn[x][0]=las, stk[++tp]=x,las=0;  
	} 
	while (tp) up(stk[tp]), --tp; 
	return stk[1]; 
} 
 
droot split(int x,int k) {
	if (!x) return mp(0,0); if (!k) return mp(0,x); 
	down(x); int lsz=sz[sn[x][0]]; droot y;  
	if (lsz>=k) { y=split(sn[x][0],k); 
	 	sn[x][0]=y.se; up(x); y.se=x; 
	} else { y=split(sn[x][1],k-1-lsz); 
	 	sn[x][1]=y.fi; up(x); y.fi=x;  
	} 
	return y; 
} 
 
int merge(int a, int b) { if(a*b==0) return a+b; 
	down(a), down(b); 
	if (ke[a]<ke[b]) { 
		sn[a][1]=merge(sn[a][1],b);
		up(a); return a; 
	} else { 
		sn[b][0]=merge(a, sn[b][0]); 
		up(b); return b; 
	} 
} 
 
void look(int x) { 
	if(!x) return; down(x); 
	if (sn[x][0]) look(sn[x][0]); 
	printf("%d ", v[x]); 
	if (sn[x][1]) look(sn[x][1]); 
} 
 
int main() { 
	in(n),in(m),rt=build();int a,b,c;
	//look(rt); //cout<<endl; 
	rep(i,m) { in(a),in(b),in(c); 
		droot A=split(rt,a); 
		droot B=split(A.se,b); 
		rt=merge(A.fi,B.se); 
		//look(rt); //cout<<endl; 
		A=split(rt, c); 
		_rev(B.fi),A.se=merge(B.fi,A.se); 
		rt=merge(A.fi,A.se); 
		//look(rt);cout<<endl; 
	} 
	look(rt); cout<<endl; return 0; 
} 
```

# 【CC TAPAIR】Counting The Important Pairs
## 题目大意

一个无向简单图， 求无序点对$(x,y)$满足同时删除这两条边，无向图不再联通的对数。

## 解题报告

在dfs树上观察， 发现如果$x$,$y$两条边有相同的编号，那么同时删除无向图不再联通。

编号是啥？ 令dfs树上的返祖边，每个反祖边有一个编号，那么树上一条边的编号是删除这条边后，使用哪些返祖边可以使图依旧联通的编号集合。

其实有两种情况：
1. 两个编号相同边一个是返祖边，一个是树边，那么这个返祖边是唯一一个使得树边删除后图依然联通的边， 同时删除满足题意。
2. 如果都是树边， 容易得到不存在返祖边使得两条树边之间的联通块与其他部分联通， 同时删除两条树边，使得中间的联通块“独立”， 满足题意。

那怎么判断两个边的编号是不是相同？用奇怪的Hash技巧（比如xor?)就可以了。

## 代码

```c++
#include <bits/stdc++.h> 
 
using namespace std; 
 
#define REP(i,a,b) for(int i=int(a),nn=int(b);i<=nn;++i)
#define VEP(i,a,b) for(int i=int(a),nn=int(b);i>=nn;--i)
#define rep(i,b) for(int i=0,nn=int(b);i<nn;++i) 
#define RAND ull((1ull*rand()<<45)|(1ull*rand()<<30)|(rand()<<15)|rand())
 
typedef long long ll; 
typedef unsigned long long ull; 
 
const int N=100010; 
const int M=300010; 
 
struct edge {
	 int nxt,to; edge(int nxt=0,int to=0)
	 	:nxt(nxt),to(to) {}
} e[M<<1]; 
 
int hed[N],tot,n,m; ull v[M];
bool vi[N],instk[N]; 
 
template <typename T> 
inline void in(T &x) { 
	char ch=getchar(); int f=1; 
	for (;ch<'0'||ch>'9';ch=getchar())
		if (ch=='-') f=-1; 
	for (x=0;ch>='0'&&ch<='9';ch=getchar()) 
		x=x*10+ch-48;  x*=f; 
} 
 
inline void add(int x,int y) { 
	e[++tot]=edge(hed[x],y), hed[x]=tot; 
} 
 
inline void dfs(int x,int fm) {
	instk[x]=1,vi[x]=1; 
//	cout<<RAND<<endl; 
	int y; for (int i=hed[x];i;i=e[i].nxt) 
		if (y=e[i].to, y!=fm) { 
			if (vi[y]) { if(instk[y]) v[++n]=RAND,v[x]^=v[n],v[y]^=v[n];}
			else dfs(y,x); 
		} 
	v[fm]^=v[x],instk[x]=0; 
} 
 
int main() { in(n),in(m); srand(23333); 
	int x,y; rep(i,m) in(x),in(y),add(x,y),add(y,x); 
	dfs(1,0),v[1]=v[n],--n, sort(v+1,v+1+n); ll t; int i=0,j=0; 
	while (i<n&&v[i+1]==0) ++i; ll as=0; 
	t=i,t*=n-1,as+=t,t=i,t=t*(t-1)/2,as-=t; 
	for (++i;i<=n;i=j+1) { j=i; while (j<n&&v[j+1]==v[i]) ++j; 
		t=j-i+1, t=t*(t-1)/2, as+=t; 
	} 
	cout<<as<<endl; 
	return 0; 
} 
```


# 【CC FARASA】Furik and Rubik and Sub Array
## 题目大意

给出一个序列，求不同的区间和的个数。 

## 解题报告

题目对$n*SUM$做了约束，如果$n$很小， 直接$O(n^2logn)$艹过去了。 

如果$n$大一点，比如$2000 < n \leqslant 20000$ , 那么数的和就可以作为数组的下标，所以$O(n^2)$艹过去了。

如果$n$再大一点？比如$20000<n\leqslant 200000$, 可以把$sum(l,r)$表示成$pre(r)-pre(l-1)$ , 可以把$pre$数组和$-pre$数组分别作为多项式的幂指数， 做多项式乘法统计系数就可以了，$O(sum\log sum)$艹过去了！

但是我发现自己的多项式代码一直WA，然后把对拍用的暴力交了上去，竟然是跑得最快的。。。 是不是用什么生日悖论可以证明暴力复杂度优越呀。。

## 代码

```c++ 
#include <bits/stdc++.h> 
 
using namespace std;
 
typedef long long ll;
 
const int N = 200010;
ll a[N]; int s[N],n;
bool f[20000010];
 
int main() {
	scanf("%d", &n);
	for(int i=1;i<=n;++i) scanf("%lld", &a[i]);
	ll ans=0;
	if(n<=2000) {
		set<ll> s;
		for(int l=1;l<=n;++l) {
			ll sum=0;
			for(int r=l;r<=n;++r) 
				sum+=a[r], s.insert(sum);
		}
		ans=s.size();
	} else {
		s[0]=0, f[0]=1;
		for(int i=1;i<=n; ++i)
			s[i]=s[i-1]+a[i], f[s[i]]=1;
		int sum=s[n];
		for(int i=1;i<=sum;++i) {
			if(!f[i]) {
				for(int j=1; j<=n; ++j) {
					int x=s[j]+i; if (x>sum) break;
					if(f[x]) { f[i]=1; break; }
				}
			}
			ans+=f[i];
		}
	}
	cout<<ans-1<<endl;	
	return 0;
}
```

FFT 做法

```c++
#include <bits/stdc++.h> 
 
using namespace std; 
 
#define REP(i,a,b) for(int i=int(a),nn=int(b);i<=nn;++i) 
#define VEP(i,a,b) for(int i=int(a),nn=int(b);i>=nn;--i) 
#define rep(i,b) for(int i=0,nn=int(b);i<nn;++i) 
 
typedef long long ll; 
typedef double ff; 
 
const ff pi=acos(-1); 
ll a[200010],n;
 
template <typename T> 
inline void in(T &x) { 
	char ch=getchar();
	for (;ch<'0'||ch>'9';ch=getchar());
	for (x=0;ch>='0'&&ch<='9';ch=getchar())
		x=x*10+ch-48; 
} 
 
namespace one { 
	void main() { 
		set<ll> s; 
		REP(i,1,n) { a[i]+=a[i-1]; 
			rep(j,i) if (s.count(a[i]-a[j])==0) 
				s.insert(a[i]-a[j]); 
		} 
		printf("%lld\n", s.size()-1); 
	} 
} 
namespace two { 
	void main() { 
		bitset<20000001> exi; 
		REP(i,1,n) { a[i]+=a[i-1]; 
			rep(j,i) exi.set(a[i]-a[j]); 
		}
		exi.reset(0); 
		printf("%lld\n", exi.count()-1) ; 
	} 
}
namespace three { 
	int _n,r[600000]; 
	struct cmx { 
		ff x,y; cmx(ff x=0,ff y=0):x(x),y(y){}
		cmx operator +(const cmx b) { return cmx(x+b.x,y+b.y);}
		cmx operator -(const cmx b) { return cmx(x-b.x,y-b.y);}
		cmx operator *(const cmx b) { return cmx(x*b.x-y*b.y,x*b.y+y*b.x);}
	} A[600000],B[600000]; 
	inline void dft(cmx *a,int f) { 
		rep(i,_n) if (i<r[i]) swap(a[i],a[r[i]]); 
		for (int m=1; m<_n; m<<=1) { 
			cmx wn=cmx(cos(pi/m),sin(pi/m)*f); 
			for (int i=0; i<_n; i+=m<<1) {
				cmx w=cmx(1,0); 
				for (int j=0;j<m; ++j) {
					cmx x=a[i+j],y=a[i+j+m]*w; 
					a[i+j]=x+y, a[i+j+m]=x-y; 
					w=w*wn; 
				} 
			} 
		} 
		if (f==-1) rep(i,_n) a[i].x/=(ff)_n; 
	} 
				
	void main() {
		REP(i,1,n) { a[i]+=a[i-1],A[a[i]].x=1;} 
		A[0].x=1; rep(i,a[n]+1) B[i].x=A[a[n]-i].x; 
		int m=a[n]+a[n]+1; for (_n=1; _n<m; _n<<=1); 
		for (int i=0,j=0;i<_n;++i) { r[i]=j; 
			for (int k=_n>>1;(j^=k)<k;k>>=1); 
		} 
		dft(A,1),dft(B,1);rep(i,_n) A[i]=A[i]*B[i]; ll as=0; 
		dft(A,-1); REP(i,a[n]+1,a[n]+a[n]) if (int(A[i].x+0.5)) ++as; 
		printf("%lld\n", as-1); 
	} 
} 
 
int main() {
//	freopen("A.in","r",stdin);
//	freopen("A.out","w",stdout); 
	scanf("%d", &n); //cout<<1<<endl;
	
	REP(i,1,n) scanf("%lld",&a[i]);
	
 
	if (n<=2000) one :: main(); 
	else if (n<=20000) two :: main(); 
		else three:: main(); 
	return 0; 
} 
```

# 【CC CHANGE】Making Change 
## 题目大意

## 解题报告

## 代码

```c++
#include <bits/stdc++.h> 
 
using namespace std; 
 
#define REP(i,a,b) for(int i=int(a),nn=int(b);i<=nn;++i) 
#define VEP(i,a,b) for(int i=int(a),nn=int(b);i>=nn;--i) 
#define rep(i,b) for (int i=0,nn=int(b);i<nn;++i) 
 
typedef unsigned int ui;
 
const int N=51; 
const int p=1000000007;
 
int d[N],c[101],n,T,le,tw[400],f[N*N*20],sm;
char s[110]; 
int main() { 
	scanf("%d", &T); 
	while (T--) { 
		sm=0,scanf("%d%s",&n,s); 
		le=strlen(s); rep(i,le) c[i]=s[le-i-1]-'0'; 
		int lf=0,tn=0; while (le>1||c[0]>=1) { 
			VEP(i,le-1,0) { int x=lf*10+c[i]; 
				lf=x%2, c[i]=x/2; 
			}
			tw[++tn]=lf, lf=0; 
			while (c[le-1]==0) le--; 
		}
//		REP(i, 1,tn) cout<<tw[i]<<' ' ;cout<<endl; 
		REP(i,1,n) scanf("%d", &d[i]),sm+=d[i];
		memset(f,0,sizeof(f)); 
		f[0]=1; int smx=0,v;  
		sort(d+1,d+1+n); 
		for (int i=1;i<=tn;++i) { 
			memset(f+1+smx,0,sm*sizeof(ui));
			for (int j=1;j<=n;++j) { v=d[j],smx+=v; 
				VEP(k,smx,v) f[k]+=f[k-v],f[k]-=p*((f[k]>>30)&3); 
			}
			int _j,j; 
			for(_j=tw[i],j=0; _j<=smx;_j+=2,++j) {  
				f[j]=f[_j]; if (f[j]>=p) f[j]-=p;
			} 
			smx=j-1; 
		}
		printf("%u\n", f[0]); 
	}
	return 0; 
} 
```