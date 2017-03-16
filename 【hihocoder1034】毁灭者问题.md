---
title: 【hihocoder 1034】毁灭者问题
date: 2017/01/01 22:00:00
tags: 
  - 平衡树
categories: 
  - 数据结构

---

> 一个暴力的数据结构题

<!--more-->

# 题目
<center><h2>hihocoder1034 : 毁灭者问题</h2>
时间限制:10000ms
单点时限:1000ms
内存限制:256MB</center>

## 描述
在 Warcraft III 之冰封王座中，毁灭者是不死族打三本后期时的一个魔法飞行单位。
毁灭者的核心技能之一，叫做魔法吸收（Absorb Mana）：
![图片标题](http://media.hihocoder.com//problem_images/20140715/14054064004625.png)
<!--more-->
现在让我们来考虑下面的问题：
假设你拥有$n$个魔法单位，他们从左到有站在一行，编号从$1$到$n$。 每个单位拥有三项属性：
$s_i$: 初始法力。
$m_i$: 最大法力上限。
$r_i$: 每秒中法力回复速度。
现在你操纵一个毁灭者，有$m$个操作,`t l r`，表示时刻$t$，毁灭者对所有编号从 $l$到 $r$ 的单位，使用了魔法吸收。操作按照时间顺序给出，计算毁灭者一共吸收了多少法力。

## 输入
输入数据的第一行有一个整数 $n(1 ≤ n ≤10^5)$ — 你的魔法单位的数目。

接下来的$n$行，每行有三个整数$ s_i,m_i,r_i(0≤s_i≤m_i≤10^5, 0≤r_i≤10^5)$ 描述一个魔法单位。

接下来一行又一个整数 $m(1 ≤ m ≤ 10^5)$， — 操作的数目。

接下来的 $m$ 行，每行描述一个操作 $t, l, r(0 ≤ t ≤ 10^9, 1 ≤ l ≤ r ≤ n)$，$t$ 非降。

## 输出
输出一行一个整数表示毁灭者一共吸收了多少法力。

# 解题报告

* 首先考虑线段树，但是每个点情况都不一样的情况下，区间修改不知道怎么搞;

* 因为每个点都是不同的，考虑分别计算每个点对最终答案的贡献：
    *  将询问区间分别按照左右端点排序，然后扫描，碰到左端点加入，右端点删除，就可以得到对每个节点询问的时间；
    *  容易发现，去除当前点第一次询问的贡献，剩余时间段中是否到达$m$可以通过讨论$\Delta t$和$\frac{m}{r}$,如果$\Delta t <= \frac{m}{r}$,对答案的贡献是$\Sigma_{\Delta t(\Delta t < = \frac{m}{r})}*r$,否则对答案的贡献是$Times_{(\Delta t>frac{m}{r})}*m$；
    *  我的做法是用两棵Splay，分别维护$t$和$\Delta t$,每次添加、删除$t$，修改相应$\Delta t$（$-1+2$或者$-2+1$或者左右端特判）；
    *  每次查询的时候，拿出$t$中第一个点特判，然后分别查询$\Delta t$中$<\frac{m}{r}+1$的$sum$和$sz$,然后贡献到答案中；

* 如果上天再给我一次机会，我一定用set+bit艹这个题！

# 代码
```c++
#include<iostream>
#include<algorithm>
#include<cstdio>
#include<cstring>
#include<cstdlib>
using namespace std;
const int N=200005;
const int inf=1500000000;
int n,m,id1[N],id2[N];//id1右端点排序，id2左端点排序 
char *cp=(char *) malloc(10000000); 
long long ans=0;
struct Query {
	int t,l, r; 
	Query(int t=0,int l=0,int r=0) 
		:t(t),l(l), r(r) {}
} query[N];
struct Obj {
	int s,m,r; 
	Obj(int s=0,int m=0,int r=0) 
		:s(s),m(m),r(r) {}
} obj[N];
struct Splay { //单点加、删、询问小于个数，pre,next 
	int son[N][2],f[N],val[N],sz[N],ts[N],root,cnt,sm[N];
	inline void newnode(int &x,int v) { 
		x=++cnt, f[x]=son[x][0]=son[x][1]=0; 
		sm[x]=val[x]=v, sz[x]=ts[x]=1; 
	}
	inline void cut(int x) {
		 sm[x]=sz[x]=f[x]=son[x][1]=son[x][0]=ts[x]=val[x]=0;
	}
	Splay() {
		root=sz[0]=son[0][0]=son[0][1]=f[0]=val[0]=sm[0]=0;
		newnode(root,0),newnode(son[root][1],inf); 
		f[cnt]=root,sz[root]=2;
	}
	inline void pushup(int x) {
		sz[x]=ts[x],sm[x]=val[x]*ts[x];
		if (son[x][0]) sz[x]+=sz[son[x][0]],sm[x]+=sm[son[x][0]]; 
		if (son[x][1]) sz[x]+=sz[son[x][1]],sm[x]+=sm[son[x][1]]; 
	}
	inline void rotate(int x) {
		int y=f[x],z=f[y],d=son[y][1]==x; 
		f[son[x][d^1]]=y,son[y][d]=son[x][d^1]; 
		if (z) son[z][son[z][1]==y]=x; 
		f[x]=z, son[x][d^1]=y,f[y]=x; 
		pushup(y),pushup(x); 
	}
	inline void splay(int x,int goal=0) {
		for (int y=f[x];y^goal;rotate(x),y=f[x]) 
			if (f[y]^goal)
			rotate((son[y][1]==x^son[f[y]][1]==y)?x:y); 
		if (!goal) root=x;
	}
	inline int findsz(int v) {//询问小于个数 
		int x=root,	ans=0;
		while (x) {
		 	if (val[x]>=v) x=son[x][0];
			else ans+=(sz[x]-sz[son[x][1]]),x=son[x][1];
		}	return ans;
	}
	inline long long findps(int v) { 
		int x=root, ans=0; 
		while (x) {
			if (val[x]>=v) x=son[x][0]; 
			else ans+=(sm[x]-sm[son[x][1]]),x=son[x][1]; 
		}	return ans; 
	}
		
	inline void pre(int goal) {
		int x=son[root][0]; 
		while (son[x][1]) x=son[x][1]; 
		splay(x,goal); 
	}
	inline void next(int goal) { 
		int x=son[root][1]; 
		while (son[x][0]) x=son[x][0]; 
		splay(x,goal); 
	}
	inline void findkth(int k,int goal=0) { 	
		int x=root; 
		while (sz[son[x][0]]^k) {
			if (sz[son[x][0]]>k) x=son[x][0]; 
			else k-=sz[son[x][0]]+ts[x],x=son[x][1];
		}
		splay(x,goal);
	}
			
	inline void findval(int v) {
		int x=root;
		while (1) { 
			if (!x) {printf("fuck\n");return;}
			if (val[x]==v) {splay(x);return;}
			if (val[x]>v) x=son[x][0]; 
			else x=son[x][1]; 
		}
	}
	inline void insert(int v) {
		int x=root,next; 
		while (1) {
			if (val[x]==v) {
				++ts[x],pushup(x),splay(x);return;
			}
			if (val[x]>v) next=0;else next=1; 
			if (!son[x][next]) {
				newnode(son[x][next],v),f[son[x][next]]=x; 
				splay(son[x][next]);return; 
			}
			x=son[x][next];
		}
	}
	inline void del(int v) { 
		findval(v); 
		if (ts[root]>1) {--ts[root];pushup(root);return;}
		if (!son[root][0]) 
			root=son[root][1],cut(f[root]),f[root]=0; 
		else if (!son[root][1]) 
				root=son[root][0],cut(f[root]),f[root]=0; 
			 else {
			 	pre(root),next(root); 
			 	son[son[root][0]][1]=son[root][1];
				f[son[root][1]]=son[root][0];
			 	pushup(son[root][0]),root=son[root][0]; 
			 	cut(f[root]), f[root]=0; 
			 }
	}
} delta,body; 		 
inline bool cmp1(int a,int b) {
	return query[a].r<query[b].r;
}
inline bool cmp2(int a,int b) { 
	return query[a].l<query[b].l; 
}
inline void insert(int t) {
	body.insert(t);if (body.ts[body.root]>1) return;body.pre(body.root),body.next(body.root); 
	if (body.val[body.son[body.root][0]]) delta.insert(t-body.val[body.son[body.root][0]]); 
	if (body.val[body.son[body.root][1]]!=inf) delta.insert(body.val[body.son[body.root][1]]-t); 
	if (body.val[body.son[body.root][0]]&&body.val[body.son[body.root][1]]!=inf)
		delta.del(body.val[body.son[body.root][1]]-body.val[body.son[body.root][0]]); 
}
inline void away(int t) {
	body.findval(t); if (body.ts[body.root]>1) {body.ts[body.root]--;return;}
	body.pre(body.root),body.next(body.root);
	if (body.val[body.son[body.root][0]]) delta.del(t-body.val[body.son[body.root][0]]); 
	if (body.val[body.son[body.root][1]]!=inf) delta.del(body.val[body.son[body.root][1]]-t); 
	if (body.val[body.son[body.root][0]]&&body.val[body.son[body.root][1]]!=inf)
		delta.insert(body.val[body.son[body.root][1]]-body.val[body.son[body.root][0]]); 
	body.del(t);
}
inline void calculate(int s,int m,int r) { 
	if (body.sz[body.root]==2) return; 
	body.findkth(1); ans+=min(1LL*m,1LL*(s+1LL*body.val[body.root]*r));
	if (delta.sz[delta.root]==2||!r) return; 
	int tmp=delta.findsz(m/r+1);
	ans+=1LL*(delta.sz[delta.root]-tmp-1)*m;
	ans+=1LL*delta.findps(m/r+1)*r;
}
inline void in(int &x) {
	for (;*cp<'0'||*cp>'9';cp++); 
	for (x=0;*cp>='0'&&*cp<='9';cp++) 
		x=x*10+*cp-48; 
}
int main() {
	freopen("hiho1034.in","r",stdin); 
	freopen("hiho1034.out","w",stdout); 
	fread(cp,1,10000000,stdin);	in(n); int i,j=0,k=0; 
	for (i=1;i<=n;++i) in(obj[i].s),in(obj[i].m),in(obj[i].r); 
	for (in(m),i=1;i<=m;++i) 
		in(query[i].t),in(query[i].l),in(query[i].r); 
	for (i=1;i<=m;++i) id1[i]=i,id2[i]=i; 
	sort(id1+1,id1+1+m,cmp1);sort(id2+1,id2+1+m,cmp2);
	for (i=1,j=0,k=0;i<=n;++i) { 
		while (j<m&&query[id2[j+1]].l==i) 
			insert(query[id2[++j]].t); 
		calculate(obj[i].s,obj[i].m,obj[i].r); 
		while (k<m&&query[id1[k+1]].r==i) 
			away(query[id1[++k]].t); 
	}
	printf("%lld\n",ans); 
}
```
![精神AC](https://leanote.com/api/file/getImage?fileId=57ec8256ab644107bb000b68)