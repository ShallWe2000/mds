---
title: 【bzoj 4383】[POI2015]Pustynia
date: 2017/01/01 22:00:00
tags: 
  - 线段树
  - 拓扑排序
categories: 
  - 数据结构

---
> 线段树优化建图的裸题？

<!--more-->

# 题目
<center><h2>4383: [POI2015]Pustynia</h2>

Time Limit: 10 Sec  Memory Limit: 128 MBSec  Special Judge</center>
## Description

给定一个长度为$n$的正整数序列$a$，每个数都在$1$到$10^9$范围内，告诉你其中$s$个数，并给出$m$条信息，每条信息包含三个数$l,r,k$以及接下来$k$个正整数，表示$a[l],a[l+1],...,a[r-1],a[r]$里这$k$个数中的任意一个都比任意一个剩下的$r-l+1-k$个数大（严格大于，即没有等号）。
请任意构造出一组满足条件的方案，或者判断无解。
<!--more-->
##Input

第一行包含三个正整数$n,s,m$($1<=s<=n<=100000$，$1<=m<=200000$)。
接下来s行，每行包含两个正整数$p[i],d[i]$($1<=p[i]<=n，1<=d[i]<=10^9$)，表示已知$a[p[i]]=d[i]$，保证$p[i]$递增。
接下来m行，每行一开始为三个正整数$l[i],r[i],k[i]$($1<=l[i] < r[i]<=n，1<=k[i]<=r[i]-l[i]$)，接下来$k[i]$个正整数$x[1],x[2],...,x[k[i]]$($l[i]<=x[1] < x[2]<... < x[k[i]]<=r[i]$)，表示这$k[i]$个数中的任意一个都比任意一个剩下的$r[i]-l[i]+1-k[i]$个数大。$\sum k <= 300000$

## Output

若无解，则输出NIE。
否则第一行输出TAK，第二行输出$n$个正整数，依次输出序列$a$中每个数。

# 解题报告
很显然的一个做法：将所有大于关系$x>y$形象为$p_{y \to x}=1$，这样就可以得到每一个$a_i$的最小值.
这样的大于关系是$O(n^2)$的时空复杂度，难以承受；

因为每次连边的点可以构成一个区间而一个区间可以被拆成$\log n$个子区间，所以可以使用线段树优化连边；

具体做法是,在建立线段树的`build`过程中，树上连接$p_{son \to father}=0$,对于每个$m$，新建一个节点`trtot`，向$k$个数连接边权为$0$的边，然后线段树的`ins`操作，把区间子对应节点向`trtot`连
边权为$1$的边；然后跑拓扑排序就可以;

# 代码
```c++
#include<iostream>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<algorithm>
using namespace std;
const int inf=1000000000; 
const int N=400001; 
const int M=2000001;
struct E {
	int next,to,v; 
	E(int next=0,int to=0,int v=0)
		:next(next),to(to),v(v){}
} e[M];
int n,m,cnt,trtot,tot,rt,head[N],enter[N];
int ls[N],rs[N],pos[N],h[N],a[N],f[N];
inline void in(int &x){
	char ch=getchar(); 
	for (;ch<'0'||ch>'9';ch=getchar()); 
	for (x=0;ch>='0'&&ch<='9';ch=getchar())
		x=x*10+ch-48; 
}
void add(int x,int y,int z) {
	e[++tot]=E(head[x],y,z); 
	head[x]=tot,enter[y]++;
}
void build(int &k,int l,int r) {
	k=++trtot; int mid=(l+r)>>1;
	if (l==r) {
		pos[l]=k; return; 
	}
	build(ls[k],l,mid),build(rs[k],mid+1,r);
	add(ls[k],k,0),add(rs[k],k,0);
}
void ins(int k,int l,int r,int x,int y) {
	if (l==x&&r==y) {
		add(k,trtot,1); return;
	}
	int mid=(l+r)>>1;
	if (y<=mid) 
		ins(ls[k],l,mid,x,y); 
	else
		if (x>mid) 
			ins(rs[k],mid+1,r,x,y); 
		else { 
			ins(ls[k],l,mid,x,mid); 
			ins(rs[k],mid+1,r,mid+1,y);
		}
}
int main() {
	in(n),in(cnt),in(m);int i,x;
	build(rt,1,n);
	for (i=1; i<=cnt; i++)
		in(x),in(a[pos[x]]);
	int l,r,t;
	for (i=1; i<=m;i++){
		in(l),l--,in(r),in(t);
		trtot++;
		while (t--) {
			in(x),add(trtot,pos[x],0);
			if (l+1<x) ins(rt,1,n,l+1,x-1); 
			l=x;
		}
		if (x<r) ins(rt,1,n,x+1,r);
	}
	int first=0,tail=0;
	for (i=1; i<=trtot; i++) 
		if (!enter[i])
			h[++tail]=i,f[i]=1;
	while (first<tail) {
		x=h[++first];
		if (f[x]>inf) {
			puts("NIE"); return 0;
		}
		if (f[x]>a[x] && a[x]) {
			puts("NIE"); return 0;
		} else 
			f[x]=max(f[x],a[x]);
		for (i=head[x];i;i=e[i].next) {
			t=e[i].to,f[t]=max(f[t],f[x]+e[i].v);
			enter[t]--; 
			if (!enter[t]) h[++tail]=t;
		}
	}
	if (tail<trtot) {
		puts("NIE"); return 0;
	}
	puts("TAK");
	for (i=1;i<=n;i++) printf("%d ",f[pos[i]]);
	return 0;
}
```