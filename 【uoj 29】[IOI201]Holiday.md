---
title: 【uoj 29】[IOI201]Holiday
date: 2017/01/01 22:00:00
tags: 
  - 主席树
  - 分治
categories: 
  - 数据结构

---
# 题目
<center><h2>29. 【IOI2014】Holiday</h2>
时间限制：$5s$ ,空间限制：$64$**MB** </center>

## 题目描述

健佳正在制定下个假期去台湾的游玩计划。在这个假期，健佳将会在城市之间奔波，并且参观这些城市的景点。
在台湾共有$n$个城市，它们全部位于一条高速公路上。这些城市连续地编号为$0$到 $n−1$。对于城市$i$ $(0<i<n−1)$而言，与其相邻的城市是$i−1$和$i+1$。但是对于城市$0$，唯一与其相邻的是城市$1$。而对于城市$n−1$，唯一与其相邻的是城市$n−2$。
<!--more-->
每个城市都有若干景点。健佳有$d$天假期并且打算要参观尽量多的景点。健佳已经选择了假期开始要到访的第一个城市。在假期的每一天，健佳可以选择去一个相邻的城市，或者参观所在城市的所有景点，但是不能同时进行。即使健佳在同一个城市停留多次，他也不会去重复参观该城市的景点。请帮助健佳策划这个假期，以便能让他参观尽可能多的景点。

## 任务

请实现函数`findMaxAttraction`，以计算健佳最多可以参观多少个景点。
`findMaxAttraction(n, start, d, attraction)`
$n$: 城市数。
$start$: 起点城市的编号。
$d$: 假期的天数。
$attraction$: 长度为$n$的数组；`attraction[i]` 表示城市$i$的景点数目，其中 $0≤i≤n−1$.
该函数应返回健佳最多可以参观的景点数。
$0≤d≤2n+n/2$ 每个城市中的景点数都是非负整数。$2≤n≤100000$

## 接口信息

`long long int findMaxAttraction(int n, int start, int d, int attraction[]);`

# 解题报告

**第一道** AC的 **交互题**，EXCITING!

这道题$fye$学姐互测中出过，当时写的是$O(n^2logn)$的暴力：

 - 分别预处理出`f[i]` `g[i]` `ff[i]` `gg[i]`, 分别表示向左走单次花费代价$i$的最大收益，向右走单向，向左走双向，向右走双向；
 - 这四个数组可以先枚举代价，再枚举最远点，在主席树上选择前$k$大的和；

需要优化这个算法：

 - 考虑**决策点的单调性**：以`g[i]`为例，如果`g[i]`的决策点为$x$，则`g[i+1]`的决策点$y>=x$，这个可以通过打表找规律/贪心思路得到，归纳下就可以得到决策点的单调性； 
 - 有了决策单调性，优化预处理的过程，就是很常见的`work(l,r,L,R)` $\to$ `work(l,mid-1,L,MID),work(mid+1,r,MID,R);` 这样 复杂度就是$O(n\log^2 n)$

调试过程中出现了极大的问题，最后竟然是**离散化**残了，不是很理解。
# 代码
四个数组一个`work`，因为**传参数量**影响了效率；

```c++
#include<iostream>
#include<cstring>
#include<algorithm>
#include<cstdio>
#include<cstdlib>
#include"holiday.h"
using namespace std;
const int N=100005;
struct pret{
	 int l,r,sz; long long sm;
} node[2000000];
long long int f[N*3],g[N*3],ff[N*3],gg[N*3];
long long ret;
int cnt,dis[N],head[N],b[N],tn,id[N],maxd;
bool flag=0;
void update(int from,int &x,int l,int r,int w,long long val) {
	x=++cnt, node[x]=node[from],node[x].sm+=val,node[x].sz+=1; 
	if (l==r) return; int mid=(l+r)>>1; 
	if (w<=mid) update(node[from].l,node[x].l,l,mid,w,val); 
	else update(node[from].r,node[x].r,mid+1,r,w,val); 
}
void query(int lson,int rson,int l,int r,int k) {
	if (k<=0) return;
	if (l==r) {
		ret+=1LL*b[l]*min(node[rson].sz-node[lson].sz,k);
		return; 
	}	
	int mid=(l+r)>>1,tmp;
	if ((tmp=node[node[rson].r].sz-node[node[lson].r].sz)>=k) 
		query(node[lson].r,node[rson].r,mid+1,r,k); 
	else ret+=node[node[rson].r].sm-node[node[lson].r].sm,
		query(node[lson].l,node[rson].l,l,mid,k-tmp);
}
void work(int l,int r,int low,int up,long long *f,bool type,int from) {
	if (l>r) return; if (low>up) return;
	int mid=(l+r)>>1,i,cho=0; 
	if (type) {
		for (i=low;i<=up;++i) {
			ret=0,query(head[from-1],head[i],1,tn,mid-dis[i]);
			if (ret>f[mid]||!cho) 
				f[mid]=ret, cho=i; 
		}
		work(l,mid-1,low,cho,f,type,from);
		work(mid+1,r,cho,up,f,type,from); 
	} else {
		for (i=up;i>=low;--i) {
			ret=0,query(head[i-1],head[from],1,tn,mid-dis[i]);
			if (ret>f[mid]||!cho) {
				f[mid]=ret, cho=i; 
			}
		}
		work(l,mid-1,cho,up,f,type,from);
		work(mid+1,r,low,cho,f,type,from); 
	}

}	
long long int findMaxAttraction(int n,int s,int d,int *a) {
	int i; long long ans=0; maxd=d; s++;	
	for (i=n;i;--i) a[i]=a[i-1];
	for (i=1;i<=n;++i) b[i]=a[i],id[i]=i; 
	sort(b+1,b+n+1); tn=unique(b+1,b+n+1)-b-1;
	for (i=1;i<=n;++i) id[i]=lower_bound(b+1,b+tn+1,a[i])-b; 
	for (i=1;i<=n;++i) 
		update(head[i-1],head[i],1,tn,id[i],a[i]); 
	for (i=s;i<=n;++i) dis[i]=i-s; 
	work(1,d,s,min(n,s+d),g,1,s); //s->n,dan
	for (i=s;i;--i) dis[i]=s-i; 
	work(1,d,max(1,s-d),s-1,f,0,s-1); //s-1->1 dan
	for (i=s+1;i<=n;++i) dis[i]=(i-s)<<1; 
	work(1,d,s,min(n,s+(d>>1)),gg,1,s);//s->n shuang
	for (i=s-1;i;--i) dis[i]=(s-i)<<1; 
	work(1,d,max(1,s-(d>>1)),s-1,ff,0,s-1);//s-1->1 shuang
	for (i=0;i<=d;++i) {
		ans=max(ans,max(gg[i]+f[d-i],ff[i]+g[d-i]));
	}
	return ans; 
}
```

 