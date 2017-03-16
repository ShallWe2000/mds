---
title: 【poj 3241】Object Clustering
date: 2017/01/01 22:00:00
tags: 
  - 曼哈顿最小距离生成树
categories: 
  - 数据结构

---

> 曼哈顿距离最小生成树： 真·板子题

<!--more-->

# 题目
<center><h2>Object Clustering</h2></center>
## Description
We have $N$ ($N≤10000$) objects, and wish to classify them into several groups by judgement of their resemblance. To simply the model, each object has $2$ indexes $a$ and $b$ ($a,b≤500$). The resemblance of object $i$ and object $j$ is defined by $d_{i,j}=|a_i-a_j|+|b_i-b_j|$, and then we say $i$ is $d_{i,j}$resemble to $j$. Now we want to find the minimum value of $X$, so that we can classify the $N$ objects into $K$ ($K < N$) groups, and in each group, one object is at most $X$ resemble to another object in the same group, i.e, for every object $i$, if $i$ is not the only member of the group, then there exists one object $j$ ($i ≠ j$) in the same group that satisfies $d_{i,j} ≤ X$
<!--more-->
# Input
The first line contains two integers $N$ and $K$. The following $N$ lines each contain two integers $a$ and $b$, which describe a object.
# Output
A single line containsthe minimum $X$.
# 题目分析
简述题意：

模板题，学习了曼哈顿最小生成树，简单地描述一下， 就是对于一个点，只要连接八象限中曼哈顿距离最近的点，就一定可以构造出曼哈顿最小生成树。

考虑对于点$i,j$,$d_{i,j}=|a_i-a_j|+|b_i-b_j|$,令$a_i>a_j,b_i>b_j$,则距离$j$最近的点$i$具有$min(a_i+b_i)$
实际上我们<strong>不需要八个象限，只需要四个象限</strong>，利用位置关系的相对性就可以连全所需的边，可以通过<strong>坐标系的旋转</strong>,这样只需要一个函数就搞定了。

现在我们把问题放在$x$正半轴与直线$y=x$之间的区域，在这个区域内可以选择的点$j$满足$a_j-b_j > a_i-b_i$,只存在添加的后缀最大值？<strong>树状数组</strong>，把$x$坐标从大到小排个序，对$a-b$做一个<string>离散化</string>，上树状数组就搞定了（树状数组要反着写，就是把修改和查询的枚举顺序倒一下就好了）
# 代码
```c++
#include<iostream>
#include<cstdio>
#include<cstring>
#include<algorithm>
#include<cstdlib>
using namespace std;
const int N=10005;
const int inf=2000000000;
struct E{
	int x,y,v;
	E(int x=0,int y=0,int v=0)
		:x(x),y(y),v(v){}
} e[N<<2];
int bit[N],pos[N],id[N],vec[N];
int x[N],y[N],tot,cnt,n,k,f[N];
inline bool cmpi(int a,int b){
	return x[a]==x[b]?y[a]<y[b]:x[a]<x[b]; 
}
inline bool cmpe(E a,E b){
	return a.v<b.v; 
}
inline int dis(int a,int b){
	return abs(x[a]-x[b])+abs(y[a]-y[b]); 
}
inline void add(int x,int y,int dis){
	e[++tot]=E(x,y,dis); 
}
inline void update(int x,int val,int pur){
	for (;x;x-=x&-x)
		if (val<bit[x]) 
			bit[x]=val,pos[x]=pur; 
}
inline int query(int x){
	int tmp=-1,mn=inf;
	for (;x<=cnt;x+=x&-x) 
		if (bit[x]<mn)
			mn=bit[x],tmp=pos[x]; 
	return tmp;
}
inline int hash(int x){
	return lower_bound(vec+1,vec+1+cnt,x)-vec;
}
inline int find(int x){
	if (f[x]!=x) f[x]=find(f[x]); 
	return f[x];
}
void Mmst(){
	tot=0; 
	for (int dir=0;dir<4;dir++){
		if (dir==1||dir==3)
			for (int i=1;i<=n;i++) 
				swap(x[i],y[i]); 
		else
			for (int i=1;i<=n;i++)
				x[i]=-x[i]; 
		for (int i=1;i<=n;i++) 
			id[i]=i; 
		sort(id+1,id+1+n,cmpi);
		for (int i=1;i<=n;i++) 
			vec[i]=y[i]-x[i]; 
		sort(vec+1,vec+1+n); 
		cnt=unique(vec+1,vec+1+n)-vec; 
		for (int i=1;i<=n;++i) 
			bit[i]=inf,pos[i]=-1;
		int u,v;
		for (int i=n;i;i--){
			u=hash(y[id[i]]-x[id[i]]); 
			v=query(u);
			if (v!=-1)
				add(id[i],v,dis(id[i],v)); 
			update(u,x[id[i]]+y[id[i]],id[i]); 
		}
	}
}
		
int main(){
	while (scanf("%d%d",&n,&k)!=EOF){
		for (int i=1;i<=n;i++)
			scanf("%d%d",&x[i],&y[i]); 
		Mmst(); 
		sort(e+1,e+1+tot,cmpe);
		for (int i=1;i<=n;++i) f[i]=i; 
		k=n-k; 
		int a,b,x,y;
		for (int i=1;i<=tot;i++){
			a=e[i].x,b=e[i].y; 
			if ((x=find(a))!=(y=find(b))){
				--k;
				f[x]=y; 
				if (!k){
					printf("%d\n",e[i].v); 
					break; 
				}
			}
		}
	}
	return 0; 
}
```

