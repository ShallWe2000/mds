---
title: 【bzoj 4422】[Cerc2015]Cow Confinement
date: 2017/01/01 22:00:00
tags: 
  - 线段树
  - 扫描线
categories: 
  - 数据结构

---

> 一个很有难度的题呦。。。

<!--more-->

# 题目

<center><h2>4422: [Cerc2015]Cow Confinement</h2> </center>

## Description

一个$10^6$行$10^6$列的网格图，上面有一些牛、花和一些矩形围栏，围栏在格子的边界上，牛和花在格子里，牛只能向下或向右走，牛也不能穿过围栏和地图边界，求每头牛它能到达的花的数量。注意栅栏不会相交

## Input

第一行一个数f表示矩形围栏的数量。
接下来f行，每行四个数$x1,y1,x2,y2$，表示$(x1,y1)$在围栏内部矩形的左上角，$(x2,y2)$在右下角。
接下来一行一个数m表示花的数量。
接下来m行每行两个数$x,y$，表示在$(x,y)$处有一朵花。
接下来一行一个数$n$表示牛的数量。
接下来n行每行两个数$x,y$，表示在$(x,y)$处有一头牛。

## Output

总共n行，每行一个数ans，第i个数表示第i头牛能到ans个花。

<h1>解题报告</h1></br>
<p>考试的时候连暴力都没有打出来，真是十分的失败。</br>
介绍一点部分分的做法：</p>

 1. 如果$x,y$比较小的话，那么很显然可以直接暴力`dfs`
 2. 如果$n$比较小的话，可以离散化，然后暴力$O \left( n^3 \right) $
 3. 同样如果要做到$O \left( n^2 \right)$ ,就需要dp，也就是令$f[i][j]=f[i][j+1]+f[i+1][j]-f[i+1][j+1]$ ，简单易懂。
 
考虑优化上述dp,离线+线段树，有这样一个事情， 就是cows只能向下或者向右走，这样对于一个点$(x,y)$，他可能拥有$(x+1,y)$ $(x,y+1)$无法到达的一些flowers，用类似于差分的想法， 令$f[i]$表示当前行$f[i+1]$无法到达的花朵，实际上cow$(x,y)$的答案就是找到下方第一个栅栏（实际上是横向的一条边）$(x_b,y)$,查询当前列$y$一个差分的和$(x,x_b)$


问题更加单纯了，就是转移，详细一点说就是从$y+1$到$y$，差分信息怎样变化

 1. 没有栅栏，只有花？只需要**单点修改**下就好
 2. ![进入栅栏](http://img.blog.csdn.net/20160810180037807) 
 出现一个栅栏? ，也就是进入栅栏的边界，设$(x_l ,x_r)$是栅栏的上下坐标，很显然，$(x_l,x_r)$部分的差分应该删掉，并且标记为被覆盖，而$x_l-1$这个位置，会获得$(x_l,x_r)$部分的差分，这点很容易，难以想到的一点是，此时需要询问一下$cow(x_r+1,y+1)$,并记录这个数值，为3做准备。
  ![出栅栏](http://img.blog.csdn.net/20160810181056615)  
 这个首先也需要做一个区间归零，然后要在$(x_l-1,y)$的位置减去这个区间加入时做的询问，原因是没有了栅栏的限制，下方的差分会计算两次。这个我感受了好久。
 
# 程序的实现
单点修改，区间覆盖（归零），区间查询，从某位置开始第一个障碍点的添加和查询
 好难的样子，但好像一个线段树就艹掉了。
 
# 代码
 

```c++
#include<cstdio>
#include<iostream>
#include<cstring>
#include<cstdlib>
#include<algorithm>
using namespace std;
const int X=1000001,Y=1000000,N=2000001;
struct FEN{
	int xl,xr,y,i;
	bool flag;
	bool operator < (const FEN &o)const{
		return y!=o.y?y>o.y:xl<o.xl;
	}
} fen[N<<1];

struct FLO{
	int x,y;
	bool operator < (const FLO &o)const{
		return y>o.y;
	}
} flo[N];
struct CS{
	int x,y,i;
	bool operator < (const CS &o)const{
		return y>o.y;
	}
} cow[N];
struct SS{
	int nm;
	bool cover,cut;
} seg[X<<2];
int ans[N],fs[N];;
char * cp=(char *)malloc(20000000);
inline void in(int &x){
	for (;*cp<'0'||*cp>'9';cp++);
	for (x=0;*cp>='0'&&*cp<='9';cp++)
		x=x*10+*cp-'0';
}
inline void pushup(int x){
	seg[x].nm=seg[x<<1].nm+seg[x<<1|1].nm;
	seg[x].cut=seg[x<<1].cut|seg[x<<1|1].cut;
}
inline void paint(int x){
	seg[x].cover=1;
	seg[x].nm=0;
}
inline void pushdown(int x){
	if(seg[x].cover){
		paint(x<<1),paint(x<<1|1);
		seg[x].cover=0;
	}
}
void add(int x,int l,int r,int pur,int val){
	seg[x].nm+=val;
	if (l==r) return; 
	int mid=(l+r)>>1;
	pushdown(x);
	if(pur<=mid)
		add(x<<1,l,mid,pur,val);
	else 
		add(x<<1|1,mid+1,r,pur,val);
	pushup(x);
}
void cover(int x,int l,int r,int L,int R){
	if(L<=l&&r<=R){
		paint(x);
		return;
	}
	pushdown(x);
	int mid=(l+r)>>1;
	if(L<=mid)cover(x<<1,l,mid,L,R);
	if(R>mid)cover(x<<1|1,mid+1,r,L,R);
	pushup(x);
}
int query(int x,int l,int r,int L,int R){
	if(L<=l&&r<=R){
		return seg[x].nm;
	}
	int mid=(l+r)>>1,ans=0;
	pushdown(x);
	if(L<=mid)ans+=query(x<<1,l,mid,L,R);
	if(R>mid)ans+=query(x<<1|1,mid+1,r,L,R);
	return ans;
}
void update(int x,int l,int r,int pur){
	if(l==r){
		seg[x].cut^=1;
		return; 
	}
	int mid=l+r>>1;
	pushdown(x);
	if(pur<=mid)
		update(x<<1,l,mid,pur);
	else 
		update(x<<1|1,mid+1,r,pur);
	pushup(x);
}
int next(int x,int l,int r,int L){
	if(l>=L){
		if(seg[x].cut){
			while(l!=r)
				if(seg[x<<1].cut)
					x<<=1,r=l+r>>1;
				else 
					x=x<<1|1,l=(l+r>>1)+1;
			return l;
		}
		else return 0;
	}
	int tmp,mid=(l+r)>>1;
	pushdown(x);
	if(L<=mid&&(tmp=next(x<<1,l,mid,L)))
		return tmp;
	else 
		return next(x<<1|1,mid+1,r,L);
}
int main(){
//	freopen("4422.in","r",stdin);
	fread(cp,1,20000000,stdin);
	int f,m,n,x1,y1,x2,y2;
	in(f);
	for(int i=f;i--;){
		in(x1),in(y1),in(x2),in(y2);
		fen[i<<1]=(FEN){x1,x2,y1-1,i,0};
		fen[i<<1|1]=(FEN){x1,x2,y2,i,1};
	}
	sort(fen,fen+(f<<1));
	in(m);
	for(int i=m;i--;)
		in(flo[i].x),in(flo[i].y);
	sort(flo,flo+m);
	in(n);
	for(int i=0;i<n;++i){
		in(cow[i].x),in(cow[i].y);
		cow[i].i=i;
	}
	sort(cow,cow+n);
	f=m=n=0;
	update(1,1,Y,Y);
	int sum,cut;
	for(int i=Y;i;--i){
		for(;fen[f].y==i;++f)
			if(fen[f].flag==0){
				cover(1,1,Y,fen[f].xl,fen[f].xr);
				if(fen[f].xl!=1)
					add(1,1,Y,fen[f].xl-1,-fs[fen[f].i]);
				if(fen[f].xl!=1)
					update(1,1,Y,fen[f].xl-1);
				if(fen[f].xr!=Y)
					update(1,1,Y,fen[f].xr);
			}else{
				cut=next(1,1,Y,fen[f].xr);
				sum=query(1,1,Y,fen[f].xl,fen[f].xr);
				fs[fen[f].i]=query(1,1,Y,fen[f].xr+1,cut);
				cover(1,1,Y,fen[f].xl,fen[f].xr);
				if(fen[f].xl>1)
					add(1,1,Y,fen[f].xl-1,sum+fs[fen[f].i]);
				if(fen[f].xl!=1)
					update(1,1,Y,fen[f].xl-1);
				if(fen[f].xr!=Y)
					update(1,1,Y,fen[f].xr);
			}
		for(;flo[m].y==i;++m){
			add(1,1,Y,flo[m].x,1);
		}
		for(;cow[n].y==i;++n){
			cut=next(1,1,Y,cow[n].x);
			ans[cow[n].i]=query(1,1,Y,cow[n].x,cut);
		}
	}
	for(int i=0;i<n;++i)
		printf("%d\n",ans[i]);
	return 0;
}
```