---
title: 【bzoj 4423】[AMPPZ2013]Bytehattan
date: 2017/01/01 23:00:00
tags: 
  - 并查集
  - 平面图
categories: 
  - 图论题

---

> 利用到了网格图的性质

<!--more-->

# 题目
<center><h2>4423: [AMPPZ2013]Bytehattan</h2>

Time Limit: 3 Sec  Memory Limit: 128 MB</center>
## Description

比特哈顿镇有$n*n$个格点，形成了一个网格图。一开始整张图是完整的。
有$k$次操作，每次会删掉图中的一条边$(u,v)$，你需要回答在删除这条边之后$u$和$v$是否仍然连通。

## Input

第一行包含两个正整数$n,k$($2<=n<=1500,1<=k<=2n(n-1)$)，表示网格图的大小以及操作的个数。
接下来$k$行，每行包含两条信息，每条信息包含两个正整数$a,b$($1<=a,b<=n$)以及一个字符$c$($c=N$或者$E$)。
如果$c=N$，表示删除$(a,b)$到$(a,b+1)$这条边；如果$c=E$，表示删除$(a,b)$到$(a+1,b)$这条边。
数据进行了加密，对于每个操作，如果上一个询问回答为TAK或者这是第一个操作，那么只考虑第一条信息，否则只考虑第二条信息。
数据保证每条边最多被删除一次。
## Output
输出$k$行，对于每个询问，如果仍然连通，输出TAK，否则输出NIE。
# 解题报告
网格图？考虑对偶图的性质，当前一条边，如果它所隔开的两个格子在对偶图中已经联通，说明两点之间的边如果删除，就会形成割，两点也就不连通了。
![图片标题](https://leanote.com/api/file/getImage?fileId=57db61c4ab6441695200d2b2)
所以并查集搞一下
# 代码
```c++
#include<iostream>
#include<cstdio>
#include<cstring>
#include<algorithm>
#include<cstdlib>
using namespace std;
const int N=1501; 
int f[N*N],lastans,n,k,id[N][N];
inline void in(int &x) {
	char ch=getchar(); 
	for (;ch<'0'||ch>'9';ch=getchar());
	for (x=0;ch>='0'&&ch<='9';ch=getchar())
		x=x*10+ch-48; 
}
inline int num(int x,int y) {
	if (x<=0||y<=0||x>=n||y>=n) return 0; 
	return id[x][y]; 
}
inline int find(int x) {
	if (f[x]!=x) f[x]=find(f[x]); 
	return f[x]; 
}
inline void work(int a,int b,char type) {
	int one,two,x,y,fx,fy; 
	if (type=='N') {
		x=a-1,y=b,one=num(x,y); 
		x=a,y=b,two=num(x,y); 
	} else {
		x=a,y=b-1,one=num(x,y); 
		x=a,y=b,two=num(x,y); 
	}
	if ((fx=find(one))==(fy=find(two)))
		lastans=0;
	else 
		lastans=1,f[fx]=fy; 
	if (lastans) printf("TAK\n"); 
	else printf("NIE\n");
}
int main() {
//	freopen("4423.in","r",stdin); 
//	freopen("4423(1).out","w",stdout);
	in(n),in(k); int x=(n-1)*(n-1),y,x1,y1;
	char s1[10],s2[10]; 
	for (int i=0;i<=x;i++) f[i]=i; 
	int tmp=0;
	for (int i=1;i<n;i++) 
		for (int j=1;j<n;j++) 
			id[i][j]=++tmp;
	lastans=1;
	for (int i=1;i<=k;i++) {
		in(x),in(y);scanf("%s",s1);
		in(x1),in(y1);scanf("%s",s2); 
		if (!lastans) 
			swap(x,x1),swap(y,y1),swap(s1,s2); 
		work(x,y,s1[0]); 
	}
	return 0; 
}
```