---
title: 【bzoj 3956】Count
date: 2017/01/01 22:00:00
tags: 
  - 倍增
  - 单调栈
  - 树状数组
categories: 
  - 数据结构

---

> 一个送分题？

<!--more-->

# 题目
<center><h2>3956: Count</h2>
Time Limit: 10 Sec  Memory Limit: 512 MB</center>

## Description

![图片标题](http://www.lydsy.com/JudgeOnline/upload/201504/11%284%29.png)
<!--more-->
## Input

![图片标题](http://www.lydsy.com/JudgeOnline/upload/201504/22%281%29.png)

## Output

![图片标题](http://www.lydsy.com/JudgeOnline/upload/201504/22%281%29.png)

## HINT
$M,N<=3*10^5,A_i<=10^9$

# 解题报告
比较简单的一道题，TA给的送分题。

重要结论：满足条件的数对不会超过$2n$个

证明：假定$i$是合法数对$(i,j)$中比较小的，那么$j$只能是$i$右侧第一个大于$i$的数，同理假定$j$是合法数对$(i,j)$中较小的，$i$只能是$j$左侧第一个大于$j$的数，所以，一个数作为较小数，最多只能向左向右产生一个合法数对，所以，满足条件的数对不会超过$2n$个

寻找合法数对的做法:使用单调栈，维护一个单调变小的数列，每次弹栈和入栈的时候记录数对的左右端点位置,这个做法正确性显然。

统计答案做法:查询区间$(l,r)$最大值的位置$x$,可以得知，在这个区间中的合法数对一定不会跨越$x$,那么我们可以确定，起点在$(l,x-1)$的合法数对$+$终点在$(x+1,r)$的合法数对就是答案！

# 代码

```c++
#include<iostream>
#include<cstdio>
#include<cstring>
#include<algorithm>
#include<cmath>
using namespace std;
const int N=1000001; 
const int K=21;
int sta[N],top,n,f[N][K],type;
int be[N<<1],en[N<<1],a[N],k,m,lg[N];
char * cp=(char *)malloc(30000000);
inline void in(int &x){
	x=strtol(cp,&cp,10);
}
char * os=(char *)malloc(20000000),*op=os;
void ou(int x){
	if(x){
		ou(x/10);
		*op++='0'+x%10;
	}
}
inline void out(int x){
	if(x)ou(x);
	else *op++='0';
	*op++='\n';
}
inline void add(int a[],int x,int val) {
	for (;x<=n;x+=x&-x) a[x]+=val; 
}
inline int query(int a[],int x) {
	if (!x) return 0; int tmp=0; 
	for (;x;x-=x&-x) tmp+=a[x]; 
	return tmp; 
}
inline int st(int x,int y) {
	int l=lg[y-x+1],one=f[x][l],two=f[y-(1<<l)+1][l]; 
	if (a[one]>a[two]) return one; 
	else return two; 
}
	
int main() {
//	freopen("cp.in","r",stdin); 
//	freopen("cp.out","w",stdout);
	fread(cp,1,30000000,stdin);
	in(n),k=log2(n),in(m),in(type),lg[1]=0;
	for (int i=2;i<=n;i++) lg[i]=lg[i>>1]+1; 
	int x,y,j,last=0; top=0;
	for (int i=1;i<=n;i++) {
		in(a[i]);
		while (a[i]>a[sta[top]]&&top) {
			j=sta[top],--top; 
			add(be,j,1),add(en,i,1); 
		}
		if (top) {
			add(be,sta[top],1),add(en,i,1);
			if (a[sta[top]]==a[i]) --top;
		}
		sta[++top]=i; 
	}
	for (int i=1;i<=n;i++) f[i][0]=i;
	for (int i=1;i<=k;i++) 
		for (int j=1;j<=n;j++)
			if (j+(1<<i-1)>n) f[j][i]=f[j][i-1];
			else 
				if (a[f[j][i-1]]>a[f[j+(1<<i-1)][i-1]])
					f[j][i]=f[j][i-1]; 
				else
				 	f[j][i]=f[j+(1<<i-1)][i-1]; 
	while (m--) {
		in(x),in(y); 
		if (type) {
			x=(x+last-1)%n+1,y=(y+last-1)%n+1; 
			if (x>y) swap(x,y); 
		}
		j=st(x,y); 
		last=query(be,j-1)-query(be,x-1);
		last=last+query(en,y)-query(en,j);
		out(last);
	}
	fwrite(os,1,op-os,stdout);
	return 0; 
}	
```
