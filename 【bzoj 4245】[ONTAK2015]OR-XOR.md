---
title: 【bzoj 4245】[ONTAK2015]OR-XOR
date: 2017/01/01 22:00:00
tags: 
  - 贪心
  - 链表
  - 位运算
categories: 
  - 数学题

---

> 花式刷水题。

<!--more-->
# 题目
<center><h4>4245: [ONTAK2015]OR-XOR</h4>
Time Limit: 10 Sec  Memory Limit: 256 MB
Submit: 424  Solved: 232</center>
## Description
给定一个长度为$n$的序列$a[1],a[2],...,a[n]$，请将它划分为$m$段连续的区间，设第$i$段的费用$c[i]$为该段内所有数字的异或和，则总费用为$c[1] \ or \ c[2] or ... or \  c[m]$。请求出总费用的最小值。

## Input
第一行包含两个正整数$n,m(1<=m<=n<=500000)$，分别表示序列的长度和需要划分的段数。
第一行包含$n$个整数，其中第$i$个数为$a[i](0<=a[i]<=10^{18})$。
## Output
输出一个整数，即总费用的最小值。
# 解题报告
看yveh弄了道题，TA说是水题，我就折腾了下，确实比较容易；
位运算，$or$的存在，说明分成的所有部分中，只要有一个部分的某个数位为$1$，那最终结果这一位就是$1$，所以需要让高位尽可能为$0$;
显然，如果某一位一共有奇数个，那这一位一定贡献答案；
偶数个的情况下，通过**异或的差分性质**，考虑可以让分成的部分某一位都为0的断点位置；
从高位向下进行，如果当前位满足条件的断点位置$>=m$个，且最后一个位置可以是断点，那就将当前位不满足条件的端点删除；否则就在答案中当前位，置$1$； 
算是一个贪心吧，我用链表维护了一下，跑得很慢； 
# 代码
```c++
#include<iostream>
#include<cstdio>
#include<cstring>
#include<algorithm>
const int N=500002;
using namespace std;
int next[N],pre[N];
long long n,m; 
long long a[N],ans=0; 
char *cp=(char *)malloc(10000000);
inline void in(long long &x) {
	for (;*cp<'0'||*cp>'9';cp++); 
	for (x=0;*cp>='0'&&*cp<='9';cp++) 
		x=x*10+*cp-48;
}
int main() {
//	freopen("a.in","r",stdin); 
//	freopen("a.out","w",stdout);
	fread(cp,1,10000000,stdin);
	in(n),in(m);
	for (int i=1;i<=n;i++) 
		in(a[i]);
	for (int i=0;i<=n;i++) 
		next[i]=i+1; 
	for (int i=1;i<=n+1;i++) 
		pre[i]=i-1;
	for (int i=1;i<=n;i++) 
		a[i]=a[i-1]^a[i];
	int tmp;
	for (int i=59;i>=0;--i) {
		tmp=0;
		for (int j=next[0];j!=n+1;j=next[j])
			if (!(a[j]&(1LL<<i))) tmp=tmp+1;
		if (tmp>=m&&!(a[n]&(1LL<<i))) {
			for (int j=1;next[j]!=n+1;j=next[j]) 
				if (a[j]&(1LL<<i)) 
					next[pre[j]]=next[j],pre[next[j]]=pre[j];
		} else
			ans|=(1LL<<i);
	}
	return printf("%lld\n",ans),0; 
}
```