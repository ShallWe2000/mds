---
title: 【bzoj 3107:】[cqoi2013]二进制a+b
date: 2017/01/01 22:00:00
tags:
  - 位运算
  - 构造
categories: 
  - 数学题

---

> 一个可以DP但是构造更优越的题？

<!--more--> 

# 题目
<center><h2>3107: [cqoi2013]二进制a+b</h2>

Time Limit: 10 Sec  Memory Limit: 128 MB</center>

## Description

输入三个整数$a$,$b$,$c$，把它们写成无前导$0$的二进制整数。比如$a=7$, $b=6$,$c=9$，写成二进制为$a=111$,$b=110$,$c=1001$。接下来以位数最多的为基准，其他整数在前面添加前导$0$，使得$a$,$b$,$c$拥有相同的位数。比如在刚才的例子中，添加完前导$0$后为$a=0111$,$b=0110$,$c=1001$。最后，把$a$, $b$, $c$的各位进行重排，得到$a’$, $b’$,$c’$，使得$a’+b’=c’$。比如在刚才的例子中，可以这样重排：$a’=0111$,$b’=0011$,$c’=1010$。
你的任务是让$c’$最小。如果无解，输出$-1$。

## Input
输入仅一行，包含三个整数$a,b,c$。

## Output 
输出仅一行，为$c’$的最小值。
## HINT
$a,b,c<=2^{30}$
解题报告
考场现场<strong>构造失败</strong>，$5$分钟修改$\to AC$;
之前观察<strong>二进制加法</strong>时得到消$1$技巧，就是在$11...1$末位加$1$,可以变化成$100...0$；
进行这样一个摆放操作，令原数为连续的$x$个$1$，消耗$1$个$1$,共使用$x+1$个$1$,可以认为收益是消除$x$个$1$,代价是$1<<x$;
现在我们把题意简化成$a$个$1+$$b$个$1$消去$c$个$1$的问题，令$a>=b$

 1. 当$a>=c$，根据上面提到的消$1$方法，可以在$a$中拿出$c$个$1$摆在连续的位置，然后在末位放置一个$1$,这样形成$1..00$共$c$个0的数段，$b$中剩余的$b-1$个$1$可以从后往前填充这个数段，$a$中剩余的$1$和$b-1$中剩余的$1$放到数段后方，这样构造出一个满足条件的$c'$，我们需要证明这个$c'$是最小的： 观察$c'$的形态 $100...011...1011...1$ 在$1$数不变的情况下，可以变小的可能是最靠后的$0$能否被前方的一个$1$取代，考虑这个$0$是消$1$过程中最末尾的$0$,如果把这个$0$提前，消$1$个数会减少,如果此时要保证消$1$个数不变，需要重启一段消$1$,那么新开启的一段消$1$开头多出的一个$1$会使结束的$0$位置后移，<strong>比原构造数劣</strong>;
 2. 当$a<c$,先将$a$个$1$消掉，形成$1..a0_s..$，然后在头部$1$处放置$b-1$中的$1$，每次可以消掉一个$1$,剩余的$b-1$中的$1$放在尾部，得到的构造数形式为$100...011...1$为最优；


时间复杂度$O(logn)$,有一种$O(log^4n)$的$DP$,见[xym的blog](http://blog.csdn.net/xym_CSDN/article/details/52301891)

# 代码
```c++
#include<iostream>
#include<cstdio>
using namespace std;
long long tmp;
long long mi[63];
int t,a,b,c,l;
inline int len(long long x) {
	for (int i=60;i>=0;i--) 
		if (x/mi[i]%2) return i+1;
	return 0;
}
inline int one(int x) {
	int tmp=0;
	while (x) {
		tmp+=x&1;
		x>>=1; 
	}
	return tmp;
}
inline long long calc(int x,int y,int z) {
	int cut=x+y-z;
	if (x==0&&cut!=0) return -1; 
	if (cut==0) return mi[x+y]-1;
	if (cut<0||cut>=x+y) return -1;
	if (cut<=y) {
		long long tmp=mi[cut],l=1; 
		for (int i=1;i<cut;i++) {
			if (l==x) break; 
			l++; tmp+=mi[i]; 
		}
		for (int i=l+1;i<=x;i++) 
			tmp=tmp*2+1; 
		for (int i=cut+1;i<=y;i++) 
			tmp=tmp*2+1; 
		return tmp;
	} else {
		long long tmp=mi[y];int cnt=1;
		for (int i=y+1;i<=cut;i++) 	
			tmp=tmp*2LL,cnt++;
		for (int i=1;i<=x-cnt;i++) 
			tmp+=mi[i]; 
		return tmp;
	}
}
int main() {
//	freopen("binary.in","r",stdin);
//	freopen("binary.out","w",stdout);
	mi[0]=1; 
	for (int i=1;i<=62;i++) 
		mi[i]=mi[i-1]*2LL;
	l=0;
	cin>>a>>b>>c; 
	l=max(l,len(a)); 
	l=max(l,len(b));
	l=max(l,len(c)); 
	a=one(a),b=one(b),c=one(c);
	if (a>b) swap(a,b);
	if ((tmp=calc(a,b,c))!=1&&len(tmp)<=l)
		cout<<tmp<<endl;
	else 
		cout<<-1<<endl; 
	return 0;
}
	
```