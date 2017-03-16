---
title: 【2012国家集训队test】calc(by WJMZBMR)
date: 2017/01/01 22:00:00
tags:
  - 倍增
  - 容斥原理
  - 拉格朗日差值

categories:
  - 数学题

---

> 好久不写blog了，屯了接近五十道题，准备好好写写；

<!--more-->

# 题目

## 题目描述

一个有序序列 $ a _1,...,a_n $ 是合法的，当且仅当 ：
- 长度为$n$；
- 其中的数两两不同；
- $a_i \in [1,A]$ 对所有的$i$成立;
一个序列的值$val(a)= \prod _{i=1} ^{i=n} a_i$;
求所有不同序列的值的和 $\pmod p$,$p$为素数.

## 数据范围
- $100 \% : A,p \leq 10^9, n \leq 500, p>A>n+1$

# 解题报告

## 第一个做法：

- 考虑数列中的数是$distinct$的，那么可以先规定一个递增的顺序，然后乘一个全排列；
- 观察数据范围，$A$很大，考虑**倍增**，令$f[A][n]$表示在$[1,A]$中选$n$个数；
	- 计算$f[2A][n]$:记从$[A+1,2A]$中选$n$个数的值为$b_i$，考虑在$[1,A]$中选$a$个数，在$[A+1,2A]$中选$n-a$个数;
	$$b\_i= \Sigma \_{j=0}^{j=i} A^{i-j} * a_j * C\_{i-j}^{A-j}$$
	- 如何快速求$C\_{i-j}^{A-j}$,化简得
	$$C\_{i-j}^{A-j}= \frac{(A-j)!}{(i-j)!(A-i)!}= \prod\_{k=i}^{k=j}(A-k) *(i-j)!^{-1}$$
    这就可以做了。
	- 然后？
	$$f[2A][n]=f[A][a]*b[n-a]$$
- 时间复杂度： $O(log A n^2)$ 可以$fft$优化下，$noip$后补一下；

---

## 第二个做法：

- 令$f[A][n]$表示$[1,A]$中选$n$个数的值和，写一个粗暴的转移：
$$f[A][n]=f[A-1][n-1]*A+f[A-1][n]$$
- 观察这个转移，$f[x][N]$是一个多项式的形式，$f[A-1][n]$相当于前缀和，$f[A-1][n-1] * A$是一个幂指数向左平移的操作，所以，对$f[x][N] \to f[x][N+1]$,多项式指数$+2$，那么$f[x][N]$指数为$2* n+1$,可以求出$f[i][n](i=1 \to 2n+1)$，然后进行**拉格朗日插值**；
$$f(x)=\sum_{i=0}^n a_ix^i=\sum_{i=0}^nf(x_i)\prod_{j=0}^n\frac{x-x_j}{x_i-x_j} [i≠j]$$

## 第三个做法：

- 容斥原理，来自$reflash$大爷，令$f[i]$表示$[1,A]$中选$i$个的方案数；
$$f[i] = g[1] * f[i - 1]  + \sum((-1)^{i - j + 1} * f[j] * C(i - 1,i - 1 - j) * (i - 1 - j)! * g[i - j])$$
- $g[n]=\Sigma_{j=1}^n j^n$，预处理伯努利数；
- 就是随便选-至少有一个重复的+至少两个重复的...;

# 代码

## 第一个做法：
```c++
#include<iostream>
#include<algorithm>
#include<cstring>
#include<cstdlib>
#include<cstdio>
using namespace std;
long long p,A,n,C[501][501],g[501],f[501],pow[501],jie[501],rej[501],ajie[501],rea[501];
inline long long rev(int x) {
	int ans=1;
	for (int k=p-2;k;k>>=1) {
		if ( k&1) ans=1LL*ans*x%p;
		x=1LL*x*x%p;
	}
	return ans;
}
int main() {
	scanf("%d%d%d",&A,&n,&p); long long i,j,k,pre;
	C[0][0]=1;
	for (i=1;i<=n;++i) {
		C[i][0]=1;
		for (j=1;j<=i;++j)
			C[i][j]=(C[i-1][j-1]+C[i-1][j])%p;
	}
	for (jie[0]=1,i=1;i<=n;++i) jie[i]=1LL*jie[i-1]*i%p;
	for (i=n-1,rej[n]=rev(jie[n]);i>=0;--i) rej[i]=1LL*(i+1)*rej[i+1]%p;
	for (j=30;j>=0;--j) if (A&(1<<j)) break;
	g[1]=1,g[0]=1;
	for (pre=1,pow[0]=1,--j;j>=0;--j) {
		memset(f,0,sizeof(f)); f[0]=1;
		for (i=1;i<=n;++i) pow[i]=1LL*pow[i-1]*pre%p;
		for (i=1,ajie[0]=1;i<=n&&i<=pre;++i)
			ajie[i]=1LL*ajie[i-1]%p*(pre-i+1)%p;
		for (i=n>pre?pre-1:n-1,rea[min(n,pre)]=rev(ajie[min(pre,n)]);i>=0;--i)
			rea[i]=1LL*rea[i+1]*max(1LL,pre-i)%p;  
		for (i=0;i<=n&&i<=pre;++i) {
			if (!g[i]) break;
			f[i]=0;
			for (k=i;k>=0;--k)
				(f[i]+=1LL*g[k]*pow[i-k]%p*rej[i-k]%p*rea[k]%p)%=p;
			f[i]=1LL*f[i]*ajie[i]%p;
		}
		for (i=n;i;--i)
			for (k=0;k<i;++k)
				(g[i]+=1LL*g[k]*f[i-k]%p)%=p;
		pre<<=1;
	 	if (A&(1<<j))
	 		for (i=n;i;--i)
	 			(g[i]+=1LL*g[i-1]*(pre+1)%p)%=p;
	 	pre+=((A>>j)&1);
	}
	printf("%d\n",1LL*g[n]*jie[n]%p);
	return 0;
}
```
## 第二种做法

```c++
#include<iostream>
#include<cstring>
#include<algorithm>
#include<cstdio>
#include<cstdlib>
#include<cmath>
using namespace std;
const int N=510;
int f[N<<2][N],n,p,A,ajie[N<<1],rea[N<<1],jie[N<<1],rej[N<<1],ans,gg[N<<1];
inline int rev(int x) {
	int ans=1;
	for (int k=p-2;k;k>>=1) {
		if (k&1) ans=1LL*ans*x%p;
		x=1LL*x*x%p;
	}
	return (ans+p)%p;
}
int main() {
	freopen("lu#4.in","r",stdin);
	freopen("lu#4.out","w",stdout);
	scanf("%d%d%d",&A,&n,&p);
	f[0][0]=1; int i,j,M,nn=n*2+1;
	for (i=1;i<=min(A,nn);++i) {
		f[i][0]=1;
		for (j=1;j<=n;++j) {
			f[i][j]=(f[i-1][j]+1LL*f[i-1][j-1]*i%p)%p;
		}
	}
	for (i=1,gg[0]=1;i<=nn;++i)
		gg[i]=1LL*gg[i-1]*i%p;
	if (A<=nn) {
		ans=1LL*f[A][n]*gg[n]%p;
		printf("%d\n",ans);
	}
	else {
		for (i=1,ajie[0]=1;i<=nn;++i)
			ajie[i]=1LL*ajie[i-1]*(A-i)%p;  
		for (i=1,jie[0]=1;i<=nn;++i)
			jie[i]=1LL*jie[i-1]*(p-i)%p;
		int tmp=rev(jie[nn-1]);
		for (i=1;i<=nn;++i) {
			(ans+=1LL*f[i][n]*rev(A-i)%p*tmp%p)%=p;
			tmp=1LL*tmp*(p-nn+i)%p*rev(i)%p;
		}
		ans=1LL*ans*ajie[nn]%p*gg[n]%p;
		printf("%d\n",ans);
	}
	return 0;
}		
```

## 第三种做法：
http://blog.csdn.net/qq_20669971/article/details/52790835

还有TA爷更科学的容斥：
http://blog.csdn.net/ta201314/article/details/52753481
