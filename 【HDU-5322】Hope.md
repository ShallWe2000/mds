---
title: 【HDU-5322】Hope
date: 2017-02-14 14:30:09
tags:
  - FFT
  - CDQ分治
  - 递推
categories:
  - 数学题
---

> 明明是一个代码超短，啥高级东西都不需要的化式子题，硬生生地上了分治FFT

<!--more-->

# 题目大意
对于一个排列$P$ , 如果存在最小的$j, j>i, p[j]>p[i]$ , 那么从$i$向$j$连一条边, 每个联通分量对答案的贡献是$sz^2$ , 问长度为$n$的所有排列的总贡献 。 


# 解题报告

从$1$到$n$枚举$i$, 枚举$i$插入的位置， 显然$i$会将排列形成的连通图分成前后两部分, 所以可以得到: 

$$f[i] = \sum_{j=1}^i C_{i-1}^{j-1} * (j-1)! * j^2 * f[i-j]$$ 

再化一次,

$$f[i] = (i-1)! \sum_{j=1}^{i}j^2 * \frac{f[i-j]}{(i-j)！}$$

这样就可以分治fft做了, 但是好像可以继续化简.  

$$f[i] = (i-1)! * \sum_{k=1}^{i} \frac{f[k]}{k!} * (i-k)^2$$
$$= (i-1)! * \sum_{k=1}^{i} \frac{f[k]}{k!} * (i^2 - 2k * i +k^2)$$

令 $A[k]=\frac{f[k]}{k!}, B[k] = \frac{f[k]*k}{k!}, C[k] = \frac{f[k]*k^2}{k!}$ , 就可以$O(1)$得到$f[i]$ .

显然， $A, B, C$可以通过$O(1)$的递推得到。

所以复杂度是$O(n)$ . 

# 代码

* FFT

```c++
#include <bits/stdc++.h> 

using namespace std;

#define FORU(i, a, b) for (int i = int(a), nn = int(b); i <= nn; ++i) 
#define FORD(i, a, b) for (int i = int(a), nn = int(b); i >= nn; --i) 
#define REP(i, b) for (int i = 0, nn = int(b); i < b; ++i) 

typedef long long ll; 
typedef double ff; 

const int N=100001,p=998244353;
int dp[N], fac[N], inv[N], n;

inline int fast(int x, int k){
	int ans = 1;
	for(; k; k >>= 1, x = 1LL*x*x%p)
		if (k & 1) ans = 1LL*ans*x%p;
	return ans;
}

void solve(){
	fac[0] = inv[0]=1;
	for(int i=1; i < N; ++i)
		fac[i] = 1LL * fac[i-1] * i % p;
	inv[N-1] = fast(fac[N-1], p-2); 
	for (int i=N-2; i; --i) 
		inv[i] = 1LL * inv[i+1] * (i+1) % p; 
	int sum1=0, sum2=0, sum3=0;
	
	for(int i = 0; i<N; ++i){
		if(!i)	dp[i]=1;
		else dp[i]=1LL*fac[i-1]*((1LL*i*i%p*sum1%p - 2ll*i*sum2%p + sum3)%p+p)%p;
		sum1 = (sum1 + 1LL*dp[i]*inv[i])%p;
		sum2 = (sum2 + 1LL*i*dp[i]%p*inv[i])%p;
		sum3 = (sum3 + 1LL*i*i%p*dp[i]%p*inv[i])%p;
	}
}

int main(){
	solve();
	while(~scanf("%d",&n))
		printf("%d\n",dp[n]);
	return 0;
}
```

* 更好的做法

```c++
#include <bits/stdc++.h> 

using namespace std;

#define FORU(i, a, b) for (int i = int(a), nn = int(b); i <= nn; ++i) 
#define FORD(i, a, b) for (int i = int(a), nn = int(b); i >= nn; --i) 
#define REP(i, b) for (int i = 0, nn = int(b); i < b; ++i) 

typedef long long ll; 
typedef double ff; 

const int N=100001,p=998244353;
int dp[N], fac[N], inv[N], n;

inline int fast(int x, int k){
	int ans = 1;
	for(; k; k >>= 1, x = 1LL*x*x%p)
		if (k & 1) ans = 1LL*ans*x%p;
	return ans;
}

void solve(){
	fac[0] = inv[0]=1;
	for(int i=1; i < N; ++i)
		fac[i] = 1LL * fac[i-1] * i % p;
	inv[N-1] = fast(fac[N-1], p-2); 
	for (int i=N-2; i; --i) 
		inv[i] = 1LL * inv[i+1] * (i+1) % p; 
	int sum1=0, sum2=0, sum3=0;
	
	for(int i = 0; i<N; ++i){
		if(!i)	dp[i]=1;
		else dp[i]=1LL*fac[i-1]*((1LL*i*i%p*sum1%p - 2ll*i*sum2%p + sum3)%p+p)%p;
		sum1 = (sum1 + 1LL*dp[i]*inv[i])%p;
		sum2 = (sum2 + 1LL*i*dp[i]%p*inv[i])%p;
		sum3 = (sum3 + 1LL*i*i%p*dp[i]%p*inv[i])%p;
	}
}

int main(){
	solve();
	while(~scanf("%d",&n))
		printf("%d\n",dp[n]);
	return 0;
}
```