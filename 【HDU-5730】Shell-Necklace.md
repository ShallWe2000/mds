---
title: 【HDU 5730】Shell Necklace
date: 2017-02-14 14:28:34
tags:
  - FFT
  - CDQ分治
categories:
  - 数学题
---

> 做一道分治FFT的板子题。

<!--more--> 

# 题目大意

将 $n$ 个数的序列分成若干个部分, 已知连续的长度为 $i$ 的部分, 方案数是 $a[i]$ , 问总方案数. 


# 解题报告

设长度为 $n$ 的方案数为 $f[n]$ , 那么 $f[i] = \sum f[j] * a[i-j]$ , 

这个是个裸题, 分治FFT, 应该说是CDQ + FFT, 因为 $f_j,(j < i)$ 对 $i$ 有等价的贡献;

在处理区间 $[l, r]$ 时, 先递归处理$[l, mid]$, 得到 $f[l...mid]$ ;

令 $A(i) = a(i)$ , $B(i) = f[i+l]$ , 做长度为 $r-l+1$ 的卷积. 

实际上应该做 $r-l+1 + mid-l + 1$ 的卷积, 因为需要得到的$f[mid+1...r]$ 对应 $F[mid -l + 1...r-l+1]$ 这部分的结果不会受到多项式乘法溢出结果的影响. 

# 代码

```c++
#include <bits/stdc++.h> 

using namespace std;

#define FORU(i, a, b) for (register int i = int(a), nn = int(b); i <= nn; ++i) 
#define FORD(i, a, b) for (int i = int(a), nn = int(b); i >= nn; --i) 
#define REP(i, b) for (register int i = 0, nn = int(b); i < b; ++i) 
#define DEBUG(x) cout << (#x) << ": " << x << endl; 

typedef long long ll; 
typedef double ff; 

const int N = 200100; 
const ff pi = acos(-1); 
const int p = 313; 

struct cmx { 
	ff x, y; 
	cmx(ff x = 0, ff y = 0) 
		:x(x), y(y) {}
	cmx operator + (const cmx & b) const { 
		return cmx(x + b.x, y + b.y); 
	} 
	cmx operator - (const cmx & b) const { 
		return cmx(x - b.x, y - b.y); 
	} 
	cmx operator * (const cmx & b) const { 
		return cmx(x*b.x - y * b.y, x*b.y + y*b.x);
	}
} _[N], __[N]; 

int n, a[N], _n, dp[N], re[N]; 

inline void dft(cmx *a, int f = 1) {
 	for (register int i =0 ; i < _n; ++i) 
 		if (i < re[i]) swap(a[i], a[re[i]]); 
 	cmx x, y; 
 	for (register int m = 1; m < _n; m <<= 1) {
 	 	cmx wn = cmx(cos(pi/m), f * sin(pi/m)); 
 	 	for (register int i = 0; i < _n; i += m << 1) { 
 	 		cmx w = cmx(1, 0); 
 	 		for (register int j = 0; j < m; ++j) { 
 	 			x = a[i+j], y = w*a[i+j+m]; 
 	 			a[i+j]=x+y, a[i+j+m]=x-y; 
 	 			w = w * wn; 
 	 		} 
 	 	}
 	} 
 	
 	if (f == -1) { 
 		REP(i, _n) a[i].x /= (ff)_n; 
 	}
}

void solve(int l, int r) { 
	if (l == r) return; 
	int mid = (l + r) >> 1; 
	
	solve(l, mid); 
	
	int len = r - l + 1; 
	_n = 1; for (;_n <= (len); _n <<= 1); 
	for (register int i=0, j=0; i < _n; ++i) { 
		re[i] = j; 
		for (register int k = _n>>1; (j^=k) < k; k>>=1);
	} 
//	DEBUG(_n);
	for (register int i=0; i<_n; ++i) 
		_[i] = __[i] = cmx(0, 0); 
	for (register int i=l; i<=mid; ++i) 
		_[i-l] = cmx(dp[i], 0); 
	for (register int i=0; i < len; ++i) 
		__[i] = cmx(a[i], 0);
//	DEBUG(_[0].x);   DEBUG(__[0].x); 
	dft(_, 1), dft(__, 1); 
	REP(i, _n) _[i] = _[i] * __[i]; 
	dft(_, -1); 
//	DEBUG(_[1].x); 
	for (register int i = mid+1; i <= r; ++i) {
		dp[i] += (int)(_[i-l].x + 0.5); 
		dp[i] %= p; 
	} 
//	DEBUG(dp[1]); 
	solve(mid + 1, r); 
} 
	 
int main() { 
	while (scanf("%d", &n)) {
	 	if (!n) break; 
	 	memset(a, 0, sizeof(a)); 
	 	memset(dp, 0, sizeof(dp)); 
	 	FORU(i, 1, n) scanf("%d", &a[i]), a[i] %= p; 
	 	dp[0] = 1; 
	 	solve(0, n);
	 	printf("%d\n",  dp[n]); 
	} 
	return 0; 
}
```
