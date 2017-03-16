---
title: 【UVA 12633】Super Rooks on Chessboard
date: 2017-02-14 14:26:44
tags:
  - FFT
categories:
  - 数学题
---

> 好像没怎么在UVA上做过题。。

<!--more-->

# 题目大意

每个黑点可以使所在的行列左上右下对角线不可用...问剩下多少点可用. 


# 解题报告

首先, 可以光看行列, 得到剩下哪些行列可用. 因为左上右下对角线可以通过 $x+y$ 表示一发, 所以记录每个 $x+y$ 有多少暂时合法的行列组合. 

因为$f[x+y] = \sum f[x]*f[y]$ , 所以用FFT大力跑一些, 然后把合法的$x+y$对应的方案数求和, 就得到答案了. 

# 代码

```c++

#include <bits/stdc++.h> 

using namespace std;

#define FORU(i, a, b) for (int i = int(a), nn = int(b); i <= nn; ++i) 
#define FORD(i, a, b) for (int i = int(a), nn = int(b); i >= nn; --i) 
#define REP(i, b) for (int i = 0, nn = int(b); i < b; ++i) 
#define DEBUG(x) cout << (#x)  << ' ' << x << endl

typedef long long ll; 
typedef double ff; 

const int N = 200100; 
const ff pi = acos(-1); 

struct cmx {
 	ff x, y; 
 	cmx(ff x = 0, ff y = 0) 
 		:x(x), y(y) {}
 	cmx operator + (const cmx &b) const {
 	 	return cmx(x + b.x, y + b.y); 
 	}
 	cmx operator - (const cmx &b) const { 
 		return cmx(x - b.x, y - b.y); 
 	}
 	cmx operator * (const cmx &b) const { 
 		return cmx(x*b.x - y*b.y, x*b.y + y*b.x);
 	} 
} A[N], B[N]; 
 
 		 
int r, c, n, re[N], T, _n; 
bool ro[N], co[N], xy[N]; 

inline void in(int &x) { 
	char ch = getchar(); 
	for (;ch < '0' || ch > '9'; ch = getchar()); 
	for (x = 0; ch >= '0' && ch <= '9'; ch = getchar())
	 	x = x * 10 + ch - 48;
} 

void fft(cmx *a, int f = 1) { 
	for (int i = 0; i < _n; ++i) 
		if (i < re[i]) swap(a[i], a[re[i]]); 
	
	for (int m = 1; m < _n; m <<= 1) {
		cmx wn = cmx(cos(pi / m), f * sin(pi/m)); 
		for (int i = 0; i < _n; i += m << 1) { 
			cmx w = cmx(1, 0); 
			for (int j = 0; j < m; ++j) { 
				cmx x = a[i+j], y = a[i+j+m]*w; 
				a[i+j]=x+y, a[i+j+m]=x-y; 
				w = w * wn; 
			} 
		} 
	} 
	if (f == -1) 
		for (int i = 0; i < _n; ++i) 
			a[i].x /= (ff)_n; 
} 

			
void solve(int test) { 
	memset(A, 0, sizeof(A)); 
	memset(B, 0, sizeof(B)); 
	
	for (int i = 0; i < r; ++i) 
		if (!ro[i]) A[i].x = 1;// cout << i << ' ';
//	cout << endl;  
	for (int i = 0; i < c; ++i) 	
		if (!co[i]) B[i].x = 1;// cout << i << ' ';
//	cout << endl; 
	 
	
	for (_n = 1; _n < (r+c); _n <<= 1); 
	for (int i=0, j=0; i < _n; ++i) { 
		re[i] = j; 
		for (int k=_n>>1; (j ^= k) < k; k>>=1);
	} 
	
//	DEBUG(_n); 
	fft(A, 1), fft(B, 1); 
	for (int i=0; i < _n; ++i) 
		A[i] = A[i] * B[i]; 
	fft(A, -1); 
	
	ll ans = 0; 
	for (int i=0; i < _n; ++i) 
		if (!xy[i]) ans += (ll)(A[i].x + 0.5); 
	
	printf("Case %d: %lld\n", test, ans); 
} 

			
int main() {
	in(T); 
	FORU(I, 1, T) {
 		memset(ro, 0, sizeof(ro)); 
 		memset(co, 0, sizeof(co)); 
 		memset(xy, 0, sizeof(xy)); 
 		in(r), in(c), in(n); 	
 		int x, y; 

 		REP(i, n) { 
 			in(x), in(y), --x, y = c-y; 
 			ro[x] = 1, co[y] = 1, xy[x+y] = 1; 
 		}
 		
 		solve(I); 
 	} 
 	return 0; 
} 
```
