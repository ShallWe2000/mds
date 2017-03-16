---
title: 【bzoj 3166】[Heoi2013]Alo
date: 2017-02-14 14:09:33
tags:
  - STL
  - 可持久化数据结构
categories:
  - 数据结构
---


> 搞了一发可持久化trie树的板子题，一A舒爽

<!--more-->

# 题目大意
区间$[l, r]$的次大值为$k$, 则该区间的值为$max(k \space xor \space a_p | a_p ≠ k , l ≤ p ≤ r)$; 
求所有区间的值的最大值； 

# 解题报告

首先可以找到每个数$a_i$作为次大值的区间，去除被包含的小区间，可以发现，每个数最多用两个区间就可以表示他作为次大值的区间;


确定每个数的区间， 需要找到他左侧第一个比他大的数$l_1$, 第二个比他大的数$l_2$， 以及他右边第一个比他大的数$r_1$ 和第二个比他大的数$r_2$, 那么两个区间就是$(l_2, r_1)$ 和 $(l_1, r_2)$, 这个很显然； 


寻找$l_1, l_2, r_1, r_2$可以通过排序+`set`搞定； 


然后问题变成对于$k$求区间$[l, r]$中$max(k \space xor \space a_i | l <= i <= r, a_i !=k )$， 这个是可持久化trie树的裸题； 


# 代码

```cpp
#include <bits/stdc++.h> 
using namespace std; 
#define fi first
#define se second
#define mp make_pair
#define pb push_back
#define FORU(i, a, b) for (int i = a, nn = int(b); i <= nn; ++i) 
#define FORD(i, a, b) for (int i = a, nn = int(b); i >= nn; --i) 
#define REP(i, b) for (int i = 0, nn = int(b); i < nn; ++i) 
#define sz(x) (int)(x).size()
#define getTw(a, k) ( (a>>k) & 1)

typedef long long ll; 

const int N = 51000; 
struct node { 
	int son[2], sz; 
} t[N * 32];

int a[N], n, f[N][4], ans, tot, head[N];
// p[x][0] left_first, p[x][1] left_second, p[x][2] right_first, p[x][3] right_first;
int aAndp[N];  
set<int> place; 
char *cp = (char *)malloc(10000000); 
inline void in(int &x) { 
	for(; *cp < '0' || *cp > '9'; cp++); 
	for(x = 0; *cp >= '0' && *cp <= '9'; cp++) 
		x = x * 10 + *cp - 48; 
} 
inline bool cmp(int x, int y) { 
	return a[x] < a[y]; 
}
void build(int &x, int f, int a, int w) { 
	x = ++tot;
	if (f) t[x] = t[f], t[x].sz++;
	else t[x].sz = 1; 
	if (w >= 0) build(t[x].son[getTw(a, w)], t[f].son[getTw(a, w)], a, w-1); 
}
inline void build_trie() {
 	head[0] = 0; 
	FORU(i, 1, n) build(head[i], head[i-1], a[i], 30);
}
inline int trie_ans(int x, int l, int r) {
	int ans = 0, d, lson, rson; 
	FORD(k, 30, 0) {
		d = getTw(x, k);
		lson = t[l].son[d^1], rson = t[r].son[d^1];
 		if (t[rson].sz - t[lson].sz > 0) 
			ans |= 1<<k, l = lson, r = rson;
		else
			l = t[l].son[d], r = t[r].son[d]; 
	}
	return ans; 
} 	
inline void maxify(int &ans,int x) { 
	if (x > ans) ans = x; 
} 

int main() { 
	fread(cp, 1, 10000000, stdin); 	 
	in(n); 
	FORU(i, 1, n) in(a[i]), aAndp[i] = i;
	sort(aAndp + 1, aAndp + 1 + n, cmp);
	set<int> :: iterator ps;
	int x, y; 
	FORD(i, n, 1) { 
		x = aAndp[i]; 
		ps = place.upper_bound(x); 
		if (ps != place.end()) { 
			f[x][2] = *ps, ++ps;
			if (ps != place.end()) f[x][3] = *ps; 
			--ps; 
		} 
		if (ps != place.begin()) { 
			--ps, f[x][0] = *ps; 
			if (ps != place.begin()) 
				--ps, f[x][1] = *ps; 
		}
		place.insert(x); 
	} 
	build_trie(); 
	FORU(i, 1, n) { 
		if (!f[i][2]) f[i][2] = n + 1; 
		if (!f[i][3]) f[i][3] = n + 1; 	
		if (f[i][0]) maxify(ans, trie_ans(a[i], head[f[i][1]], head[f[i][2] - 1]) ); 
		if (f[i][2] != n+1) maxify(ans, trie_ans(a[i], head[f[i][0]], head[f[i][3] - 1]) ); 
	}
	cout << ans << endl; 
	return 0; 
} 
```


------------

# 原题

[题目链接](http://www.lydsy.com/JudgeOnline/problem.php?id=3166)
<body>
<title>Problem 3166. -- [Heoi2013]Alo</title><center><h2>3166: [Heoi2013]Alo</h2><span class="green">Time Limit: </span>20 Sec&nbsp;&nbsp;<span class="green">Memory Limit: </span>256 MB<br><span class="green">Submit: </span>874&nbsp;&nbsp;<span class="green">Solved: </span>416<br>[<a href="submitpage.php?id=3166">Submit</a>][<a href="problemstatus.php?id=3166">Status</a>][<a href="bbs.php?id=3166">Discuss</a>]</center><h2>Description</h2><div class="content"><p><span style="font-size: medium">Welcome to ALO ( Arithmetic and Logistic Online)。这是一个VR MMORPG&nbsp;，<br>
如名字所见，到处充满了数学的谜题。<br>
现在你拥有n颗宝石，每颗宝石有一个能量密度，记为ai，这些宝石的能量<br>
密度两两不同。现在你可以选取连续的一些宝石（必须多于一个）进行融合，设为 &nbsp;ai, ai+1, …, a&nbsp;j，则融合而成的宝石的能量密度为这些宝石中能量密度的次大值<br>
与其他任意一颗宝石的能量密度按位异或的值，即，设该段宝石能量密度次大值<br>
为k，则生成的宝石的能量密度为max{k xor&nbsp;ap | ap ≠ k , i ≤ p ≤ j}。 <br>
现在你需要知道你怎么选取需要融合的宝石，才能使生成的宝石能量密度最大。 <br>
</span></p></div><h2>Input</h2><div class="content"><p><font size="4">第一行，一个整数&nbsp;n，表示宝石个数。&nbsp;<br>
第二行，&nbsp;n个整数，分别表示a1至an，表示每颗宝石的能量密度，保证对于i ≠ j有&nbsp;ai ≠ aj。
</font></p></div><h2>Output</h2><div class="content"><p><font size="4">输出一行一个整数，表示最大能生成的宝石能量密度。 <br>
</font></p></div><h2>Sample Input</h2>
			<div class="content"><span class="sampledata">5 <br>
9  2 1 4 7
</span></div><h2>Sample Output</h2>
			<div class="content"><span class="sampledata">14 
 </span></div><h2>HINT</h2>
【样例解释】 
选择区间[1,5]，最大值为 7 xor 9。 
对于 100%的数据有 1 ≤ n ≤ 50000, 0 ≤ ai ≤ 10^9</p>
<hr>
</body>