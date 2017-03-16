---
title: 【bzoj 1604】[Usaco2008 Open]Cow Neighborhoods 奶牛的邻居
date: 2017/01/01 22:00:00
tags: 
  - 曼哈顿距离最小生成树
  - 平衡树
  - STL
categories: 
  - 数据结构

---

> 曼哈顿距离最小生成树的模板+奇怪的做法

<!--more-->

# 题目
<center><h2>1604: [Usaco2008 Open]Cow Neighborhoods 奶牛的邻居</h2>
Time Limit: 5 Sec  Memory Limit: 64 MB</center>

## Description

了解奶牛们的人都知道，奶牛喜欢成群结队．观察约翰的$N(1≤N≤100000)$只奶牛，你会发现她们已经结成了几个“群”．每只奶牛在吃草的时候有一个独一无二的位置坐标$X_i,Y_i(l≤Xi,Yi≤[1..10^9]$；$X_i,Y_i∈Z$．当满足下列两个条件之一，两只奶牛$i$和$j$是属于同一个群的：
  1．两只奶牛的曼哈顿距离不超过$C(1≤C≤10^9)$，即$|X_i-X_j|+|Y_i-Y_j|<=c$.
  2．两只奶牛有共同的邻居．即，存在一只奶牛$k$，使$i$与$k$，$j$与$k$均同属一个群．
给出奶牛们的位置，请计算草原上有多少个牛群，以及最大的牛群里有多少奶牛
####Input
第$1$行输入$N$和$C$，之后$N$行每行输入一只奶牛的坐标．
####Output
仅一行，先输出牛群数，再输出最大牛群里的牛数，用空格隔开．
###解题报告
题目中隐含的<strong>关键字</strong>：<strong>曼哈顿距离</strong>，<strong>并查集可以维护的信息</strong>
加上点数很多，两两连边是$O(n^2)$级别，而最终答案最最理想的状态就是<strong>生成树</strong>的形态。
直接使用<strong>曼哈顿距离最小生成树</strong>，连好$4n$左右条边后，按照边权排序，选取可以使用的，用并查集维护信息。
<big>提交->AC,<strong>速度倒数！！</strong></big>
上网上搜了<strong>题解</strong>：标签是<strong>set+特殊的技巧</strong>（orz iwtwiioi)
做法是这样的，设曼哈顿距离为$M$ ,则有
$$M=|x_i-x_j|+|y_i-y_j|$$
拆绝对值得，
$$M=(x_i+y_i)-(x_j+y_j) (x_i>x_j,y_i>y_j)①$$
$$M=(x_i-y_i)-(x_j-y_j) (x_i>x_j,y_i<y_j)②$$
$$M=-(x_i-y_i)+(x_j-y_j) (x_i<x_j,y_i>y_j)③$$
$$M=-(x_i+y_i)+(x_j+y_j) (x_i<x_j,y_i<y_j)④$$
其中①④，②③分别可以看成
$$M_1=|(x_i+y_i)-(x_j+y_j)|$$
$$M_2=|(x_i-y_i)-(x_j-y_j)|$$
两点的曼哈顿距离实际上是$M_1$和$M_2$的<strong>较大值</strong>
然后$M<=c <=> M_1<=c且M_2<=c$
我们考虑对满足这个条件的点对连边；
具体可以将点坐标重设为$(x_i+y_i,x_i-y_i)$按照横坐标排序后，维护一个队列，队列中的点横坐标之差小于等于$c$,并将这些点的纵坐标加入`multiset`，每次加入的时候，如果与两侧的点差$<=c$,就连边。
<strong>细节</strong>：
1. 为什么只判断两侧的点： 设当前点在`multiset`中下标为$x$,如果$x$与$x+2$差在$c$之内，那么一定存在$x \to x+1$和$x+1 \to x+2$,连通性上没有影响；
2. `multiset`: 可重元素集合，结构体需要重载运算符，删除时，`erase(data)`会删除与`data`相等的所有元素，`erase(find(data))`只删除一个；
###代码
####曼哈顿最小生成树
```c++
#include<iostream>
#include<cstdio>
#include<cstring>
#include<algorithm>
#include<cstdlib>
using namespace std;
const int N=100005;
const int inf=2000000000;
struct E {
    int x,y,v;
    E(int x=0,int y=0,int v=0)
        :x(x),y(y),v(v) {}
} e[N<<2];
int bit[N],pos[N],id[N],vec[N],sz[N];
int x[N],y[N],tot,cnt,n,k,f[N];
inline bool cmpi(int a,int b) {
    return x[a]==x[b]?y[a]<y[b]:x[a]<x[b]; 
}
inline bool cmpe(E a,E b) {
    return a.v<b.v; 
}
inline int dis(int a,int b) {
    return abs(x[a]-x[b])+abs(y[a]-y[b]); 
}
inline void add(int x,int y,int dis) {
    e[++tot]=E(x,y,dis); 
}
inline void update(int x,int val,int pur) {
    for (;x;x-=x&-x)
        if (val<=bit[x]) 
            bit[x]=val,pos[x]=pur; 
}
inline int query(int x) {
    int tmp=-1,mn=inf;
    for (;x<=cnt;x+=x&-x) 
        if (bit[x]<mn)
            mn=bit[x],tmp=pos[x]; 
    return tmp;
}
inline int hash(int x) {
    return lower_bound(vec+1,vec+1+cnt,x)-vec;
}
inline int find(int x) {
    if (f[x]!=x) f[x]=find(f[x]); 
    return f[x];
}
void Mmst() {
    tot=0; 
    for (int dir=0;dir<4;dir++) {
        if (dir==1||dir==3)
            for (int i=1;i<=n;i++) 
                swap(x[i],y[i]); 
        else
            for (int i=1;i<=n;i++)
                x[i]=-x[i]; 
        for (int i=1;i<=n;i++) 
            id[i]=i; 
        sort(id+1,id+1+n,cmpi);
        for (int i=1;i<=n;i++) 
            vec[i]=y[i]-x[i]; 
        sort(vec+1,vec+1+n); 
        cnt=unique(vec+1,vec+1+n)-vec-1; 
        for (int i=1;i<=n;++i) 
            bit[i]=inf,pos[i]=-1;
        int u,v;
        for (int i=n;i;i--) {
            u=hash(y[id[i]]-x[id[i]]); 
            v=query(u);
            if (v!=-1)
                add(id[i],v,dis(id[i],v)); 
            update(u,x[id[i]]+y[id[i]],id[i]); 
        }
    }
}
int main() {
//	freopen("graph.in","r",stdin); 
	scanf("%d%d",&n,&k);
    for (int i=1;i<=n;i++)
        scanf("%d%d",&x[i],&y[i]); 
    Mmst(); 
    sort(e+1,e+1+tot,cmpe);
    for (int i=1;i<=n;++i) f[i]=i;  
    int a,b,x,y,tmp=0;
    for (int i=1;i<=tot;i++) {
    	if (e[i].v>k) break;
        a=e[i].x,b=e[i].y; 
        if ((x=find(a))!=(y=find(b))) { 
			f[x]=y; 
		}
 	}
 	int ans1=0,ans2=0; 
	for (int i=1;i<=n;i++) { 
		if (!sz[find(i)]) ans1++;
		sz[find(i)]++; 
	}
	for (int i=1;i<=n;i++) 
		if (f[i]==i) ans2=max(ans2,sz[i]);
	printf("%d %d",ans1,ans2);
    return 0; 
}
```
####标算
```c++
#include<iostream>
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<algorithm>
#include<set>
using namespace std;
const int N=100001; 
const long long inf=10000000000000LL;
int n,c,ans,mx;
int f[N],sz[N];
struct point {
	long long x,y;int id;
	point(long long x=0,long long y=0,int id=0) 
		:x(x),y(y),id(id){}
} p[N];
multiset <point> s;
set <point>::iterator it;
inline bool operator<(point a,point b) {
    return a.y<b.y;
}
inline bool cmpx(point a,point b) {
    return a.x<b.x;
}
int find(int x) {
    return x==f[x]?x:f[x]=find(f[x]);
}
inline void merge(int x,int y) {
    int p=find(x),q=find(y);
    if(p!=q) f[p]=q,ans--;
}
void solve() {
    s.insert(point(0,inf,0));
	s.insert(point(0,-inf,0));   
    int now=1; s.insert(p[1]);
    for(int i=2;i<=n;i++) {
        while(p[i].x-p[now].x>c) {
            s.erase(s.find(p[now]));
            now++;
        }
        it=s.lower_bound(p[i]);
        point r=*it,l=*--it;
        if(p[i].y-l.y<=c)
            merge(p[i].id,l.id);
        if(r.y-p[i].y<=c)
            merge(p[i].id,r.id);
        s.insert(p[i]);
    }
}
inline void in(long long &x) {
	char ch=getchar(); int f=1;
	for (;ch<'0'||ch>'9';ch=getchar())
		if (ch=='-') f=-1;
	for (x=0;ch>='0'&&ch<='9';ch=getchar())
		x=x*10+ch-48; 
	x*=f;
}
int main() {
	scanf("%d%d",&n,&c);ans=n;
    for(int i=1;i<=n;i++) f[i]=i;
    long long x,y;
    for(int i=1;i<=n;i++) {
    	in(x),in(y); 
    	p[i]=point(x+y,x-y,i);
    }
    sort(p+1,p+n+1,cmpx);
    solve();
    for(int i=1;i<=n;i++) sz[find(i)]++;
    for(int i=1;i<=n;i++) mx=max(mx,sz[i]);
    printf("%d %d\n",ans,mx);
    return 0;
}
```
