---
title: 【bzoj 1178】[Apio2009]CONVENTION会议中心
date: 2017/01/01 22:00:00
tags: 
  - 倍增
  - 线段树
  - 平衡树
  - STL
categories: 
  - 数据结构

---

> 经典贪心算法变得这么难？

<!--more-->

# 题目
<center><h2>1178: [Apio2009]CONVENTION会议中心</h2>

Time Limit: 15 Sec  Memory Limit: 162 MB</center>

## Description

Siruseri政府建造了一座新的会议中心。许多公司对租借会议中心的会堂很感兴趣，他们希望能够在里面举行会议。 对于一个客户而言，仅当在开会时能够独自占用整个会堂，他才会租借会堂。会议中心的销售主管认为：最好的策略应该是将会堂租借给<strong>尽可能多</strong>的客户。显然，有可能存在不止一种满足要求的策略。
<!--more--> 例如下面的例子。总共有$4$个公司。他们对租借会堂发出了请求，并提出了他们所需占用会堂的起止日期（如下表所示）。

| 公司   | 开始日期 | 结束日期   |
|-------|:------:|-----------:|
|   1   |   4    |     9     |
|   2   |   9    |     11    |
|   3   |   13   |     19    |
|   4   |   10   |     17    |

上例中，最多将会堂租借给两家公司。租借策略分别是租给公司$1$和公司$3$，或是公司$2$和公司$3$，也可以是公司$1$和公司$4$。注意会议中心一天最多租借给一个公司，所以公司$1$和公司$2$不能同时租借会议中心，因为他们在第九天重合了。 销售主管为了公平起见，决定按照如下的程序来确定选择何种租借策略：首先，将租借给客户数量最多的策略作为候选，将所有的公司按照他们发出请求的顺序编号。对于候选策略，将策略中的每家公司的编号按升序排列。最后，选出其中<strong>字典序最小</strong>的候选策略作为最终的策略。 例中，会堂最终将被租借给公司$1$和公司$3$：$3$个候选策略是$ \{(1,3),(2,3),(1,4)\}$。而在字典序中$(1,3) < (1,4) < (2,3)$。 你的任务是帮助销售主管确定应该将会堂租借给哪些公司。

## Input

输入的第一行有一个整数$N$，表示发出租借会堂申请的公司的个数。第$2$到第$N+1$行每行有$2$个整数。第$i+1$行的整数表示第$i$家公司申请租借的起始和终止日期。对于每个公司的申请，起始日期为不小于$1$的整数，终止日期为不大于$10^9$的整数。$N≤200000$

## Output

输出的第一行应有一个整数$M$，表示最多可以租借给多少家公司。第二行应列出$M$个数，表示最终将会堂租借给哪些公司。

## HINT
修复后数据:JudgeOnline/upload/201605/dd.rar
# 解题报告
前段时间做的题，感谢当时<strong>Rivendell学长</strong>给予的指导
有一个<strong>众所周知</strong>的<strong>贪心策略</strong>用来求最多选择多少个线段使他们两两不相交：去包含后，按照<strong>右端点</strong>排序，每次选择可行的右端点最靠左的线段;

 -  <big>简略</big>证明：①去包含，如果选择包含其他线段的线段，那么将这条线段替换为它所包含的线段，一定不会变劣；②选择右端点最靠左的方案，在<strong>方案数相同情况下，对剩余线段影响最小</strong>，与剩余线段交换不会变优；

现在要求<strong>字典序最小</strong>，考虑字典序概念的特殊性：从前向后出现第一个不等号即可确定。
所以按照字典序从小到大的顺序，如果①能够选择当前线段，并且②不会影响最终的答案，那么一定选择当前线段。

 -  能够选择：可以使用<strong>线段树</strong>将已经选择的线段染色，查询区间是否存在染色；也可以使用<strong>平衡树</strong>（`set`即可）,将选择的线段的左右端点插入，每次查询区间内是否有端点，区间是否在一对左右端点之间。
 -  不影响最终答案：摘出当前线段$[s,e]$所在的最大未选择区间$[ll,rr]$,当前线段将该区间划分为$[ll,s-1],[s,e],[e+1,rr]$ 设不强制放任何线段情况下，区间$[l,r]$中，最多选择的线段数为$Q_{l,r}$ ,那么当前线段不影响最终答案的条件就是 $Q_{ll,s-1}+1+Q_{e+1,rr}=Q_{ll,rr}$这很显然。

问题变成如何处理$Q_{l,r}$直接预处理不可能，但很好的一点是从$i$开始选择$j$个线段最近的结束点这个信息是可合并信息，可以<strong>倍增!</strong>
令$f[i][j]$表示从$i$点开始选择$2^j$个线段最近的结束位置，那么$f[i][j]=f[f[i][j-1]+1][j-1]$,有了这个$st$表，就可以$log^n$查询$Q_{l,r}$了。
总时间复杂度$O(log^n)$

# 代码
我以后会好好起变量名的...
```c++
#include <bits/stdc++.h>
#define M 200001
#define inf 0x3f3f3f3f
using namespace std;
int f[M<<1][19],l[M<<1],L[M<<1],N,cnt,m,n;
struct data {
    int s,e,id;
}a[M],b[M],c[M];
struct rec {
    int x,k;
}le,ri;
set<rec> s;
set<rec>::iterator itl,itr;
inline bool cmp1(data a,data b) {
    if (a.s==b.s) return a.e>b.e; 
    return a.s<b.s;
}
inline bool cmp2(data a,data b) {
    return a.id<b.id;
}
inline bool operator <(rec a,rec b) {
    return a.x<b.x;
}
inline void ST() {
    sort(l+1,l+1+m);
    for (int i=1;i<=m;i++)
        if (i==1||l[i]!=l[i-1])
            L[++N]=l[i];
    for (int i=1;i<=n;i++) {
        a[i].s=lower_bound(L+1,L+1+N,a[i].s)-L;
        a[i].e=lower_bound(L+1,L+1+N,a[i].e)-L;
    }
    sort(a+1,a+1+n,cmp1);
    int la=inf;
    for (int i=n;i;i--)
        if (a[i].e<la) b[++cnt]=a[i],la=a[i].e;
    for (int i=1;i<=cnt;i++)
        c[i]=b[cnt-i+1];
    memcpy(b,c,sizeof(c));
    memset(f,0x3f,sizeof(f));
    for (int i=N,j=cnt;i;i--) {
        f[i][0]=f[i+1][0];
        if (b[j].s==i) f[i][0]=min(f[i][0],b[j].e);
        for (int k=1;k<=17;k++)
            if (f[i][k-1]!=inf)
                f[i][k]=f[f[i][k-1]+1][k-1];
        while (b[j].s==i)
            j--;
    }
}
inline int Calc(int l,int r) {
    int ans=0;
    for (int i=17;i>=0;i--)
        if (f[l][i]<=r)
            l=f[l][i]+1,ans+=1<<i;
    return ans;
}
inline int read() {
	char ch=getchar();
	int x=0,f=1; 
	for (;ch<'0'||ch>'9';ch=getchar())
		if (ch=='-') f=-1; 
	for (;ch>='0'&&ch<='9';ch=getchar())
		x=x*10+ch-48;
	return x*f;
}
int main() {
	n=read();
    for (int i=1;i<=n;i++) {
    	a[i].s=read(),a[i].e=read();
	    a[i].id=i;
        l[++m]=a[i].s,l[++m]=a[i].e;
    }
    ST();int nn;
    printf("%d\n",nn=Calc(1,N));
    sort(a+1,a+1+n,cmp2);
    le.x=0,le.k=2,ri.x=N+1,ri.k=1;
    s.insert(le),s.insert(ri);
    for (int i=1;i<=n;i++) {
        le.x=a[i].s,le.k=1,ri.x=a[i].e,ri.k=2;
        itl=s.lower_bound(le);
        itr=s.upper_bound(ri);
        if (itl!=itr||itr->k==2) continue;
        itl--;
        int ll=itl->x+1,rr=itr->x-1;
        if (Calc(ll,a[i].s-1)+Calc(a[i].e+1,rr)+1!=Calc(ll,rr))
            continue;
        s.insert(le),s.insert(ri);
        nn--;
        if (nn)
        	printf("%d ",a[i].id);
        else 
        	printf("%d\n",a[i].id);
    }
    return 0;
}
```