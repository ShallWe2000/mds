---
title: 【bzoj 2555】SubString
date: 2017/01/01 22:00:00
tags: 
  - 后缀自动机
  - LCT
categories: 
  - 字符串题

---

> 一个LCT+SAM的大板子题

<!--more-->

# 题目
<center><h2>2555: SubString</h2>
Time Limit: $30$ **Sec**  Memory Limit: $512$ **MB**</center>

## Description
懒得写背景了，给你一个字符串$init$，要求你支持两个操作
(1):在当前字符串的后面插入一个字符串
(2):询问字符串$s$在当前字符串中出现了几次？(作为连续子串)
你必须在线支持这些操作。

## Input
第一行一个数$Q$表示操作个数
第二行一个字符串表示初始字符串`init`
接下来$Q$行，每行2个字符串`Type,Str`
`Type`是`ADD`的话表示在后面插入字符串。
`Type`是`QUERY`的话表示询问某字符串在当前字符串中出现了几次。
为了体现在线操作，你需要维护一个变量`mask`，初始值为$0$
读入串`Str`之后，使用这个过程将之解码成真正询问的串`TrueStr`。
询问的时候，对`TrueStr`询问后输出一行答案`Result`
然后`mask = mask xor Result`  
插入的时候，将TrueStr插到当前字符串后面即可。
### HINT:`ADD`和`QUERY`操作的字符串都需要解压
$100\%$ 的数据字符串最终长度 $<=600000$，询问次数$<=10000$,询问总长度$<= 3000000$
## 数据范围
$100\%$ 的数据字符串最终长度$<=600000$，询问次数$<=10000$,询问总长度$<= 3000000$
# 解题报告
算法想起来真是不麻烦，在线插入字符用**后缀自动机**解决很简便，出现次数在**后缀自动机**中子串出现次数用`right`集合大小表示，`right`可以在`parent`树上统计，实际上是**查询子树和**操作；
**parent**树需要**删边+加边**，考虑用**LCT**搞，但是**LCT**不能统计子树信息，无所谓，当点$sm[x]+=1$,令$x$到根路径上所有点$sm[y]+=1$，这样需要**LCT**是一棵有根**LCT**,**splay**中的信息不从左右儿子转移，而是由**标记**从上向下修改；
上午调了好久，因为lct的**板子**有个地方打错了，就是`rotate`中所有循环终止条件都应该是`root(x)`，昨天打残`splay`，复习板子又把`lct`搞坏了，真是省选滚粗节奏呀!
# 代码
不是很长，但可以用**后缀平衡树或二分+普通平衡树+sa**取代，想学习一下；
```c++
#include<iostream>
#include<algorithm>
#include<cstdio>
#include<cstring>
#include<cstdlib>
using namespace std;
const int N=1200002; 
int mask=0,Q; 
char s[3000001],type[20];
string chars; 
struct LCT {
	int son[N][2],f[N],mark[N],sm[N],sta[N]; 
	bool root(int x) {
		return (son[f[x]][0]!=x&&son[f[x]][1]!=x); 
	}
	void add(int x,int y) {
		if (!x) return;
		sm[x]+=y, mark[x]+=y; 
	}
	void down(int x) {
		if (mark[x]) { 
			if (son[x][1]) add(son[x][1],mark[x]); 
			if (son[x][0]) add(son[x][0],mark[x]); 
			mark[x]=0; 
		}
	}
	void rotate(int x) {
		int y=f[x],z=f[y],d=(son[y][1]==x); 
		f[son[x][d^1]]=y,son[y][d]=son[x][d^1]; 
		if (!root(y)) son[z][son[z][1]==y]=x; 
		f[x]=z,f[y]=x,son[x][d^1]=y; 
	}
	void splay(int x) {
		int top=1,i=x; sta[top]=x; 
		while (!root(i)) sta[++top]=i=f[i]; 
		while (top) down(sta[top]),top--; 
		for (i=f[x];!root(x);rotate(x),i=f[x]) {
			if (root(i)) continue; 
			if ((son[i][1]==x)^(son[f[i]][1]==i)) rotate(x); 
			else rotate(i);
		}
	}
	void access(int x) {
		 for (int t=0;x;t=x,x=f[x]) 
		 	splay(x),son[x][1]=t;
	}
	void link(int x,int y) {
		f[x]=y,access(y),splay(y),add(y,sm[x]); 
	} 
	void cut(int x) {
		access(x),splay(x),add(son[x][0],-sm[x]); 
		f[son[x][0]]=0, son[x][0]=0; 
	}
} lct; 
struct SAM {
	int son[N][100],val[N],par[N]; 
	int p,q,np,nq,root,last,cnt; 
	SAM() { root=last=cnt=1; }
	void getstring(int mask) {
		scanf("%s",s),chars=s; 
		for (int j=0;j<chars.length();++j) {
        	mask=(mask*131+j)%chars.length();
        	char t=chars[j];
        	chars[j]=chars[mask],chars[mask]=t;
    	}
	}
	void insert(int x) {
		p=last,np=++cnt,last=np,lct.sm[np]=1,val[np]=val[p]+1; 
		while (p&&!son[p][x]) son[p][x]=np,p=par[p]; 
		if (!p) par[np]=root, lct.link(np,root); 
		else {
			q=son[p][x]; 
			if (val[p]+1==val[q]) par[np]=q,lct.link(np,q); 
			else {
				 nq=++cnt, val[nq]=val[p]+1,lct.sm[nq]=0; 
				 memcpy(son[nq],son[q],sizeof(son[q])); 
				 par[nq]=par[q],lct.link(nq,par[q]); 
				 par[q]=par[np]=nq,lct.cut(q); 
				 lct.link(q,nq),lct.link(np,nq);
				 while (p&&son[p][x]==q) son[p][x]=nq,p=par[p];
			}
		}
	}
	void build() { 
		scanf("%s",s); int l=strlen(s); 
		for (int i=0;i<l;++i) insert(s[i]-'A'); 
	}
	void add() {
		getstring(mask); int l=chars.length();
		for (int i=0;i<l;++i) 
			insert(chars[i]-'A'); 
	}
	int query() {
		getstring(mask); int l=chars.length(),i,x;
		for (i=0,x=root;i<l;++i)
			if (son[x][chars[i]-'A']) x=son[x][chars[i]-'A']; 
			else return 0; 
		lct.splay(x);  
		return lct.sm[x];
	}
} sam;
int main() {
	freopen("bzoj2555.in","r",stdin); 
	freopen("bzoj2555.out","w",stdout); 
	scanf("%d",&Q); sam.build(),mask=0;
	for (int i=1,tmp;i<=Q;++i) {
		scanf("%s",type); 
		if (type[0]=='A') sam.add(); 
		else {	
			tmp=sam.query();
			printf("%d\n",tmp),mask^=tmp; 
		}
	}
	return 0; 
}
```