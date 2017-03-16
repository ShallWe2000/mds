---
title: 【bzoj 1453】[Wc]Dface双面棋盘
date: 2017-02-14 14:20:27
tags:
  - LCT
  - 线段树
  - 并查集
categories:
  - 数据结构
---

> 康复数据结构---线段树+lCT 

<!--more-->

# 题目大意
黑白棋盘, 操作使 $(i, j)$ 变色, 然后查询两种颜色的联通块个数.
# 解题报告
Part One

经典动态图问题, 每次操作可以看做将 $(i, j)$ 与四周颜色相同的点间的边删除, 并与颜色不相同的点连边;

动态图的做法是使用LCT维护一棵删除时间最大生成树, 因为在一个联通块中, 删除时间较早的边, 在删除时间靠后的边加入后变得没有任何意义;

这个方法的前提是离线, 这个题就是一个离线问题;

先将相邻两点之间的边编号, 然后预处理每条边删除的时间点(我使用了队列), 然后每部操作, 分别进行删边和加边, 维护最大时间生成树;

太久没有打LCT, 使得我的`rotate`按照`splay`的套路打的, 结果使得虚父亲向儿子连了边... 
<img src="1.jpg" width="100" height="100">

Part Two

我自己想的naive and slow 的做法...  

把每一列看做一个节点, 建线段树, 这个节点要维护一个两侧位置并查集信息和两个颜色的联通个数,  因为容易发现, 相邻节点合并的时候, 联通块合并的情况只和两侧的联通性有关

合并的时候, 新的联通块个数等于两个节点的联通块个数和-合并的联通块个数... $O(n)$ 进行并查集的合并后, 对两侧的并查集信息进行重标号... 总时间复杂度是 $O(n^2logn)$ ;

# 代码
* LCT

```c++
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

typedef long long ll;

typedef pair<ll, ll> pll;
typedef vector<int> vi;
typedef pair<int, int> pii;

const int N = 210;
const int dx[4] = {1, -1, 0, 0};
const int dy[4] = {0, 0, 1, -1};

int a[N][N], n, m, ans[2], poi_edge[N * N][2], poi_id[N][N], cnt;
vector< pii > opera;
queue<int> on_edge[N * N * 2];
inline void in(int & x) {
    char ch = getchar(); int f = 1;
    for (; ch < '0' || ch > '9'; ch = getchar())
        if (ch == '-') f = -1;
    for (x = 0; ch >= '0' && ch <= '9'; ch = getchar())
        x = x * 10 + ch - 48;
    x *= f;
}
namespace LCT {
    int mn[N * N * 3];
    bool rev[N * N * 3];
    int f[N * N * 3], son[N * N * 3][2], stack[N * N * 2];
    void build(int x) {
        mn[x] = x;
        son[x][0] = son[x][1] = f[x] = rev[x] = 0;
    }
    inline bool compare(int x, int y) {
        x -= n * n, y -= n * n;
        if (x <= 0) return false;
        if (y <= 0) return true;
        if (on_edge[x].empty()) return false;
        if (on_edge[y].empty()) return true;
        return on_edge[x].front() < on_edge[y].front();
    }
    inline bool is_root(int x) {
        return (f[x] == 0) || (son[f[x]][0] != x && son[f[x]][1] != x);
    }
    inline void _rev(int x) {
        if (!x) return ;
        swap(son[x][0], son[x][1]), rev[x] ^= 1;
    }
    inline void push_down(int x) {
        if (rev[x]) {
            if (son[x][0]) _rev( son[x][0] );
            if (son[x][1]) _rev( son[x][1] );
            rev[x] = 0;
        }
    }
    inline void push_up(int x) {
        mn[x] = x;
        if (son[x][0])
            if (compare(mn[son[x][0]], mn[x]))
                mn[x] = mn[son[x][0]];
        if (son[x][1])
            if (compare(mn[son[x][1]], mn[x]))
                mn[x] = mn[son[x][1]];
    }
    inline void rotate(int x) {
        int y = f[x], z = f[y], d = (son[y][1] == x);
        son[y][d] = son[x][d^1];
        if (son[x][d ^ 1]) f[ son[x][d ^ 1] ] = y;
        if (!is_root(y)) son[z][ son[z][1] == y] = x;
        f[x] = z, son[x][d ^ 1] = y, f[y] = x;
        push_up(y);
    }

    inline void splay(int x) {
        int top=0, i = x; stack[++top] = x;
        while (!is_root(i)) stack[++top] = i = f[i];
        while (top) push_down(stack[top]), --top;
        for ( ; !is_root(x); rotate(x)) {
            int y = f[x];
            if (!is_root(y)) {
                if ((son[y][1] == x) ^ (son[f[y]][1] == y))
                    rotate(x);
                else rotate(y);
            }
        }
        push_up(x);
    }
    inline void access(int x) {
        for (int r = 0; x; r = x, x = f[x])
            splay(x), son[x][1] = r, push_up(x);
    }
    inline void mkroot(int x) {
        access(x), splay(x), _rev(x);
    }
    inline int root(int x) {
        access(x), splay(x);
        while (son[x][0])  x = son[x][0];
        return x;
    }
    inline void link(int x, int y) {
        mkroot(x), f[x] = y, access(x);
    }
    inline bool cut(int x, int y) {
        mkroot(x),  access(y), splay(y);
        if (son[y][0] != x || son[x][1] != 0) return 0;
        son[y][0] = f[x] = 0, push_up(y);
        return 1;
    }
    inline bool connect(int x, int  y) {
        int rtx = root(x), rty = root(y);
        return rtx == rty;
    }
}
inline bool legal(int x, int y) {
     return x >= 1 && x <= n && y >= 1 && y <= n;
}

inline int _edge(int x, int y, int d) {
    int nx = x + dx[d], ny = y + dy[d];
    return min(poi_id[x][y], poi_id[nx][ny])+ n * n + n * n * (d >> 1);
}
inline pii get_poi(int x) {
    int d = x > (n * n * 2);
    x = x - (1 + d) * n * n;
    return mp(x, (d? x + 1: x + n));
}

inline void _link(int x, int color) {
    using namespace LCT;
    pii tmp = get_poi(x);
    int a = tmp.fi,  b = tmp.se;
    if (!connect(a, b)) {
        link(a, x);
        link(b, x);
        -- ans[color];

    } else {
        mkroot(a), access(b), splay(b);
        int _x = mn[b];
        if (compare(_x, x)) {
            pii _tmp = get_poi(_x);
            int _a = _tmp.fi, _b = _tmp.se;
            cut(_a, _x), cut(_b, _x);
            link(a, x), link(b, x);
        }
    }
}
inline void _cut(int x, int color) {
    using namespace LCT;
    pii tmp = get_poi(x);
    int a = tmp.fi, b = tmp.se;
    if (cut(a, x) && cut(b, x))
        ++ ans[color];
    on_edge[x - n * n].pop();
}



int main() {
    in(n);
    FORU(i, 1, n) FORU(j, 1, n) {
        in(a[i][j]), poi_id[i][j] = ++cnt;
    }
    FORU(i, 1, n) FORU(j, 1, n) {
        if (i != n) {
            poi_edge[ poi_id[i][j] ][0] = poi_id[i][j] + n * n;
        }
        if (j != n) {
            poi_edge[ poi_id[i][j] ][1] = poi_id[i][j] + n * n * 2;
        }
    }// 点-> 边, 边-> 点
    FORU(i, 1, n * n * 3) LCT :: build(i);
    cin >> m; int x, y, nx, ny;
    FORU(i, 1, m) {
        in(x), in(y);
        opera.pb ( mp(x, y) );
        REP(d, 4) {
            nx = x + dx[d], ny = y + dy[d];
            if (legal(nx, ny))
                if (a[nx][ny] == a[x][y])
                    on_edge[_edge(x, y, d) - n * n].push(i);
        }
        a[x][y] ^= 1;
    }// 离线, 记录每条边的删除时间;


    FORU(i, 1, n) FORU(j, 1, n) REP(d, 4) {
        x = i + dx[d], y = j + dy[d];
        if (legal(x, y) && a[i][j] == a[x][y])
            on_edge[_edge(i, j, d) - n * n].push(m + 1);
        ++ d;
    }// 逮捕漏网之鱼
    FORD(i, m, 1) {
        x = opera[i - 1].fi, y = opera[i - 1].se;
        a[x][y] ^= 1;
    } // 还原悲剧现场

    FORU(i, 1, n) FORU(j, 1, n) {
        ans[ a[i][j] ] += 1;
        REP(d, 4) {
            x = i + dx[d], y = j + dy[d];
            if (legal(x, y) && a[x][y] == a[i][j])
                _link(_edge(i, j, d), a[i][j]);
            ++ d;
        }
    }// 初始连边
    int _x, _y;
    FORU(i, 1, m) {
        x = opera[i - 1].fi, y = opera[i - 1].se;
        REP(d, 4) {
            _x = x + dx[d], _y = y + dy[d];
            if (legal(_x, _y) && a[_x][_y] == a[x][y])
                _cut(_edge(x, y, d), a[x][y]);
        }
        -- ans[ a[x][y] ], a[x][y] ^= 1, ++ ans[ a[x][y] ];
        REP(d, 4) {
            _x = x + dx[d], _y = y + dy[d];
            if (legal(_x, _y) && a[_x][_y] == a[x][y])
                _link(_edge(x, y, d), a[x][y]);
        }
        printf("%d %d\n", ans[1], ans[0]);
    }
    return 0;
}
```
* 线段树

```c++
#include <bits/stdc++.h>

const int N = 210;

struct P {
    int c[2], f[N*2];
} T[N*4];
int n, m, i, j, a[N][N], f[N*4], t[N*4];
int F(int x){
    return f[x]==x ? x : f[x]=F(f[x]);
}
inline void cal(int x,int p) {
    int i, j = 1;
    T[x].c[a[p][1]] = 1, T[x].c[a[p][1]^1] = 0;
    T[x].f[1] = T[x].f[1+n] = 1;
    for(i = 2; i <= n; i++){
        if(a[p][i] != a[p][j])
            T[x].c[ a[p][j = i] ]++;
        T[x].f[i] = T[x].f[i+n] = j;
    }
}
inline void up(int x,int p){
    int l=x<<1, r=x<<1|1, i;
    for(i=0;i<2;i++)
        T[x].c[i] = T[l].c[i]+T[r].c[i];
    for(i = 1; i <= n*2; i++)
        f[i] = T[l].f[i];
    for(i = 1; i <= n*2; i++)
        f[i+n*2] = T[r].f[i] + n*2;
    for(i = 1; i <= n; i++)
        if(a[p][i] == a[p+1][i] && F(i+n) != F(i+n*2))
            T[x].c[a[p][i]]--,f[f[i+n]] = f[i+n*2];
    for(i = 1; i <= n*4; i++){
        f[i] = F(i);
        if (i <= n) t[f[i]] = i;
        if (i > n*3) t[f[i]] = i-n*2;
    }
    for(i = 1; i <= n; i++) T[x].f[i] = t[f[i]];
    for(i = 1; i <= n; i++) T[x].f[i+n] = t[f[i+n*3]];
}
void build(int x, int a, int b) {
    if(a==b) {cal(x,a);return;}
    int mid = (a+b) >> 1;
    build(x << 1, a, mid);
    build(x << 1 | 1, mid+1, b);
    up(x , mid);
}
void change(int x, int a, int b, int c) {
    if(a==b) {cal(x,a);return;}
    int mid = (a+b)>>1;
    if(c <= mid) change(x << 1, a, mid, c);
    else change(x << 1 | 1, mid + 1 , b , c);
    up(x, mid);
}
int main() {
    scanf("%d",&n);
    for(i = 1; i <= n; i++)
        for(j = 1; j <= n; j++)
            scanf("%d", &a[i][j]);
    build(1, 1, n);
    scanf("%d", &m);
    while (m--) {
        scanf("%d%d", &i, &j);
        a[i][j] ^= 1,change(1 , 1 , n , i);
        printf("%d %d\n", T[1].c[1], T[1].c[0]);
    }
    return 0;
}
```
