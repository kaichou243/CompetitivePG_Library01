//ABC010-D
#include<bits/stdc++.h>
using namespace std;
using FLOW=int;
const FLOW INF = 1000000000;

struct Edge {
    int rev, from, to;
    FLOW cap, icap;
    Edge(int r, int f, int t, FLOW c) : rev(r), from(f), to(t), cap(c), icap(c) {}
    friend ostream& operator << (ostream& s, const Edge& E) {
        if (E.cap > 0) return s << E.from << "->" << E.to << '(' << E.cap << ')';
        else return s;
    }
};

struct Graph {
    int V;
    vector<vector<Edge>> list;

    Graph(int n = 0) : V(n), list(n) { for (int i = 0; i < n; ++i) list[i].clear(); }
    void init(int n = 0) { V = n; list.resize(n); for (int i = 0; i < n; ++i) list[i].clear(); }
    void resize(int n = 0) { V = n; }
    void reset() { for (int i = 0; i < V; ++i) for (int j = 0; j < list[i].size(); ++j) list[i][j].cap = list[i][j].icap; }
    inline vector<Edge>& operator [] (int i) { return list[i]; }

    Edge &redge(Edge e) {
        if (e.from != e.to) return list[e.to][e.rev];
        else return list[e.to][e.rev + 1];
    }

    void addedge(int from, int to, FLOW cap) {
        list[from].push_back(Edge((int)list[to].size(), from, to, cap));
        list[to].push_back(Edge((int)list[from].size() - 1, to, from, 0));
    }
};

void dibfs(Graph &G, int s,vector<int> &level) {
    for (int i = 0; i < G.V; ++i) level[i] = -1;
    level[s] = 0;
    queue<int> que;
    que.push(s);
    while (!que.empty()) {
        int v = que.front();
        que.pop();
        for (int i = 0; i < G[v].size(); ++i) {
            Edge &e = G[v][i];
            if (level[e.to] < 0 && e.cap > 0) {
                level[e.to] = level[v] + 1;
                que.push(e.to);
            }
        }
    }
}

FLOW didfs(Graph &G, int v, int t, FLOW f,vector<int> &level,vector<int> &iter) {
    if (v == t) return f;
    for (int &i = iter[v]; i < G[v].size(); ++i) {
        Edge &e = G[v][i], &re = G.redge(e);
        if (level[v] < level[e.to] && e.cap > 0) {
            FLOW d = didfs(G, e.to, t, min(f, e.cap),level,iter);
            if (d > 0) {
                e.cap -= d;
                re.cap += d;
                return d;
            }
        }
    }
    return 0;
}
FLOW Dinic(Graph &G, int s, int t) {
    FLOW res = 0;
    vector<int> level(G.V),iter(G.V);
    while (true) {
        dibfs(G, s,level);
        if (level[t] < 0) return res;
        iter.assign(G.V,0);
        FLOW flow;
        while ((flow = didfs(G, s, t, INF,level,iter)) > 0) {
            res += flow;
        }
    }
}


int main() {
  int N,girl,E;
  cin>>N>>girl>>E;
  Graph G(N+1);
  int t=N;
  for(int i=0;i<girl;i++){
    int p;
    cin>>p;
    G.addedge(p,t,1);
  }
  for(int i=0;i<E;i++){
    int a,b;
    cin>>a>>b;
    G.addedge(a,b,1);
    G.addedge(b,a,1);
  }
  cout<<Dinic(G,0,t)<<endl;
}
