#include<bits/stdc++.h>
using ll=long long;
ll INF=1e18;
struct Edge{
  ll to;
  ll cost;
};
using Graph=vector<vector<Edge>>;
template<typename T>
void dijkstra(const Graph &G,int s,vector<ll> &dist,vector<T> &cnt){
  int N = G.size();
  dist.assign(N, INF);
  cnt.assign(N,0);
  priority_queue<P, vector<P>, greater<P>> pq;
  dist[s] = 0;
  cnt[s] = 1;
  pq.emplace(dist[s], s);
  while (!pq.empty()){
    P p = pq.top();
    pq.pop();
    int v = p.second;
    if (dist[v] < p.first){
      continue;
    }
    for (auto e : G[v]){
      if (dist[e.to] > dist[v] + e.cost){
        dist[e.to] = dist[v] + e.cost;
        pq.emplace(dist[e.to], e.to);
        cnt[e.to]=cnt[v];
      }else if (dist[e.to] == dist[v] + e.cost){
        cnt[e.to]+=cnt[v];
      }
    }
  }
}
int main(){
}
