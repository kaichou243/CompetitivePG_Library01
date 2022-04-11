using ll=long long;
struct Edge{
  ll to;
  ll cost;
};
using Graph=vector<vector<Edge>>;
template<typename T>
vector<T> shortest_path_faster_algorithm(const Graph &G, int s) {
  vector<T> dist(G.size(), INF);
  vector<int> pending(G.size(), 0), times(G.size(), 0);
  queue<int> que;
  que.emplace(s);
  pending[s]=true;
  ++times[s];
  dist[s]=0;
  while(!que.empty()) {
    int p = que.front();
    que.pop();
    pending[p] = false;
    for(auto &e : G[p]) {
      T next_cost = dist[p] + e.cost;
      if(next_cost >= dist[e.to]) continue;
      dist[e.to] = next_cost;
      if(!pending[e.to]) {
        if(++times[e.to] >= (int)G.size()) return vector<T>();
        pending[e.to] = true;
        que.emplace(e.to);
      }
    }
  }
  return dist;
}
