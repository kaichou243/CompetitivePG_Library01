struct MST{
  struct MSTEdge{
    ll u,v;
    ll cost;
    bool used;
    bool operator<(const MSTEdge& o) const {
      return cost < o.cost;
    }
  };
  int n;
  vector<MSTEdge> edges;
  MST(int sz) : n(sz), edges(sz) {}
  void add_edge(int u,int v,ll c){
    edges.push_back({u,v,c,false});
  }
  ll kruskal(){
    UnionFind uf(n);
    sort(all(edges));
    ll min_cost=0;
    for(int i=0;i<sz(edges);i++){
      auto& [u,v,c,used]=edges[i];
      if(!uf.same(u,v)){
        uf.unite(u,v);
        used=true;
        min_cost+=c;
      }
    }
    return min_cost;
  }
  Graph Tree(){
    kruskal();
    Graph G(n);
    for(int i=0;i<sz(edges);i++){
      if(edges[i].used){
        G[edges[i].u].push_back({edges[i].v,edges[i].cost});
        G[edges[i].v].push_back({edges[i].u,edges[i].cost});
      }
    }
    return G;
  }
};
