using ll=long long;
struct Edge{
  ll to;
  ll cost;
};
using Graph=vector<vector<Edge>>;
struct HLD{
  Graph G;
  vector<int> siz,a,ind,edge,dep,vertex,par;
  void dfs_sz(int v,int p){
    int ret=1;
    for(auto nv : G[v]){
      if(nv.to==p) continue;
      dfs_sz(nv.to,v);
      ret+=siz[nv.to];
    }
    siz[v]=ret;
  }
  void Ddfs(int v,int p){
    for(auto nv : G[v]){
      if(nv.to==p) continue;
      dep[nv.to]=dep[v]+1;
      Ddfs(nv.to,v);
    }
  }
  void hld(int v,int k,int p){
    ind[v]=vertex.size();
    a[v]=k;
    vertex.push_back(v);
    if(G[v].size()==1&&v!=0) return;
    int mxind=-1;
    int mx=0;
    for(int i=0;i<(int)G[v].size();i++){
      auto nv=G[v][i];
      if(nv.to==p) continue;
      if(siz[nv.to]>mx){
        mx=siz[nv.to];
        mxind=i;
      }
    }
    hld(G[v][mxind].to,k,v);
    for(int i=0;i<(int)G[v].size();i++){
      auto nv=G[v][i];
      if(nv.to==p) continue;
      if(i!=mxind) hld(nv.to,nv.to,v);
    }
  }
  void getpar(int v,int p){
    par[v]=p;
    for(auto nv : G[v]){
      if(nv.to==p) continue;
      getpar(nv.to,v);
      edge[ind[nv.to]]=nv.cost;
    }
  }
  HLD(const Graph &g) : siz((int)g.size()), a((int)g.size()), ind((int)g.size()), edge((int)g.size()), dep((int)g.size()), par((int)g.size()){
    G=g;
    dfs_sz(0,-1);
    hld(0,0,-1);
    getpar(0,-1);
    Ddfs(0,-1);
  }
  //[l,r]
  pair<int,int> getsubtree(int v){
    return make_pair(ind[v],ind[v]+siz[v]-1);
  }
  //[l,r]...
  vector<pair<int,int>> getpath(int u,int v){
    vector<pair<int,int>> ret;
    while(a[u]!=a[v]){
      if(dep[a[u]]<=dep[a[v]]){
        ret.push_back({ind[a[v]],ind[v]});
        v=par[a[v]];
      }else{
        ret.push_back({ind[a[u]],ind[u]});
        u=par[a[u]];
      }
    }
    ret.push_back({min(ind[u],ind[v]),max(ind[u],ind[v])});
    return ret;
  }
  int lca(int u,int v){
    vector<pair<int,int>> path=getpath(u,v);
    return vertex[path[path.size()-1].first];
  }
};
