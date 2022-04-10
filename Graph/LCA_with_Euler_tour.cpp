struct SegP{
    int n;
    vector<P> dat;
    SegP(){}
    SegP(vector<P> v){
        int sz=v.size();
        n=1;
        while(n<sz) n*=2;
        dat.assign(2*n-1,{1e9,1e9});
        for(int i=0;i<sz;i++) dat[i+n-1]=v[i];
        for(int i=n-2;i>=0;i--) dat[i]=min(dat[2*i+1],dat[2*i+2]);
    }
    void resize(int N){
        n=1;
        while(n<N) n*=2;
        dat.assign(2*n-1,{1e9,1e9});
    }
    void vecset(vector<P> v){
        for(int i=0;i<v.size();i++){
            dat[i+n-1]=v[i];
        }
        for(int i=n-2;i>=0;i--) dat[i]=min(dat[2*i+1],dat[2*i+2]);
    }
    void update(int k,ll x){
        k+=(n-1);
        dat[k]={x,k};
        while(k>0){
            k=(k-1)/2;
            dat[k]=min(dat[2*k+1],dat[2*k+2]);
        }
    }
    P query(int l,int r,int k=0,int a=0,int b=-1){
        if(b<0) b=n;
        if(b<=l||r<=a) return make_pair(1e9,1e9);
        if(l<=a&&b<=r) return dat[k];
        P vl=query(l,r,2*k+1,a,(a+b)/2);
        P vr=query(l,r,2*k+2,(a+b)/2,b);
        return min(vl,vr);
    }
};
void Edfs(const Graph &G,int v,int p,vector<int> &ET){
  ET.push_back(v);
  for(auto nv : G[v]){
    if(nv.to==p) continue;
    Edfs(G,nv.to,v,ET);
    ET.push_back(v);
  }
}
void Ddfs(const Graph &G,int v,int p,vector<int> &depth){
  for(auto nv : G[v]){
    if(nv.to==p) continue;
    depth[nv.to]=depth[v]+1;
    Ddfs(G,nv.to,v,depth);
  }
}
//s=0
void Tbfs(const Graph &G,int s,vector<ll> &dist){
  int N=G.size();
  dist.assign(N,-1);
  dist[s]=0;
  deque<int> dq;
  dq.push_back(s);
  while(dq.size()){
    int v=dq[0];
    dq.pop_front();
    for(auto nv : G[v]){
      if(dist[nv.to]!=-1) continue;
      dist[nv.to]=dist[v]+nv.cost;
      dq.push_back(nv.to);
    }
  }
}
//first : p=-1
void dfs_sz(const vector<vector<int>> &G,int v,int p,vector<int> &sz){
  int ret=1;
  for(int nv : G[v]){
    if(nv==p) continue;
    dfs_sz(G,nv,v,sz);
    ret+=sz[nv];
  }
  sz[v]=ret;
}
struct LCA{
  vector<int> fst,ET,depth;
  vector<P> pv;
  SegP seg;
  vector<ll> dist;
  LCA(const Graph G){
    int n=G.size();
    fst.assign(n,-1);
    depth.resize(n);
    Edfs(G,0,-1,ET);
    Ddfs(G,0,-1,depth);
    for(int i=0;i<ET.size();i++){
      if(fst[ET[i]]==-1) fst[ET[i]]=i;
    }
    pv.resize(ET.size());
    for(int i=0;i<ET.size();i++){
      pv[i].first=depth[ET[i]];
      pv[i].second=ET[i];
    }
    seg.resize(pv.size());
    seg.vecset(pv);
    Tbfs(G,0,dist);
  }
  P lca(int u,int v){
    if(fst[u]>fst[v]) swap(u,v);
    P ret=seg.query(fst[u],fst[v]+1);
    return ret;
  }
  ll dis(int u,int v){
    P ca=lca(u,v);
    return dist[u]+dist[v]-2*dist[ca.second];
  }
};
