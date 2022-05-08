#include<bits/stdc++.h>
#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
#define FOR(i,n) for(int i = 0; i < (n); i++)
#define sz(c) ((int)(c).size())
#define ten(x) ((int)1e##x)
#define all(v) (v).begin(), (v).end()
using namespace std;
using ll=long long;
using P = pair<ll,ll>;
const long double PI=acos(-1);
const ll INF=1e18;
const int inf=1e9;
struct Edge {
    ll to;
    ll cost;
};
using Graph=vector<vector<Edge>>;
template <typename T>
bool chmax(T &a,const T& b){
  if (a<b){
    a=b;
    return true;
  }
  return false;
}
template <typename T>
bool chmin(T &a,const T& b){
  if (a>b){
    a=b;
    return true;
  }
  return false;
}
template<int MOD> struct Fp{
  ll val;
  constexpr Fp(long long v = 0) noexcept : val(v % MOD) {
    if (val < 0) val += MOD;
  }
  static constexpr int getmod() { return MOD; }
  constexpr Fp operator - () const noexcept {
    return val ? MOD - val : 0;
  }
  constexpr Fp operator + (const Fp& r) const noexcept { return Fp(*this) += r; }
  constexpr Fp operator - (const Fp& r) const noexcept { return Fp(*this) -= r; }
  constexpr Fp operator * (const Fp& r) const noexcept { return Fp(*this) *= r; }
  constexpr Fp operator / (const Fp& r) const noexcept { return Fp(*this) /= r; }
  constexpr Fp& operator += (const Fp& r) noexcept {
    val += r.val;
    if (val >= MOD) val -= MOD;
    return *this;
  }
  constexpr Fp& operator -= (const Fp& r) noexcept {
    val -= r.val;
    if (val < 0) val += MOD;
    return *this;
  }
  constexpr Fp& operator *= (const Fp& r) noexcept {
    val = val * r.val % MOD;
    return *this;
  }
  constexpr Fp& operator /= (const Fp& r) noexcept {
    ll a = r.val, b = MOD, u = 1, v = 0;
    while (b) {
      ll t = a / b;
      a -= t * b, swap(a, b);
      u -= t * v, swap(u, v);
    }
    val = val * u % MOD;
    if (val < 0) val += MOD;
    return *this;
  }
  constexpr bool operator == (const Fp& r) const noexcept {
    return this->val == r.val;
  }
  constexpr bool operator != (const Fp& r) const noexcept {
    return this->val != r.val;
  }
  constexpr bool operator < (const Fp& r) const noexcept {
    return this->val < r.val;
  }
  friend constexpr istream& operator >> (istream& is, Fp<MOD>& x) noexcept {
    is >> x.val;
    x.val %= MOD;
    if (x.val < 0) x.val += MOD;
    return is;
  }
  friend constexpr ostream& operator << (ostream& os, const Fp<MOD>& x) noexcept {
    return os << x.val;
  }
  friend constexpr Fp<MOD> modpow(const Fp<MOD>& a, long long n) noexcept {
    Fp<MOD> res=1,r=a;
    while(n){
      if(n&1) res*=r;
      r*=r;
      n>>=1;
    }
    return res;
  }
  friend constexpr Fp<MOD> modinv(const Fp<MOD>& r) noexcept {
        long long a = r.val, b = MOD, u = 1, v = 0;
        while (b) {
            long long t = a / b;
            a -= t * b, swap(a, b);
            u -= t * v, swap(u, v);
        }
        return Fp<MOD>(u);
  }
  explicit operator bool()const{
		return val;
  }
};
struct SCCgraph{
  // input
  vector<vector<int>> G, rG;
  // result
  vector<vector<int>> scc;
  vector<int> cmp;
  vector<vector<int>> dag;
  // constructor
  SCCgraph(int N) : G(N), rG(N) {}
  // add edge
  void add_edge(int u, int v) {
    G[u].push_back(v);
    rG[v].push_back(u);
  }
  // decomp
  vector<bool> seen;
  vector<int> vs, rvs;
  void dfs(int v) {
    seen[v] = true;
    for (auto e : G[v]) if (!seen[e]) dfs(e);
    vs.push_back(v);
  }
  void rdfs(int v, int k) {
    seen[v] = true;
    cmp[v] = k;
    for (auto e : rG[v]) if (!seen[e]) rdfs(e, k);
    rvs.push_back(v);
  }
  // reconstruct
  set<pair<int,int>> newEdges;
  void reconstruct() {
    int N = (int)G.size();
    int dV = (int)scc.size();
    dag.resize(dV);
    newEdges.clear();
    for (int i = 0; i < N; ++i) {
      int u = cmp[i];
      for (auto e : G[i]) {
        int v = cmp[e];
        if (u == v) continue;
        if (!newEdges.count({u, v})) {
          dag[u].push_back(v);
          newEdges.insert({u, v});
        }
      }
    }
  }
  void solve() {
    // first dfs
    int N = (int)G.size();
    seen.assign(N, false);
    vs.clear();
    for (int v = 0; v < N; ++v) if (!seen[v]) dfs(v);
    // back dfs
    int k = 0;
    scc.clear();
    cmp.assign(N, -1);
    seen.assign(N, false);
    for (int i = N - 1; i >= 0; --i) {
      if (!seen[vs[i]]) {
        rvs.clear();
        rdfs(vs[i], k++);
        scc.push_back(rvs);
      }
    }
    // reconstruct
    reconstruct();
  }
};
template<class T,T (*op)(T,T),T (*e)()> struct SegmentTree{
  int n;
  vector<T> dat;
  SegmentTree(int N){
    n=1;
    while(n<N)n*=2;
    dat.assign(2*n,e());
  }
  void add(int k,T x){
    k+=n;
    dat[k]+=x;
    while(k){
      k>>=1;
      dat[k]=op(dat[k*2],dat[k*2+1]);
    }
  }
  void apply(int k,T x){
    k+=n;
    dat[k]=op(dat[k],x);
    while(k){
      k>>=1;
      dat[k]=op(dat[k*2],dat[k*2+1]);
    }
  }
  void set(int k,T x){
    k+=n;
    dat[k]=x;
    while(k){
      k>>=1;
      dat[k]=op(dat[k*2],dat[k*2+1]);
    }
  }
  T query(int l,int r){
    T prodl=e(),prodr=e();
    l+=n;
    r+=n;
    while(l<r){
      if(l&1) prodl=op(prodl,dat[l++]);
      if(r&1) prodr=op(dat[--r],prodr);
      l>>=1;
      r>>=1;
    }
    return op(prodl,prodr);
  }
};
struct FenwickTree{
    int n;
    vector<ll> dat;
    FenwickTree(int N){
        n=1;
        while(n<N) n*=2;
        dat.assign(n,0);
    }
    void add(int k,ll x){
    k+=1;
    while(k<=n){
      dat[k]+=x;
      k+=k&-k;
    }
  }
  //sum[0,k)
  ll sum(int k){
    ll ans=0;
    while(k){
      ans+=dat[k];
      k-=k&-k;
    }
    return ans;
  }
  //sum[l,r)
  ll sum(int l,int r){
    return sum(r)-sum(l);
  }
};
template<class S,S (*op)(S,S),S (*e)(),class F,S (*mapping)(F,S),F (*composition)(F,F),F (*id)()>
struct LazySegTree{
  private:
  int _n,size=1,idx=0;
  vector<S>seq;
  vector<F>lazy;
  void update(int k){seq[k]=op(seq[2*k],seq[2*k+1]);}
  void all_apply(int k,F f){
    seq[k]=mapping(f,seq[k]);
    if(k<size)lazy[k]=composition(f,lazy[k]);
  }
  void eval(int k){
    all_apply(2*k,lazy[k]);
    all_apply(2*k+1,lazy[k]);
    lazy[k]=id();
  }
  public:
  LazySegTree(int n):LazySegTree(vector<S>(n,e())){}
  LazySegTree(const vector<S>&v):_n(int(v.size())){
    while(size<_n)size<<=1,idx++;
    seq=vector<S>(2*size,e());
    lazy=vector<F>(2*size,id());
    for(int i=0;i<_n;i++)seq[size+i]=v[i];
    for(int i=size-1;i>=1;i--)update(i);
  }
  void set(int p,S x){
    p+=size;
    for(int i=idx;i>=1;i--)eval(p>>i);
    seq[p]=x;
    for(int i=1;i<=idx;i++)update(p>>i);
  }
  S operator[](int p){
    p+=size;
    for(int i=idx;i>=1;i--)eval(p>>i);
    return seq[p];
  }
  S query(int l,int r){
    if(l==r)return e();
    S sml=e(),smr=e();
    l+=size,r+=size;
    for(int i=idx;i>=1;i--){
      if(((l>>i)<<i)!=l)eval(l>>i);
      if(((r>>i)<<i)!=r)eval(r>>i);
    }
    while(l<r){
      if(l&1)sml=op(sml,seq[l++]);
      if(r&1)smr=op(seq[--r],smr);
      l>>=1,r>>=1;
    }
    return op(sml,smr);
  }
  S all_query()const{return seq[1];}
  void apply(int p,F f){
    p+=size;
    for(int i=idx;i>=1;i--)eval(p>>i);
    seq[p]=mapping(f,seq[p]);
    for(int i=1;i<=idx;i++)update(p>>i);
  }
  void apply(int l,int r,F f){
    if(l==r)return ;
    l+=size;
    r+=size;
    for(int i=idx;i>=1;i--){
      if(((l>>i)<<i)!=l)eval(l>>i);
      if(((r>>i)<<i)!=r)eval((r-1)>>i);
    }
    int l2=l,r2=r;
    while(l<r){
      if(l&1)all_apply(l++,f);
      if(r&1)all_apply(--r,f);
      l>>=1;
      r>>=1;
    }
    l=l2,r=r2;
    for(int i=1;i<=idx;i++){
      if(((l>>i)<<i)!=l)update(l>>i);
      if(((r>>i)<<i)!=r)update((r-1)>>i);
    }
  }
};
struct UnionFind{
  int n;
  vector<int> data;
  vector<int> edge_num;
  UnionFind(int N) : n(N) , data(N,-1) , edge_num(N,0){}
  int root(int x){ // データxが属する木の根を再帰で得る：root(x) = {xの木の根}
    return data[x]<0?x:data[x]=root(data[x]);
  }
  void unite(int x, int y) {
    x=root(x);
    y=root(y);
    if(x!=y){
      if (data[x]>data[y]) swap(x, y);
      data[x]+=data[y];
      data[y]=x;
    }
    if(x!=y){
      edge_num[x]+=edge_num[y];
      edge_num[y]=0;
    }
    edge_num[x]+=1;
  }
  bool same(int x, int y){ // 2つのデータx, yが属する木が同じならtrueを返す
    return root(x)==root(y);
  }
  int size(int x){
    return -data[root(x)];
  }
  int edge(int x){
    return edge_num[root(x)];
  }
  vector<int> get_roots(){
    vector<int> res;
    for(int i=0;i<n;i++){
      if(data[i]<0) res.push_back(i);
    }
    return res;
  }
};
template< typename T >
struct WeightedUnionFind {
  vector<int> data;
  vector<T> we;
  WeightedUnionFind(int sz) : data(sz, -1), we(sz) {}
  int root(int k) {
    if(data[k] < 0) return k;
    int par = root(data[k]);
    we[k] += we[data[k]];
    return data[k] = par;
  }
  T weight(int t) {
    root(t);
    return we[t];
  }
  bool unite(int x, int y, T w) {
    w += weight(x);
    w -= weight(y);
    x = root(x), y = root(y);
    if(x == y) return false;
    if(data[x] > data[y]) {
      swap(x, y);
      w *= -1;
    }
    data[x] += data[y];
    data[y] = x;
    we[y] = w;
    return true;
  }
  bool same(int x,int y){
    return root(x)==root(y);
  }
  int size(int x){
    return -data[root(x)];
  }
  T diff(int x, int y) {
    return weight(y) - weight(x);
  }
};
template< typename flow_t >
struct MFGraph{
  static_assert(is_integral< flow_t >::value, "template parameter flow_t must be integral type");
  const flow_t INF;
  struct edge {
    int to;
    flow_t cap;
    int rev;
    bool isrev;
    int idx;
  };

  vector< vector< edge > > graph;
  vector< int > min_cost, iter;
  flow_t max_cap;
  int SZ;

  explicit MFGraph(int V) : INF(numeric_limits< flow_t >::max()), graph(V), max_cap(0), SZ(V) {}

  void add_edge(int from, int to, flow_t cap, int idx = -1) {
    max_cap = max(max_cap, cap);
    graph[from].emplace_back((edge) {to, cap, (int) graph[to].size(), false, idx});
    graph[to].emplace_back((edge) {from, 0, (int) graph[from].size() - 1, true, idx});
  }

  bool build_augment_path(int s, int t, const flow_t &base) {
    min_cost.assign(graph.size(), -1);
    queue< int > que;
    min_cost[s] = 0;
    que.push(s);
    while(!que.empty() && min_cost[t] == -1) {
      int p = que.front();
      que.pop();
      for(auto &e : graph[p]) {
        if(e.cap >= base && min_cost[e.to] == -1) {
          min_cost[e.to] = min_cost[p] + 1;
          que.push(e.to);
        }
      }
    }
    return min_cost[t] != -1;
  }

  flow_t find_augment_path(int idx, const int t, flow_t base, flow_t flow) {
    if(idx == t) return flow;
    flow_t sum = 0;
    for(int &i = iter[idx]; i < (int)graph[idx].size(); i++) {
      edge &e = graph[idx][i];
      if(e.cap >= base && min_cost[idx] < min_cost[e.to]) {
        flow_t d = find_augment_path(e.to, t, base, min(flow - sum, e.cap));
        if(d > 0) {
          e.cap -= d;
          graph[e.to][e.rev].cap += d;
          sum += d;
          if(flow - sum < base) break;
        }
      }
    }
    return sum;
  }

  flow_t max_flow(int s, int t) {
    if(max_cap == flow_t(0)) return flow_t(0);
    flow_t flow = 0;
    for(int i = 63 - __builtin_clzll(max_cap); i >= 0; i--) {
      flow_t now = flow_t(1) << i;
      while(build_augment_path(s, t, now)) {
        iter.assign(graph.size(), 0);
        flow += find_augment_path(s, t, now, INF);
      }
    }
    return flow;
  }

  //after max_flow(s,t)
  vector<bool> min_cut(int s) {
        vector<bool> visited(SZ);
        queue<int> que;
        que.push(s);
        while (!que.empty()) {
            int p = que.front();
            que.pop();
            visited[p] = true;
            for (auto e : graph[p]) {
                if (e.cap && !visited[e.to]) {
                    visited[e.to] = true;
                    que.push(e.to);
                }
            }
        }
        return visited;
    }
  
  void output() {
    for(int i = 0; i < graph.size(); i++) {
      for(auto &e : graph[i]) {
        if(e.isrev) continue;
        auto &rev_e = graph[e.to][e.rev];
        cout << i << "->" << e.to << " (flow: " << rev_e.cap << "/" << e.cap + rev_e.cap << ")" << endl;
      }
    }
  }
};
template< typename flow_t, typename cost_t >
struct MCFGraph{
  struct edge {
    int to;
    flow_t cap;
    cost_t cost;
    int rev;
    bool isrev;
  };

  vector< vector< edge > > graph;
  vector< cost_t > potential, min_cost;
  vector< int > prevv, preve;
  const cost_t INF;

  MCFGraph(int V) : graph(V), INF(numeric_limits< cost_t >::max()) {}

  void add_edge(int from, int to, flow_t cap, cost_t cost) {
    graph[from].emplace_back((edge) {to, cap, cost, (int) graph[to].size(), false});
    graph[to].emplace_back((edge) {from, 0, -cost, (int) graph[from].size() - 1, true});
  }

  cost_t min_cost_flow(int s, int t, flow_t f) {
    int V = (int) graph.size();
    cost_t ret = 0;
    using Pi = pair< cost_t, int >;
    priority_queue< Pi, vector< Pi >, greater< Pi > > que;
    potential.assign(V, 0);
    preve.assign(V, -1);
    prevv.assign(V, -1);

    while(f > 0) {
      min_cost.assign(V, INF);
      que.emplace(0, s);
      min_cost[s] = 0;
      while(!que.empty()) {
        Pi p = que.top();
        que.pop();
        if(min_cost[p.second] < p.first) continue;
        for(int i = 0; i < (int)graph[p.second].size(); i++) {
          edge &e = graph[p.second][i];
          cost_t nextCost = min_cost[p.second] + e.cost + potential[p.second] - potential[e.to];
          if(e.cap > 0 && min_cost[e.to] > nextCost) {
            min_cost[e.to] = nextCost;
            prevv[e.to] = p.second, preve[e.to] = i;
            que.emplace(min_cost[e.to], e.to);
          }
        }
      }
      if(min_cost[t] == INF) return -1;
      for(int v = 0; v < V; v++) potential[v] += min_cost[v];
      flow_t addflow = f;
      for(int v = t; v != s; v = prevv[v]) {
        addflow = min(addflow, graph[prevv[v]][preve[v]].cap);
      }
      f -= addflow;
      ret += addflow * potential[t];
      for(int v = t; v != s; v = prevv[v]) {
        edge &e = graph[prevv[v]][preve[v]];
        e.cap -= addflow;
        graph[v][e.rev].cap += addflow;
      }
    }
    return ret;
  }

  void output() {
    for(int i = 0; i < graph.size(); i++) {
      for(auto &e : graph[i]) {
        if(e.isrev) continue;
        auto &rev_e = graph[e.to][e.rev];
        cout << i << "->" << e.to << " (flow: " << rev_e.cap << "/" << rev_e.cap + e.cap << ")" << endl;
      }
    }
  }
};
ll mod(ll a, ll mod) {
    return (a%mod+mod)%mod;
}
ll modpow(ll a,ll n,ll mod){
  ll res=1;
  a%=mod;
  while (n>0){
    if (n & 1) res*=a;
    a *= a;
    a%=mod;
    n >>= 1;
    res%=mod;
  }
  return res;
}
vector<P> prime_factorize(ll N) {
  vector<P> res;
  for (ll a = 2; a * a <= N; ++a) {
    if (N % a != 0) continue;
    ll ex = 0;
    while(N % a == 0){
      ++ex;
      N /= a;
    }
    res.push_back({a, ex});
  }
  if (N != 1) res.push_back({N, 1});
  return res;
}
ll modinv(ll a, ll mod) {
    ll b = mod, u = 1, v = 0;
    while (b) {
        ll t = a/b;
        a -= t * b, swap(a, b);
        u -= t * v, swap(u, v);
    }
    u %= mod;
    if (u < 0) u += mod;
    return u;
}
ll extGcd(ll a, ll b, ll &p, ll &q) {  
    if (b == 0) { p = 1; q = 0; return a; }  
    ll d = extGcd(b, a%b, q, p);  
    q -= a/b * p;  
    return d;  
}
P ChineseRem(const vector<ll> &b, const vector<ll> &m) {
  ll r = 0, M = 1;
  for (int i = 0; i < (int)b.size(); ++i) {
    ll p, q;
    ll d = extGcd(M, m[i], p, q);
    if ((b[i] - r) % d != 0) return make_pair(0, -1);
    ll tmp = (b[i] - r) / d * p % (m[i]/d);
    r += M * tmp;
    M *= m[i]/d;
  }
  return make_pair(mod(r, M), M);
}
template<int m> struct Comb{
  unordered_map<int,tuple<ll,vector<ll>,vector<ll>>> mp;
  int n_;
  ll p_, pm_;
  vector<ll> ord_, fact_;
  vector<P> pf;
  Comb(int n) : n_(n), ord_(n), fact_(n){ 
    pf=prime_factorize(m);
    COMinit();
  }
  void init(int n) {
    ord_.resize(n);
    fact_.resize(n);
  }
  void init(long long p, long long pm) {
    p_=p,pm_=pm;
    ord_[0]=ord_[1]=0;
    fact_[0]=fact_[1]=1;
    auto&[pms,ord,fac]=mp[p];
    pms=pm;
    ord.resize(n_);
    fac.resize(n_);
    ord[0]=ord[1]=0;
    fac[0]=fac[1]=1;
    for (int i=2;i<n_;i++) {
      long long add=0;
      long long ni=i;
      while (ni % p == 0) add++,ni/=p;
      ord_[i]=ord_[i-1]+add;
      fact_[i]=fact_[ni-1]*ni%pm;
      ord[i]=ord_[i];
      fac[i]=fact_[i];
    }
  }
  void init(long long p, long long pm, int n) {
    init(n);
    init(p, pm);
  }
  void COMinit(){
    for(auto p : pf){
      ll ps=p.first,pfs=pow(p.first,p.second);
      init(n_);
      init(ps,pfs);
    }
  }
  ll com(ll n, ll r,int p) {
    if (n<0 || r<0 || n<r) return 0;
    auto&[pms,ord,fac]=mp[p];
    ll e=ord[n]-ord[r]-ord[n-r];
    ll res=fac[n]*modinv(fac[r]*fac[n-r]%pms,pms)%pms;
    res=res*modpow(p,e,pms)%pms;
    return res;
  }
  ll operator()(int n, int k){
    if(n<0 || k<0 || n<k) return 0;
    int sz=pf.size();
    vector<long long> vb(sz), vm(sz);
    for (int i=0;i<sz;i++) {
        long long p = pf[i].first, e = pf[i].second;
        long long pm = pow(p,e);
        long long b = 1;
        b *= com(n, k ,p) % pm;
        b %= pm;
        vm[i]=pm;
        vb[i]=b;
    }
    auto res = ChineseRem(vb,vm);
    return res.first;
  }
};
void dfs(const vector<vector<int>> &G,int v,vector<bool> &seen){
  seen[v]=true;
  for(auto nv : G[v]){
    if(seen[nv]) continue;
    dfs(G,nv,seen);
  }
}
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
void bfs(const vector<vector<int>> &G,int s,vector<int> &dist){
  int N=G.size();
  dist.assign(N,-1);
  dist[s]=0;
  deque<int> dq;
  dq.push_back(s);
  while(dq.size()){
    int v=dq[0];
    dq.pop_front();
    for(auto nv : G[v]){
      if(dist[nv]!=-1) continue;
      dist[nv]=dist[v]+1;
      dq.push_back(nv);
    }
  }
}
vector<ll> enum_divisors(ll N){
  vector<ll> res;
  for (ll i = 1; i * i <= N; ++i){
    if (N % i == 0){
       res.push_back(i);
       // 重複しないならば i の相方である N/i も push
       if (N/i != i) res.push_back(N/i);
    }
  }
  // 小さい順に並び替える
  sort(res.begin(), res.end());
  return res;
}
ll Euler(ll n){
  vector<pair<ll,ll>> pf=prime_factorize(n);
  ll res=n;
  for(auto p : pf){
    res*=(p.first-1);
    res/=p.first;
  }
  return res;
}
//dist must be initialized with INF
void warshall_floyd(vector<vector<ll>> &dist,int n){
  for(int i=0;i<n;i++) dist[i][i]=0;
  for (int k = 0; k < n; k++){
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
      }
    }
  }
}
//[1,N]
vector<int> Eratosthenes(int N) {
  vector<int> isprime(N+1,1);
  isprime[1]=0;
  for (int p=2;p<=N;++p){
    if (!isprime[p]) continue;
    for (int q=p*2;q<=N;q+=p) {
      isprime[q]=0;
    }
  }
  return isprime;
}
//log X 素因数分解
struct PrimeFact{
  vector<int> era;
  PrimeFact(int maxA) : era(maxA+1,-1){
    for(int p=2;p<=maxA;++p){
      if(era[p]!=-1) continue;
      for(int q=p;q<=maxA;q+=p){
        era[q]=p;
      }
    }
  }
  vector<P> prime_fact(int x){
    vector<P> res;
    while(1<x){
      int p=era[x],ex=0;
      while(x%p==0){
        x/=p;
        ex++;
      }
      res.push_back({p,ex});
    }
    reverse(res.begin(),res.end());
    return res;
  }
};
namespace FFT {
  using Real = long double;
  struct Comp{
    Real re, im;
    Comp(){};
    Comp(Real re, Real im) : re(re), im(im) {}
    Real slen() const { return re*re+im*im; }
    Real real() { return re; }
    inline Comp operator+(const Comp &c) { return Comp(re+c.re, im+c.im); }
    inline Comp operator-(const Comp &c) { return Comp(re-c.re, im-c.im); }
    inline Comp operator*(const Comp &c) { return Comp(re*c.re-im*c.im, re*c.im+im*c.re); }
    inline Comp operator/(const Real &c) { return Comp(re/c, im/c); } 
    inline Comp conj() const { return Comp(re, -im); }
  };
  void fft(vector<Comp>& a, bool inverse){
    int n = a.size();
    static bool prepared = false;
    static vector<Comp> z(30);
    if(!prepared){
      prepared = true;
      for(int i=0; i<30; i++){
        Real ang = 2*PI/(1 << (i+1));
        z[i] = Comp(cos(ang), sin(ang));
      }
    }
    for (size_t i = 0, j = 1; j < n; ++j) {
      for (size_t k = n >> 1; k > (i ^= k); k >>= 1);
      if (i > j) swap(a[i], a[j]);
    }
    for (int i = 0, t = 1; t < n; ++i, t <<= 1) {
      Comp bw = z[i];
      if(inverse) bw = bw.conj();
      for (int i = 0; i < n; i += t * 2) {
        Comp w(1,0);
        for (int j = 0; j < t; ++j) {
          int l = i + j, r = i + j + t;
          Comp c = a[l], d = a[r] * w;
          a[l] = c + d;
          a[r] = c - d;
          w = w * bw;
        }
      }
    }
    if(inverse)
      for(int i=0; i<n; i++) a[i] = a[i]/(Real)n;
  }
  template<typename T>
  vector<T> mul(vector<T> &a, vector<T> &b, bool issquare){
    int deg = a.size() + b.size();
    int n = 1;
    while(n < deg) n <<= 1;
    vector<Comp> fa,fb;
    fa.resize(n); fb.resize(n);
    for(int i=0;i<n;i++){
      if(i < a.size()) fa[i] = Comp(a[i], 0);
      else fa[i] = Comp(0,0);
      if(i < b.size()) fb[i] = Comp(b[i], 0);
      else fb[i] = Comp(0,0);
    }
    for(int i=0;i<b.size();i++) fb[i] = Comp(b[i], 0);
    fft(fa, false); 
    if(issquare) fb = fa;
    else fft(fb, false);
    for(int i=0;i<n;i++) fa[i] = fa[i] * fb[i];
    fft(fa, true);
    vector<T> res(n);
    for(int i=0;i<n;i++) res[i] = round(fa[i].re);
    return res;
  }
};
namespace NTT {
    int calc_primitive_root(int mod) {
        if (mod == 2) return 1;
        if (mod == 167772161) return 3;
        if (mod == 469762049) return 3;
        if (mod == 754974721) return 11;
        if (mod == 998244353) return 3;
        int divs[20] = {};
        divs[0] = 2;
        int cnt = 1;
        long long x = (mod - 1) / 2;
        while (x % 2 == 0) x /= 2;
        for (long long i = 3; i * i <= x; i += 2) {
            if (x % i == 0) {
                divs[cnt++] = i;
                while (x % i == 0) x /= i;
            }
        }
        if (x > 1) divs[cnt++] = x;
        for (int g = 2;; g++) {
            bool ok = true;
            for (int i = 0; i < cnt; i++) {
                if (modpow(g, (mod - 1) / divs[i], mod) == 1) {
                    ok = false;
                    break;
                }
            }
            if (ok) return g;
        }
    }
 
    int get_fft_size(int N, int M) {
        int size_a = 1, size_b = 1;
        while (size_a < N) size_a <<= 1;
        while (size_b < M) size_b <<= 1;
        return max(size_a, size_b) << 1;
    }
    constexpr int bsf_constexpr(unsigned int n) {
      int x = 0;
      while (!(n & (1 << x))) x++;
      return x;
    }
    int bsf(unsigned int n) {
      #ifdef _MSC_VER
      unsigned long index;
      _BitScanForward(&index, n);
      return index;
      #else
      return __builtin_ctz(n);
      #endif
    }
    template <class mint>
    struct fft_info{
      static constexpr int rank2 = bsf_constexpr(mint::getmod() - 1);
      std::array<mint, rank2 + 1> root;   // root[i]^(2^i) == 1
      std::array<mint, rank2 + 1> iroot;  // root[i] * iroot[i] == 1
      std::array<mint, std::max(0, rank2 - 2 + 1)> rate2;
      std::array<mint, std::max(0, rank2 - 2 + 1)> irate2;
 
      std::array<mint, std::max(0, rank2 - 3 + 1)> rate3;
      std::array<mint, std::max(0, rank2 - 3 + 1)> irate3;
      int g;
      fft_info(){
        int MOD=mint::getmod();
        g=calc_primitive_root(MOD);
        root[rank2] = modpow(mint(g),(MOD - 1) >> rank2);
        iroot[rank2] = modinv(root[rank2]);
        for (int i = rank2 - 1; i >= 0; i--) {
            root[i] = root[i + 1] * root[i + 1];
            iroot[i] = iroot[i + 1] * iroot[i + 1];
        }
 
        {
            mint prod = 1, iprod = 1;
            for (int i = 0; i <= rank2 - 2; i++) {
                rate2[i] = root[i + 2] * prod;
                irate2[i] = iroot[i + 2] * iprod;
                prod *= iroot[i + 2];
                iprod *= root[i + 2];
            }
        }
        {
            mint prod = 1, iprod = 1;
            for (int i = 0; i <= rank2 - 3; i++) {
                rate3[i] = root[i + 3] * prod;
                irate3[i] = iroot[i + 3] * iprod;
                prod *= iroot[i + 3];
                iprod *= root[i + 3];
            }
        }
      }
    };
    int ceil_pow2(int n) {
      int x = 0;
      while ((1U << x) < (unsigned int)(n)) x++;
      return x;
    }
    // number-theoretic transform
    template <class mint>
    void trans(std::vector<mint>& a) {
      int n = int(a.size());
      int h = ceil_pow2(n);
      int MOD=a[0].getmod();
      static const fft_info<mint> info;
 
      int len = 0;  // a[i, i+(n>>len), i+2*(n>>len), ..] is transformed
      while (len < h) {
        if (h - len == 1) {
            int p = 1 << (h - len - 1);
            mint rot = 1;
            for (int s = 0; s < (1 << len); s++) {
                int offset = s << (h - len);
                for (int i = 0; i < p; i++) {
                    auto l = a[i + offset];
                    auto r = a[i + offset + p] * rot;
                    a[i + offset] = l + r;
                    a[i + offset + p] = l - r;
                }
                if (s + 1 != (1 << len))
                    rot *= info.rate2[bsf(~(unsigned int)(s))];
            }
            len++;
        } else {
            // 4-base
            int p = 1 << (h - len - 2);
            mint rot = 1, imag = info.root[2];
            for (int s = 0; s < (1 << len); s++) {
                mint rot2 = rot * rot;
                mint rot3 = rot2 * rot;
                int offset = s << (h - len);
                for (int i = 0; i < p; i++) {
                    auto mod2 = 1ULL * MOD * MOD;
                    auto a0 = 1ULL * a[i + offset].val;
                    auto a1 = 1ULL * a[i + offset + p].val * rot.val;
                    auto a2 = 1ULL * a[i + offset + 2 * p].val * rot2.val;
                    auto a3 = 1ULL * a[i + offset + 3 * p].val * rot3.val;
                    auto a1na3imag =
                        1ULL * mint(a1 + mod2 - a3).val * imag.val;
                    auto na2 = mod2 - a2;
                    a[i + offset] = a0 + a2 + a1 + a3;
                    a[i + offset + 1 * p] = a0 + a2 + (2 * mod2 - (a1 + a3));
                    a[i + offset + 2 * p] = a0 + na2 + a1na3imag;
                    a[i + offset + 3 * p] = a0 + na2 + (mod2 - a1na3imag);
                }
                if (s + 1 != (1 << len))
                    rot *= info.rate3[bsf(~(unsigned int)(s))];
            }
            len += 2;
        }
      }
    }
    template <class mint>
    void trans_inv(std::vector<mint>& a) {
      int n = int(a.size());
      int h = ceil_pow2(n);
 
      static const fft_info<mint> info;
      int MOD=a[0].getmod();
      int len = h;  // a[i, i+(n>>len), i+2*(n>>len), ..] is transformed
      while (len) {
        if (len == 1) {
            int p = 1 << (h - len);
            mint irot = 1;
            for (int s = 0; s < (1 << (len - 1)); s++) {
                int offset = s << (h - len + 1);
                for (int i = 0; i < p; i++) {
                    auto l = a[i + offset];
                    auto r = a[i + offset + p];
                    a[i + offset] = l + r;
                    a[i + offset + p] =
                        (unsigned long long)(MOD + l.val - r.val) *
                        irot.val;
                    ;
                }
                if (s + 1 != (1 << (len - 1)))
                    irot *= info.irate2[bsf(~(unsigned int)(s))];
            }
            len--;
        } else {
            // 4-base
            int p = 1 << (h - len);
            mint irot = 1, iimag = info.iroot[2];
            for (int s = 0; s < (1 << (len - 2)); s++) {
                mint irot2 = irot * irot;
                mint irot3 = irot2 * irot;
                int offset = s << (h - len + 2);
                for (int i = 0; i < p; i++) {
                    auto a0 = 1ULL * a[i + offset + 0 * p].val;
                    auto a1 = 1ULL * a[i + offset + 1 * p].val;
                    auto a2 = 1ULL * a[i + offset + 2 * p].val;
                    auto a3 = 1ULL * a[i + offset + 3 * p].val;
 
                    auto a2na3iimag =
                        1ULL *
                        mint((MOD + a2 - a3) * iimag.val).val;
 
                    a[i + offset] = a0 + a1 + a2 + a3;
                    a[i + offset + 1 * p] =
                        (a0 + (MOD - a1) + a2na3iimag) * irot.val;
                    a[i + offset + 2 * p] =
                        (a0 + a1 + (MOD - a2) + (MOD - a3)) *
                        irot2.val;
                    a[i + offset + 3 * p] =
                        (a0 + (MOD - a1) + (MOD - a2na3iimag)) *
                        irot3.val;
                }
                if (s + 1 != (1 << (len - 2)))
                    irot *= info.irate3[bsf(~(unsigned int)(s))];
            }
            len -= 2;
        }
      }
    }
    // for garner
    static constexpr int MOD0 = 754974721;
    static constexpr int MOD1 = 167772161;
    static constexpr int MOD2 = 469762049;
    using mint0 = Fp<MOD0>;
    using mint1 = Fp<MOD1>;
    using mint2 = Fp<MOD2>;
    static const mint1 imod0 = 95869806; // modinv(MOD0, MOD1);
    static const mint2 imod1 = 104391568; // modinv(MOD1, MOD2);
    static const mint2 imod01 = 187290749; // imod1 / MOD0;
 
    // small case (T = mint, long long)
    template<class T> vector<T> naive_mul 
    (const vector<T> &A, const vector<T> &B) {
        if (A.empty() || B.empty()) return {};
        int N = (int)A.size(), M = (int)B.size();
        vector<T> res(N + M - 1);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < M; ++j)
                res[i + j] += A[i] * B[j];
        return res;
    }
 
    // mint
    template<class mint>
    vector<mint> mul(vector<mint> A,vector<mint> B) {
        if (A.empty() || B.empty()) return {};
        int n = int(A.size()), m = int(B.size());
        if (min(n, m) < 30) return naive_mul(A, B);
        int MOD = A[0].getmod();
        int z = 1 << ceil_pow2(n + m - 1);
        if (MOD == 998244353) {
          A.resize(z);
          trans(A);
          B.resize(z);
          trans(B);
          for (int i = 0; i < z; i++) {
            A[i] *= B[i];
          }
          trans_inv(A);
          A.resize(n + m - 1);
          mint iz = modinv(mint(z));
          for (int i = 0; i < n + m - 1; i++) A[i] *= iz;
          return A;
        }
        vector<mint0> a0(z, 0), b0(z, 0);
        vector<mint1> a1(z, 0), b1(z, 0);
        vector<mint2> a2(z, 0), b2(z, 0);
        for (int i = 0; i < n; ++i)
            a0[i] = A[i].val, a1[i] = A[i].val, a2[i] = A[i].val;
        for (int i = 0; i < m; ++i)
            b0[i] = B[i].val, b1[i] = B[i].val, b2[i] = B[i].val;
        trans(a0), trans(a1), trans(a2), trans(b0), trans(b1), trans(b2);
        for (int i = 0; i < z; ++i) {
            a0[i] *= b0[i];
            a1[i] *= b1[i];
            a2[i] *= b2[i];
        }
        trans_inv(a0), trans_inv(a1), trans_inv(a2);
        static const mint mod0 = MOD0, mod01 = mod0 * MOD1;
        mint0 i0=modinv(mint0(z));
        mint1 i1=modinv(mint1(z));
        mint2 i2=modinv(mint2(z));
        vector<mint> res(n + m - 1);
        for (int i = 0; i < n + m - 1; ++i) {
            a0[i]*=i0;
            a1[i]*=i1;
            a2[i]*=i2;
            int y0 = a0[i].val;
            int y1 = (imod0 * (a1[i] - y0)).val;
            int y2 = (imod01 * (a2[i] - y0) - imod1 * y1).val;
            res[i] = mod01 * y2 + mod0 * y1 + y0;
        }
        return res;
    }
};
// Formal Power Series
template <typename mint> struct FPS : vector<mint> {
    using vector<mint>::vector;
 /*
    template<class...Args>
    FPS(Args...args) : vector<mint>(args...){}
  */
    // constructor
    FPS(const vector<mint>& r) : vector<mint>(r) {}
 
    // core operator
    inline FPS pre(int siz) const {
        return FPS(begin(*this), begin(*this) + min((int)this->size(), siz));
    }
    inline FPS rev() const {
        FPS res = *this;
        reverse(begin(res), end(res));
        return res;
    }
    inline FPS& normalize() {
        while (!this->empty() && this->back() == 0) this->pop_back();
        return *this;
    }
 
    // basic operator
    inline FPS operator - () const noexcept {
        FPS res = (*this);
        for (int i = 0; i < (int)res.size(); ++i) res[i] = -res[i];
        return res;
    }
    inline FPS operator + (const mint& v) const { return FPS(*this) += v; }
    inline FPS operator + (const FPS& r) const { return FPS(*this) += r; }
    inline FPS operator - (const mint& v) const { return FPS(*this) -= v; }
    inline FPS operator - (const FPS& r) const { return FPS(*this) -= r; }
    inline FPS operator * (const mint& v) const { return FPS(*this) *= v; }
    inline FPS operator * (const FPS& r) const { return FPS(*this) *= r; }
    inline FPS operator / (const mint& v) const { return FPS(*this) /= v; }
    inline FPS operator << (int x) const { return FPS(*this) <<= x; }
    inline FPS operator >> (int x) const { return FPS(*this) >>= x; }
    inline FPS& operator += (const mint& v) {
        if (this->empty()) this->resize(1);
        (*this)[0] += v;
        return *this;
    }
    inline FPS& operator += (const FPS& r) {
        if (r.size() > this->size()) this->resize(r.size());
        for (int i = 0; i < (int)r.size(); ++i) (*this)[i] += r[i];
        return this->normalize();
    }
    inline FPS& operator -= (const mint& v) {
        if (this->empty()) this->resize(1);
        (*this)[0] -= v;
        return *this;
    }
    inline FPS& operator -= (const FPS& r) {
        if (r.size() > this->size()) this->resize(r.size());
        for (int i = 0; i < (int)r.size(); ++i) (*this)[i] -= r[i];
        return this->normalize();
    }
    inline FPS& operator *= (const mint& v) {
        for (int i = 0; i < (int)this->size(); ++i) (*this)[i] *= v;
        return *this;
    }
    inline FPS& operator *= (const FPS& r) {
        return *this = NTT::mul((*this), r);
    }
    inline FPS& operator /= (const mint& v) {
        assert(v != 0);
        mint iv = modinv(v);
        for (int i = 0; i < (int)this->size(); ++i) (*this)[i] *= iv;
        return *this;
    }
    inline FPS& operator <<= (int x) {
        FPS res(x, 0);
        res.insert(res.end(), begin(*this), end(*this));
        return *this = res;
    }
    inline FPS& operator >>= (int x) {
        FPS res;
        res.insert(res.end(), begin(*this) + x, end(*this));
        return *this = res;
    }
    inline mint eval(const mint& v){
        mint res = 0;
        for (int i = (int)this->size()-1; i >= 0; --i) {
            res *= v;
            res += (*this)[i];
        }
        return res;
    }
    inline friend FPS gcd(const FPS& f, const FPS& g) {
        if (g.empty()) return f;
        return gcd(g, f % g);
    }

    // advanced operation
    // df/dx
    inline friend FPS diff(const FPS& f) {
        int n = (int)f.size();
        FPS res(n-1);
        for (int i = 1; i < n; ++i) res[i-1] = f[i] * i;
        return res;
    }

    // \int f dx
    inline friend FPS integral(const FPS& f) {
        int n = (int)f.size();
        FPS res(n+1, 0);
        for (int i = 0; i < n; ++i) res[i+1] = f[i] / (i+1);
        return res;
    }

    // inv(f), f[0] must not be 0
    inline friend FPS inv(const FPS& f, int deg) {
        assert(f[0] != 0);
        if (deg < 0) deg = (int)f.size();
        FPS res({mint(1) / f[0]});
        for (int i = 1; i < deg; i <<= 1) {
            res = (res + res - res * res * f.pre(i << 1)).pre(i << 1);
        }
        res.resize(deg);
        return res;
    }
    inline friend FPS inv(const FPS& f) {
        return inv(f, f.size());
    }

    // division, r must be normalized (r.back() must not be 0)
    inline FPS& operator /= (const FPS& r) {
        const int n=(*this).size(),m=r.size();
        if(n<m){
            (*this).clear();
            return *this;
        }
        assert(r.back() != 0);
        this->normalize();
        if (this->size() < r.size()) {
            this->clear();
            return *this;
        }
        int need = (int)this->size() - (int)r.size() + 1;
        *this = ((*this).rev().pre(need) * inv(r.rev(), need)).pre(need).rev();
        return *this;
    }
    inline FPS& operator %= (const FPS &r) {
        const int n=(*this).size(),m=r.size();
        if(n<m) return (*this);
        assert(r.back() != 0);
        this->normalize();
        FPS q = (*this) / r;
        return *this -= q * r;
    }
    inline FPS operator / (const FPS& r) const { return FPS(*this) /= r; }
    inline FPS operator % (const FPS& r) const { return FPS(*this) %= r; }

    // log(f) = \int f'/f dx, f[0] must be 1
    inline friend FPS log(const FPS& f, int deg) {
        assert(f[0] == 1);
        FPS res = integral((diff(f) * inv(f, deg)).pre(deg-1));
        return res;
    }
    inline friend FPS log(const FPS& f) {
        return log(f, f.size());
    }

    // exp(f), f[0] must be 0
    inline friend FPS exp(const FPS& f, int deg) {
        assert(f[0] == 0);
        FPS res(1, 1);
        for (int i = 1; i < deg; i <<= 1) {
            res = res * (f.pre(i<<1) - log(res, i<<1) + 1).pre(i<<1);
        }
        res.resize(deg);
        return res;
    }
    inline friend FPS exp(const FPS& f) {
        return exp(f, f.size());
    }

    // pow(f) = exp(e * log f)
    inline friend FPS pow(const FPS& f, long long e, int deg) {
        long long i = 0;
        while (i < (int)f.size() && f[i] == 0) ++i;
        if (i == (int)f.size()) return FPS(deg, 0);
        if (i * e >= deg) return FPS(deg, 0);
        mint k = f[i];
        FPS res = exp(log((f >> i) / k, deg) * e, deg) * modpow(k, e) << (e * i);
        res.resize(deg);
        return res;
    }
    inline friend FPS pow(const FPS& f, long long e) {
        return pow(f, e, f.size());
    }

    // sqrt(f), f[0] must be 1
    inline friend FPS sqrt_base(const FPS& f, int deg) {
        assert(f[0] == 1);
        mint inv2 = mint(1) / 2;
        FPS res(1, 1);
        for (int i = 1; i < deg; i <<= 1) {
            res = (res + f.pre(i << 1) * inv(res, i << 1)).pre(i << 1);
            for (mint& x : res) x *= inv2;
        }
        res.resize(deg);
        return res;
    }
    inline friend FPS sqrt_base(const FPS& f) {
        return sqrt_base(f, f.size());
    }
    FPS taylor_shift(mint c) const {
      int n = (int) this->size();
      vector<mint> fact(n), rfact(n);
      fact[0] = rfact[0] = mint(1);
      for(int i = 1; i < n; i++) fact[i] = fact[i - 1] * mint(i);
      rfact[n - 1] = mint(1) / fact[n - 1];
      for(int i = n - 1; i > 1; i--) rfact[i - 1] = rfact[i] * mint(i);
      FPS p(*this);
      for(int i = 0; i < n; i++) p[i] *= fact[i];
      p = p.rev();
      FPS bs(n, mint(1));
      for(int i = 1; i < n; i++) bs[i] = bs[i - 1] * c * rfact[i] * fact[i - 1];
      p = (p * bs).pre(n);
      p = p.rev();
      for(int i = 0; i < n; i++) p[i] *= rfact[i];
      return p;
    }
};
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
vector<int> topo_sort(const vector<vector<int>> &G) { 
  vector<int> ret;
  int n = (int)G.size();
  vector<int> ind(n);      
  for (int i = 0; i < n; i++) {  
    for (auto e : G[i]) {
      ind[e]++;
    }
  }
  priority_queue<int,vector<int>,greater<int>> que;
  for (int i = 0; i < n; i++) {
    if (ind[i] == 0) {
      que.push(i);
    }
  }
  while (!que.empty()) { 
    int now = que.top();
    ret.push_back(now);
    que.pop();
    for (auto e : G[now]) {
      ind[e]--;
      if (ind[e] == 0) {
        que.push(e);
      }
    }
  }
  return ret;
}
template <typename mint> FPS<mint> product(vector<FPS<mint>> a){
  int siz=1;
  while(siz<int(a.size())){
    siz<<=1;
  }
  vector<FPS<mint>> res(siz*2-1,{1});
  for(int i=0;i<int(a.size());++i){
    res[i+siz-1]=a[i];
  }
  for(int i=siz-2;i>=0;--i){
    res[i]=res[2*i+1]*res[2*i+2];
  }
  return res[0];
}
template<typename mint> FPS<mint> inv_sum(vector<FPS<mint>> f){
  int siz=1;
  while(siz<int(f.size())){
    siz<<=1;
  }
  vector<FPS<mint>> mol(siz*2-1),dem(siz*2-1,{1});
  for(size_t i=0;i<f.size();++i){
    mol[i+siz-1]={1};
    dem[i+siz-1]=f[i];
  }
  for(int i=siz-2;i>=0;--i){
    dem[i]=dem[2*i+1]*dem[2*i+2];
    mol[i]=mol[2*i+1]*dem[2*i+2]+mol[2*i+2]*dem[2*i+1];
  }
  mol[0]*=inv(dem[0]);
  return mol[0];
}
template <typename mint> FPS<mint> rev(FPS<mint> p) { reverse(p.begin(),p.end()); return p; }
template <typename mint> FPS<mint> RSZ(FPS<mint> p, int x) { p.resize(x); return p; }
template<int m> struct Perm{
  unordered_map<int,tuple<ll,vector<ll>,vector<ll>>> mp;
  int n_;
  ll p_, pm_;
  vector<ll> ord_, fact_;
  vector<P> pf;
  Perm(int n) : n_(n), ord_(n), fact_(n) { 
    pf=prime_factorize(m); 
    PERMinit();
  }
  void init(int n) {
    ord_.resize(n);
    fact_.resize(n);
  }
  void init(long long p, long long pm) {
    p_=p,pm_=pm;
    ord_[0]=ord_[1]=0;
    fact_[0]=fact_[1]=1;
    auto&[pms,ord,fac]=mp[p];
    pms=pm;
    ord.resize(n_);
    fac.resize(n_);
    ord[0]=ord[1]=0;
    fac[0]=fac[1]=1;
    for (int i=2;i<n_;i++) {
      long long add=0;
      long long ni=i;
      while (ni % p == 0) add++,ni/=p;
      ord_[i]=ord_[i-1]+add;
      fact_[i]=fact_[ni-1]*ni%pm;
      ord[i]=ord_[i];
      fac[i]=fact_[i];
    }
  }
  void init(long long p, long long pm, int n) {
    init(n);
    init(p, pm);
  }
  void PERMinit(){
    for(auto p : pf){
      ll ps=p.first,pfs=pow(p.first,p.second);
      init(n_);
      init(ps,pfs);
    }
  }
  ll perm(ll n, ll r,int p) {
    if (n<0 || r<0 || n<r) return 0;
    auto&[pms,ord,fac]=mp[p];
    ll e=ord[n]-ord[n-r];
    ll res=fac[n]*modinv(fac[n-r]%pms,pms)%pms;
    res=res*modpow(p,e,pms)%pms;
    return res;
  }
  ll operator()(int n, int k){
    if(n<0 || k<0 || n<k) return 0;
    vector<long long> vb, vm;
    for (auto ps : pf) {
        long long p = ps.first, e = ps.second;
        long long pm = pow(p,e);
        long long b = 1;
        b *= perm(n, k ,p) % pm;
        b %= pm;
        vm.push_back(pm);
        vb.push_back(b);
    }
    auto res = ChineseRem(vb,vm);
    return res.first;
  }
};
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
template <typename T>
vector<int> compress(vector<T> &x){
  vector<T> vals=x;
  vector<int> res;
  sort(vals.begin(),vals.end());
  auto it=unique(vals.begin(),vals.end());
  vals.erase(it,vals.end());
  for(auto xx : x){
    int ret=lower_bound(vals.begin(),vals.end(),xx) - vals.begin();
    res.push_back(ret);
  }
  return res;
}
/**
 * @brief Succinct Indexable Dictionary(完備辞書)
 */
struct SuccinctIndexableDictionary {
  size_t length;
  size_t blocks;
  vector< unsigned > bit, sum;

  SuccinctIndexableDictionary() = default;

  SuccinctIndexableDictionary(size_t length) : length(length), blocks((length + 31) >> 5) {
    bit.assign(blocks, 0U);
    sum.assign(blocks, 0U);
  }

  void set(int k) {
    bit[k >> 5] |= 1U << (k & 31);
  }

  void build() {
    sum[0] = 0U;
    for(int i = 1; i < blocks; i++) {
      sum[i] = sum[i - 1] + __builtin_popcount(bit[i - 1]);
    }
  }

  bool operator[](int k) {
    return (bool((bit[k >> 5] >> (k & 31)) & 1));
  }

  int rank(int k) {
    return (sum[k >> 5] + __builtin_popcount(bit[k >> 5] & ((1U << (k & 31)) - 1)));
  }

  int rank(bool val, int k) {
    return (val ? rank(k) : k - rank(k));
  }
};

/*
 * @brief Wavelet Matrix(ウェーブレット行列)
 * @docs docs/wavelet-matrix.md
 */
template< typename T, int MAXLOG >
struct WaveletMatrix {
  size_t length;
  SuccinctIndexableDictionary matrix[MAXLOG];
  int mid[MAXLOG];

  WaveletMatrix() = default;

  WaveletMatrix(vector< T > v) : length(v.size()) {
    vector< T > l(length), r(length);
    for(int level = MAXLOG - 1; level >= 0; level--) {
      matrix[level] = SuccinctIndexableDictionary(length + 1);
      int left = 0, right = 0;
      for(int i = 0; i < length; i++) {
        if(((v[i] >> level) & 1)) {
          matrix[level].set(i);
          r[right++] = v[i];
        } else {
          l[left++] = v[i];
        }
      }
      mid[level] = left;
      matrix[level].build();
      v.swap(l);
      for(int i = 0; i < right; i++) {
        v[left + i] = r[i];
      }
    }
  }

  pair< int, int > succ(bool f, int l, int r, int level) {
    return {matrix[level].rank(f, l) + mid[level] * f, matrix[level].rank(f, r) + mid[level] * f};
  }

  // v[k]
  T access(int k) {
    T ret = 0;
    for(int level = MAXLOG - 1; level >= 0; level--) {
      bool f = matrix[level][k];
      if(f) ret |= T(1) << level;
      k = matrix[level].rank(f, k) + mid[level] * f;
    }
    return ret;
  }

  T operator[](const int &k) {
    return access(k);
  }

  // count i s.t. (0 <= i < r) && v[i] == x
  int rank(const T &x, int r) {
    int l = 0;
    for(int level = MAXLOG - 1; level >= 0; level--) {
      tie(l, r) = succ((x >> level) & 1, l, r, level);
    }
    return r - l;
  }

  // k-th(0-indexed) smallest number in v[l,r)
  T kth_smallest(int l, int r, int k) {
    assert(0 <= k && k < r - l);
    T ret = 0;
    for(int level = MAXLOG - 1; level >= 0; level--) {
      int cnt = matrix[level].rank(false, r) - matrix[level].rank(false, l);
      bool f = cnt <= k;
      if(f) {
        ret |= T(1) << level;
        k -= cnt;
      }
      tie(l, r) = succ(f, l, r, level);
    }
    return ret;
  }

  // k-th(0-indexed) largest number in v[l,r)
  T kth_largest(int l, int r, int k) {
    return kth_smallest(l, r, r - l - k - 1);
  }

  // count i s.t. (l <= i < r) && (v[i] < upper)
  int range_freq(int l, int r, T upper) {
    int ret = 0;
    for(int level = MAXLOG - 1; level >= 0; level--) {
      bool f = ((upper >> level) & 1);
      if(f) ret += matrix[level].rank(false, r) - matrix[level].rank(false, l);
      tie(l, r) = succ(f, l, r, level);
    }
    return ret;
  }

  // count i s.t. (l <= i < r) && (lower <= v[i] < upper)
  int range_freq(int l, int r, T lower, T upper) {
    return range_freq(l, r, upper) - range_freq(l, r, lower);
  }

  // max v[i] s.t. (l <= i < r) && (v[i] < upper)
  T prev_value(int l, int r, T upper) {
    int cnt = range_freq(l, r, upper);
    return cnt == 0 ? T(-1) : kth_smallest(l, r, cnt - 1);
  }

  // min v[i] s.t. (l <= i < r) && (lower <= v[i])
  T next_value(int l, int r, T lower) {
    int cnt = range_freq(l, r, lower);
    return cnt == r - l ? T(-1) : kth_smallest(l, r, cnt);
  }
};

template< typename T, int MAXLOG >
struct CompressedWaveletMatrix {
  WaveletMatrix< int, MAXLOG > mat;
  vector< T > ys;

  CompressedWaveletMatrix(const vector< T > &v) : ys(v) {
    sort(begin(ys), end(ys));
    ys.erase(unique(begin(ys), end(ys)), end(ys));
    vector< int > t(v.size());
    for(int i = 0; i < v.size(); i++) t[i] = get(v[i]);
    mat = WaveletMatrix< int, MAXLOG >(t);
  }

  inline int get(const T& x) {
    return lower_bound(begin(ys), end(ys), x) - begin(ys);
  }

  T access(int k) {
    return ys[mat.access(k)];
  }

  T operator[](const int &k) {
    return access(k);
  }

  int rank(const T &x, int r) {
    auto pos = get(x);
    if(pos == ys.size() || ys[pos] != x) return 0;
    return mat.rank(pos, r);
  }

  T kth_smallest(int l, int r, int k) {
    return ys[mat.kth_smallest(l, r, k)];
  }

  T kth_largest(int l, int r, int k) {
    return ys[mat.kth_largest(l, r, k)];
  }

  int range_freq(int l, int r, T upper) {
    return mat.range_freq(l, r, get(upper));
  }

  int range_freq(int l, int r, T lower, T upper) {
    return mat.range_freq(l, r, get(lower), get(upper));
  }

  T prev_value(int l, int r, T upper) {
    auto ret = mat.prev_value(l, r, get(upper));
    return ret == -1 ? T(-1) : ys[ret];
  }

  T next_value(int l, int r, T lower) {
    auto ret = mat.next_value(l, r, get(lower));
    return ret == -1 ? T(-1) : ys[ret];
  }
};
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
template<typename mint>
struct subproduct_tree{
  using poly=FPS<mint>;
  vector<poly> tree;
  int n,siz;
  subproduct_tree(const vector<mint> &x){
    n=1;
    siz=sz(x);
    while(n<siz) n*=2;;
    tree.resize(2*n,{mint(1)});
    for(int i=0;i<siz;i++) tree[i+n]={-x[i],mint(1)};
    for(int i=n-1;i>0;i--) tree[i]=tree[2*i]*tree[2*i+1];
  }
  vector<mint> multieval(const poly &f){
    vector<poly> remainder(2*n);
    remainder[1]=f%tree[1];
    for(int i=1;i<n;i++){
      remainder[2*i]=remainder[i]%tree[2*i];
      remainder[2*i+1]=remainder[i]%tree[2*i+1];
    }
    vector<mint> ret(siz);
    for(int i=0;i<siz;i++){
        if(remainder[i+n].empty()) ret[i]=0;
        else ret[i]=remainder[i+n][0];
    }
    return ret;
  }
  poly interpolate(const vector<mint> &y){
    poly g=diff(tree[1]);
    vector<mint> evaled=multieval(g);
    vector<poly> mol(2*n),dem(2*n,{1});
    for(int i=0;i<siz;++i){
      mol[i+n]={y[i]};
      dem[i+n]=tree[i+n]*evaled[i];
    }
    for(int i=n-1;i>0;--i){
      dem[i]=dem[2*i]*dem[2*i+1];
      mol[i]=mol[2*i]*dem[2*i+1]+mol[2*i+1]*dem[2*i];
    }
    mol[1]*=inv(dem[1]);
    return RSZ(tree[1]*mol[1],siz);
  }
};
template <typename mint> vector<mint> multieval(const FPS<mint> &f,const vector<mint> &x){
  subproduct_tree<mint> tree(x);
  return tree.multieval(f);
}
template <typename mint> FPS<mint> interpolate(const vector<mint> &x,const vector<mint> &y){
  assert(sz(x)==sz(y));
  if(sz(x)==1) return {y[0]};
  subproduct_tree<mint> tree(x);
  return tree.interpolate(y);
}
template <typename mint> struct PartitionNumber{
    using poly=FPS<mint>;
    poly f;
    PartitionNumber(int maxN,vector<int> st) : f(maxN+1){
        vector<mint> invs(maxN+1);
        for(int i=1;i<=maxN;i++) invs[i]=modinv(mint(i));
        for(auto& e:st){
            for(int j=1;j<=maxN/e;j++){
                f[e*j]+=invs[j];
            }
        }
        f=exp(f);
    }
    PartitionNumber(int maxN) : f(maxN+1){
        f[0]=mint(1);
        for(int k=1;k<=maxN;k++){
            if(1ll*k*(3*k+1)/2<=maxN) f[k*(3*k+1)/2]+=(k % 2 ? mint(-1) : mint(1));
            if(1ll*k*(3*k-1)/2<=maxN) f[k*(3*k-1)/2]+=(k % 2 ? mint(-1) : mint(1));
        }
        f=inv(f);
    }
    mint operator()(int n){
        return f[n];
    }
};
template< typename T >
struct Combination {
  vector< T > _fact, _rfact, _inv;

  Combination(int sz) : _fact(sz + 1), _rfact(sz + 1), _inv(sz + 1) {
    _fact[0] = _rfact[sz] = _inv[0] = 1;
    for(int i = 1; i <= sz; i++) _fact[i] = _fact[i - 1] * i;
    _rfact[sz] /= _fact[sz];
    for(int i = sz - 1; i >= 0; i--) _rfact[i] = _rfact[i + 1] * (i + 1);
    for(int i = 1; i <= sz; i++) _inv[i] = _rfact[i] * _fact[i - 1];
  }

  inline T fact(int k) const { return _fact[k]; }

  inline T rfact(int k) const { return _rfact[k]; }

  inline T inv(int k) const { return _inv[k]; }

  T P(int n, int r) const {
    if(r < 0 || n < r) return 0;
    return fact(n) * rfact(n - r);
  }

  T C(int p, int q) const {
    if(q < 0 || p < q) return 0;
    return fact(p) * rfact(q) * rfact(p - q);
  }

  T H(int n, int r) const {
    if(n < 0 || r < 0) return (0);
    return r == 0 ? 1 : C(n + r - 1, r);
  }
};
template< typename T >
T lagrange_poly(const vector< T > &y,ll t) {
  int N = y.size() - 1;
  Combination< T > comb(N);
  if(t <= N) return y[t];
  T ret(0);
  vector< T > dp(N + 1, 1), pd(N + 1, 1);
  for(int i = 0; i < N; i++) dp[i + 1] = dp[i] * (t - i);
  for(int i = N; i > 0; i--) pd[i - 1] = pd[i] * (t - i);
  for(int i = 0; i <= N; i++) {
    T tmp = y[i] * dp[i] * pd[i] * comb.rfact(i) * comb.rfact(N - i);
    if((N - i) & 1) ret -= tmp;
    else ret += tmp;
  }
  return ret;
}
//fast Input by yosupo
#include <unistd.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cstring>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>
namespace fastio{
/*
  quote from yosupo's submission in Library Checker
*/
int bsr(unsigned int n) {
    return 8 * (int)sizeof(unsigned int) - 1 - __builtin_clz(n);
}
// @param n `1 <= n`
// @return maximum non-negative `x` s.t. `(n & (1 << x)) != 0`
int bsr(unsigned long n) {
    return 8 * (int)sizeof(unsigned long) - 1 - __builtin_clzl(n);
}
// @param n `1 <= n`
// @return maximum non-negative `x` s.t. `(n & (1 << x)) != 0`
int bsr(unsigned long long n) {
    return 8 * (int)sizeof(unsigned long long) - 1 - __builtin_clzll(n);
}
// @param n `1 <= n`
// @return minimum non-negative `x` s.t. `(n & (1 << x)) != 0`
int bsr(unsigned __int128 n) {
    unsigned long long low = (unsigned long long)(n);
    unsigned long long high = (unsigned long long)(n >> 64);
    return high ? 127 - __builtin_clzll(high) : 63 - __builtin_ctzll(low);
}
 
namespace internal {
 
template <class T>
using is_signed_int128 =
    typename std::conditional<std::is_same<T, __int128_t>::value ||
                                  std::is_same<T, __int128>::value,
                              std::true_type,
                              std::false_type>::type;
 
template <class T>
using is_unsigned_int128 =
    typename std::conditional<std::is_same<T, __uint128_t>::value ||
                                  std::is_same<T, unsigned __int128>::value,
                              std::true_type,
                              std::false_type>::type;
 
template <class T>
using make_unsigned_int128 =
    typename std::conditional<std::is_same<T, __int128_t>::value,
                              __uint128_t,
                              unsigned __int128>;
 
template <class T>
using is_integral =
    typename std::conditional<std::is_integral<T>::value ||
                                  internal::is_signed_int128<T>::value ||
                                  internal::is_unsigned_int128<T>::value,
                              std::true_type,
                              std::false_type>::type;
 
template <class T>
using is_signed_int = typename std::conditional<(is_integral<T>::value &&
                                                 std::is_signed<T>::value) ||
                                                    is_signed_int128<T>::value,
                                                std::true_type,
                                                std::false_type>::type;
 
template <class T>
using is_unsigned_int =
    typename std::conditional<(is_integral<T>::value &&
                               std::is_unsigned<T>::value) ||
                                  is_unsigned_int128<T>::value,
                              std::true_type,
                              std::false_type>::type;
 
template <class T>
using to_unsigned = typename std::conditional<
    is_signed_int128<T>::value,
    make_unsigned_int128<T>,
    typename std::conditional<std::is_signed<T>::value,
                              std::make_unsigned<T>,
                              std::common_type<T>>::type>::type;
 
template <class T>
using is_integral_t = std::enable_if_t<is_integral<T>::value>;
 
template <class T>
using is_signed_int_t = std::enable_if_t<is_signed_int<T>::value>;
 
template <class T>
using is_unsigned_int_t = std::enable_if_t<is_unsigned_int<T>::value>;
 
template <class T> using to_unsigned_t = typename to_unsigned<T>::type;
 
}  // namespace internal
struct Scanner {
  public:
    Scanner(const Scanner&) = delete;
    Scanner& operator=(const Scanner&) = delete;
 
    Scanner(FILE* fp) : fd(fileno(fp)) {}
 
    void read() {}
    template <class H, class... T> void read(H& h, T&... t) {
        bool f = read_single(h);
        assert(f);
        read(t...);
    }
 
    int read_unsafe() { return 0; }
    template <class H, class... T> int read_unsafe(H& h, T&... t) {
        bool f = read_single(h);
        if (!f) return 0;
        return 1 + read_unsafe(t...);
    }
 
    int close() { return ::close(fd); }
 
  private:
    static constexpr int SIZE = 1 << 15;
 
    int fd = -1;
    std::array<char, SIZE + 1> line;
    int st = 0, ed = 0;
    bool eof = false;
 
    bool read_single(std::string& ref) {
        if (!skip_space()) return false;
        ref = "";
        while (true) {
            char c = top();
            if (c <= ' ') break;
            ref += c;
            st++;
        }
        return true;
    }
    bool read_single(double& ref) {
        std::string s;
        if (!read_single(s)) return false;
        ref = std::stod(s);
        return true;
    }
 
    template <class T,
              std::enable_if_t<std::is_same<T, char>::value>* = nullptr>
    bool read_single(T& ref) {
        if (!skip_space<50>()) return false;
        ref = top();
        st++;
        return true;
    }
 
    template <class T,
              internal::is_signed_int_t<T>* = nullptr,
              std::enable_if_t<!std::is_same<T, char>::value>* = nullptr>
    bool read_single(T& sref) {
        using U = internal::to_unsigned_t<T>;
        if (!skip_space<50>()) return false;
        bool neg = false;
        if (line[st] == '-') {
            neg = true;
            st++;
        }
        U ref = 0;
        do {
            ref = 10 * ref + (line[st++] & 0x0f);
        } while (line[st] >= '0');
        sref = neg ? -ref : ref;
        return true;
    }
    template <class U,
              internal::is_unsigned_int_t<U>* = nullptr,
              std::enable_if_t<!std::is_same<U, char>::value>* = nullptr>
    bool read_single(U& ref) {
        if (!skip_space<50>()) return false;
        ref = 0;
        do {
            ref = 10 * ref + (line[st++] & 0x0f);
        } while (line[st] >= '0');
        return true;
    }
 
    bool reread() {
        if (ed - st >= 50) return true;
        if (st > SIZE / 2) {
            std::memmove(line.data(), line.data() + st, ed - st);
            ed -= st;
            st = 0;
        }
        if (eof) return false;
        auto u = ::read(fd, line.data() + ed, SIZE - ed);
        if (u == 0) {
            eof = true;
            line[ed] = '\0';
            u = 1;
        }
        ed += int(u);
        line[ed] = char(127);
        return true;
    }
 
    char top() {
        if (st == ed) {
            bool f = reread();
            assert(f);
        }
        return line[st];
    }
 
    template <int TOKEN_LEN = 0>
    bool skip_space() {
        while (true) {
            while (line[st] <= ' ') st++;   
            if (ed - st > TOKEN_LEN) return true;
            if (st > ed) st = ed;
            for (auto i = st; i < ed; i++) {
                if (line[i] <= ' ') return true;
            }
            if (!reread()) return false;
        }
    }
};

//fast Output by ei1333
/**
 * @brief Printer(高速出力)
 */
struct Printer {
public:
  explicit Printer(FILE *fp) : fp(fp) {}

  ~Printer() { flush(); }

  template< bool f = false, typename T, typename... E >
  void write(const T &t, const E &... e) {
    if(f) write_single(' ');
    write_single(t);
    write< true >(e...);
  }

  template< typename... T >
  void writeln(const T &...t) {
    write(t...);
    write_single('\n');
  }

  void flush() {
    fwrite(line, 1, st - line, fp);
    st = line;
  }

private:
  FILE *fp = nullptr;
  static constexpr size_t line_size = 1 << 16;
  static constexpr size_t int_digits = 20;
  char line[line_size + 1] = {};
  char small[32] = {};
  char *st = line;

  template< bool f = false >
  void write() {}

  void write_single(const char &t) {
    if(st + 1 >= line + line_size) flush();
    *st++ = t;
  }

  template< typename T, enable_if_t< is_integral< T >::value, int > = 0 >
  void write_single(T s) {
    if(st + int_digits >= line + line_size) flush();
    if(s == 0) {
      write_single('0');
      return;
    }
    if(s < 0) {
      write_single('-');
      s = -s;
    }
    char *mp = small + sizeof(small);
    typename make_unsigned< T >::type y = s;
    size_t len = 0;
    while(y > 0) {
      *--mp = y % 10 + '0';
      y /= 10;
      ++len;
    }
    memmove(st, mp, len);
    st += len;
  }

  void write_single(const string &s) {
    for(auto &c : s) write_single(c);
  }

  void write_single(const char *s) {
    while(*s != 0) write_single(*s++);
  }

  template< typename T >
  void write_single(const vector< T > &s) {
    for(size_t i = 0; i < s.size(); i++) {
      if(i) write_single(' ');
      write_single(s[i]);
    }
  }
};

}; //namespace fastio
using mint=Fp<998244353>;
int main(){
    fastio::Scanner sc(stdin);
    fastio::Printer pr(stdout);
}
