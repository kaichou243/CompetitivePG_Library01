#include<bits/stdc++.h>
#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
using namespace std;
using ll=long long;
using FLOW=int;
using P = pair<ll,ll>;
const long double PI=acos(-1);
const ll INF=1e18;
const int inf=1e9;
struct Edge {
    ll to;
    ll cost;
};
struct edge {
    int rev, from, to;
    FLOW cap, icap;
    edge(int r, int f, int t, FLOW c) : rev(r), from(f), to(t), cap(c), icap(c) {}
    friend ostream& operator << (ostream& s, const edge& E) {
        if (E.cap > 0) return s << E.from << "->" << E.to << '(' << E.cap << ')';
        else return s;
    }
};
using Graph=vector<vector<Edge>>;
struct MaxFlowGraph{
    int V;
    vector<vector<edge>> list;
 
    MaxFlowGraph(int n = 0) : V(n), list(n) { for (int i = 0; i < n; ++i) list[i].clear(); }
    void init(int n = 0) { V = n; list.resize(n); for (int i = 0; i < n; ++i) list[i].clear(); }
    void resize(int n = 0) { V = n; }
    void reset() { for (int i = 0; i < V; ++i) for (int j = 0; j < list[i].size(); ++j) list[i][j].cap = list[i][j].icap; }
    inline vector<edge>& operator [] (int i) { return list[i]; }
 
    edge &redge(edge e) {
        if (e.from != e.to) return list[e.to][e.rev];
        else return list[e.to][e.rev + 1];
    }
 
    void addedge(int from, int to, FLOW cap) {
        list[from].push_back(edge((int)list[to].size(), from, to, cap));
        list[to].push_back(edge((int)list[from].size() - 1, to, from, 0));
    }
};
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
  constexpr int getmod() const { return MOD; }
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
    if (n == 0) return 1;
    auto t = modpow(a, n / 2);
    t = t * t;
    if (n & 1) t = t * a;
    return t;
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
  void addedge(int u, int v) {
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
struct SegmentTree{
  ll g(ll a,ll b){
    return a+b; //g(a,b) : aとbで二項演算を行った結果を返す関数
  }
  int n;
  vector<ll> dat;
  SegmentTree(int N){
    n=1;
    while(n<N)n*=2;
    dat.assign(2*n-1,0);
  }
  void add(int k,ll x){
    k+=n-1;
    dat[k]+=x;
    while(k){
      k=(k-1)/2;
      dat[k]=g(dat[k*2+1],dat[k*2+2]);
    }
  }
  void update(int k,ll x){
    k+=n-1;
    dat[k]=x;
    while(k){
      k=(k-1)/2;
      dat[k]=g(dat[k*2+1],dat[k*2+2]);
    }
  }
  ll query(int l,int r){
    return query_sub(l,r,0,0,n); //区間[l,r)に対する二項演算の総積を返す
  }
  ll query_sub(int l,int r,int k=0,int a=0,int b=-1){
    if(b<0)b=n;
    if(r<=a || b<=l)return 0; //使う二項演算に合わせて単位元を変える。要注意!
    if(l<=a && b<=r)return dat[k];
    return g(query_sub(l,r,k*2+1,a,(a+b)/2),query_sub(l,r,k*2+2,(a+b)/2,b));
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
template <typename X, typename M>
struct LazySegTree {
    using FX = function<X(X, X)>;
    using FA = function<X(X, M)>;
    using FM = function<M(M, M)>;
    using FP = function<M(M, int)>;
    int n;
    FX fx;
    FA fa;
    FM fm;
    FP fp;
    const X ex;
    const M em;
    vector<X> dat;
    vector<M> lazy;
    LazySegTree(int n_, FX fx_, FA fa_, FM fm_, FP fp_, X ex_, M em_)
        : n(), fx(fx_), fa(fa_), fm(fm_), fp(fp_), ex(ex_), em(em_), dat(n_ * 4, ex), lazy(n_ * 4, em) {
        int x = 1;
        while (n_ > x) x *= 2;
        n = x;
    }
 
    void set(int i, X x) { dat[i + n - 1] = x; }
    void build() {
        for (int k = n - 2; k >= 0; k--) dat[k] = fx(dat[2 * k + 1], dat[2 * k + 2]);
    }
 
    /* lazy eval */
    void eval(int k, int len) {
        if (lazy[k] == em) return;  // 更新するものが無ければ終了
        if (k < n - 1) {            // 葉でなければ子に伝搬
            lazy[k * 2 + 1] = fm(lazy[k * 2 + 1], lazy[k]);
            lazy[k * 2 + 2] = fm(lazy[k * 2 + 2], lazy[k]);
        }
        // 自身を更新
        dat[k] = fa(dat[k], fp(lazy[k], len));
        lazy[k] = em;
    }
 
    void update(int a, int b, M x, int k, int l, int r) {
        eval(k, r - l);
        if (a <= l && r <= b) {  // 完全に内側の時
            lazy[k] = fm(lazy[k], x);
            eval(k, r - l);
        } else if (a < r && l < b) {                     // 一部区間が被る時
            update(a, b, x, k * 2 + 1, l, (l + r) / 2);  // 左の子
            update(a, b, x, k * 2 + 2, (l + r) / 2, r);  // 右の子
            dat[k] = fx(dat[k * 2 + 1], dat[k * 2 + 2]);
        }
    }
    //[a,b)
    void update(int a, int b, M x) { update(a, b, x, 0, 0, n); }
    X query_sub(int a, int b, int k, int l, int r) {
        eval(k, r - l);
        if (r <= a || b <= l) {  // 完全に外側の時
            return ex;
        } else if (a <= l && r <= b) {  // 完全に内側の時
            return dat[k];
        } else {  // 一部区間が被る時
            X vl = query_sub(a, b, k * 2 + 1, l, (l + r) / 2);
            X vr = query_sub(a, b, k * 2 + 2, (l + r) / 2, r);
            return fx(vl, vr);
        }
    }
    //[a,b)
    X query(int a, int b) { return query_sub(a, b, 0, 0, n); }
    X get(int i){
      return query(i,i+1);
    }
};
// RAQ+RSQ
using X = ll; //dat
using M = ll; //lazy
auto fx = [](X x1, X x2) -> X { return x1 + x2; }; //dat×dat 行う二項演算
auto fa = [](X x, M m) -> X { return x + m; }; //datをlazyで更新する演算
auto fm = [](M m1, M m2) -> M { return m1 + m2; }; //lazyの伝達を行う演算
auto fp = [](M m, int n) -> M { return m * n; }; //p(m,n)=m＊m＊...＊m　子のdatをlazyで更新した時のdatの作用
X ex = 0; //x＊ex=x (Xの単位元) dat
M em = 0; //x＊em=x (Mの単位元) lazy
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
    unordered_set<int> st;
    for(int i=0;i<n;i++){
      st.insert(root(i));
    }
    vector<int> res(st.begin(),st.end());
    return res;
  }
};
void dibfs(MaxFlowGraph &G, int s,vector<int> &level) {
    for (int i = 0; i < G.V; ++i) level[i] = -1;
    level[s] = 0;
    queue<int> que;
    que.push(s);
    while (!que.empty()) {
        int v = que.front();
        que.pop();
        for (int i = 0; i < G[v].size(); ++i) {
            edge &e = G[v][i];
            if (level[e.to] < 0 && e.cap > 0) {
                level[e.to] = level[v] + 1;
                que.push(e.to);
            }
        }
    }
}
 
FLOW didfs(MaxFlowGraph &G, int v, int t, FLOW f,vector<int> &level,vector<int> &iter) {
    if (v == t) return f;
    for (int &i = iter[v]; i < G[v].size(); ++i) {
        edge &e = G[v][i], &re = G.redge(e);
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
FLOW Dinic(MaxFlowGraph &G, int s, int t) {
    FLOW res = 0;
    vector<int> level(G.V),iter(G.V);
    while (true) {
        dibfs(G, s,level);
        if (level[t] < 0) return res;
        iter.assign(G.V,0);
        FLOW flow;
        while ((flow = didfs(G, s, t, inf,level,iter)) > 0) {
            res += flow;
        }
    }
}
ll mod(ll a, ll mod) {
    return (a%mod+mod)%mod;
}
ll modpow(ll a,ll n,ll mod){
  ll res=1;
  while (n>0){
    if (n & 1) res = res * a % mod;
    a = a * a % mod;
    n >>= 1;
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
struct Comb{
  unordered_map<int,tuple<ll,vector<ll>,vector<ll>>> mp;
  int n_,m;
  ll p_, pm_;
  vector<ll> ord_, fact_;
  vector<P> pf;
  Comb(int n,int M) : n_(n), ord_(n), fact_(n) { 
    m=M;
    pf=prime_factorize(M); 
  }
  Comb(ll p, ll pm, int n) :
    n_(n), p_(p), pm_(pm), ord_(n), fact_(n) {
    init(p, pm);
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
  ll COM(int n, int k){
    if(n<0 || k<0 || n<k) return 0;
    vector<long long> vb, vm;
    for (auto ps : pf) {
        long long p = ps.first, e = ps.second;
        long long pm = pow(p,e);
        long long b = 1;
        b *= com(n, k ,p) % pm;
        b %= pm;
        vm.push_back(pm);
        vb.push_back(b);
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
  cnt.assign(N,T(0));
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
    long long modpow(long long a, long long n, int mod) {
        long long res = 1;
        while (n > 0) {
            if (n & 1) res = res * a % mod;
            a = a * a % mod;
            n >>= 1;
        }
        return res;
    }
 
    long long modinv(long long a, int mod) {
        long long b = mod, u = 1, v = 0;
        while (b) {
            long long t = a / b;
            a -= t * b, swap(a, b);
            u -= t * v, swap(u, v);
        }
        u %= mod;
        if (u < 0) u += mod;
        return u;
    }
 
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
 
    // number-theoretic transform
    template<class mint> void trans(vector<mint> &v, bool inv = false) {
        if (v.empty()) return;
        int N = (int)v.size();
        int MOD = v[0].getmod();
        int PR = calc_primitive_root(MOD);
        static bool first = true;
        static vector<long long> vbw(30), vibw(30);
        if (first) {
            first = false;
            for (int k = 0; k < 30; ++k) {
                vbw[k] = modpow(PR, (MOD - 1) >> (k + 1), MOD);
                vibw[k] = modinv(vbw[k], MOD);
            }
        }
        for (int i = 0, j = 1; j < N - 1; j++) {
            for (int k = N >> 1; k > (i ^= k); k >>= 1);
            if (i > j) swap(v[i], v[j]);
        }
        for (int k = 0, t = 2; t <= N; ++k, t <<= 1) {
            long long bw = vbw[k];
            if (inv) bw = vibw[k];
            for (int i = 0; i < N; i += t) {
                mint w = 1;
                for (int j = 0; j < t/2; ++j) {
                    int j1 = i + j, j2 = i + j + t/2;
                    mint c1 = v[j1], c2 = v[j2] * w;
                    v[j1] = c1 + c2;
                    v[j2] = c1 - c2;
                    w *= bw;
                }
            }
        }
        if (inv) {
            long long invN = modinv(N, MOD);
            for (int i = 0; i < N; ++i) v[i] = v[i] * invN;
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
    vector<mint> mul(const vector<mint> &A, const vector<mint> &B) {
        if (A.empty() || B.empty()) return {};
        int N = (int)A.size(), M = (int)B.size();
        if (min(N, M) < 30) return naive_mul(A, B);
        int MOD = A[0].getmod();
        int size_fft = get_fft_size(N, M);
        if (MOD == 998244353) {
            vector<mint> a(size_fft), b(size_fft), c(size_fft);
            for (int i = 0; i < N; ++i) a[i] = A[i];
            for (int i = 0; i < M; ++i) b[i] = B[i];
            trans(a), trans(b);
            vector<mint> res(size_fft);
            for (int i = 0; i < size_fft; ++i) res[i] = a[i] * b[i];
            trans(res, true);
            res.resize(N + M - 1);
            return res;
        }
        vector<mint0> a0(size_fft, 0), b0(size_fft, 0), c0(size_fft, 0);
        vector<mint1> a1(size_fft, 0), b1(size_fft, 0), c1(size_fft, 0);
        vector<mint2> a2(size_fft, 0), b2(size_fft, 0), c2(size_fft, 0);
        for (int i = 0; i < N; ++i)
            a0[i] = A[i].val, a1[i] = A[i].val, a2[i] = A[i].val;
        for (int i = 0; i < M; ++i)
            b0[i] = B[i].val, b1[i] = B[i].val, b2[i] = B[i].val;
        trans(a0), trans(a1), trans(a2), trans(b0), trans(b1), trans(b2);
        for (int i = 0; i < size_fft; ++i) {
            c0[i] = a0[i] * b0[i];
            c1[i] = a1[i] * b1[i];
            c2[i] = a2[i] * b2[i];
        }
        trans(c0, true), trans(c1, true), trans(c2, true);
        static const mint mod0 = MOD0, mod01 = mod0 * MOD1;
        vector<mint> res(N + M - 1);
        for (int i = 0; i < N + M - 1; ++i) {
            int y0 = c0[i].val;
            int y1 = (imod0 * (c1[i] - y0)).val;
            int y2 = (imod01 * (c2[i] - y0) - imod1 * y1).val;
            res[i] = mod01 * y2 + mod0 * y1 + y0;
        }
        return res;
    }
 
    // long long
    vector<long long> mul_ll(const vector<long long> &A, const vector<long long> &B) {
        if (A.empty() || B.empty()) return {};
        int N = (int)A.size(), M = (int)B.size();
        if (min(N, M) < 30) return naive_mul(A, B);
        int size_fft = get_fft_size(N, M);
        vector<mint0> a0(size_fft, 0), b0(size_fft, 0), c0(size_fft, 0);
        vector<mint1> a1(size_fft, 0), b1(size_fft, 0), c1(size_fft, 0);
        vector<mint2> a2(size_fft, 0), b2(size_fft, 0), c2(size_fft, 0);
        for (int i = 0; i < N; ++i)
            a0[i] = A[i], a1[i] = A[i], a2[i] = A[i];
        for (int i = 0; i < M; ++i)
            b0[i] = B[i], b1[i] = B[i], b2[i] = B[i];
        trans(a0), trans(a1), trans(a2), trans(b0), trans(b1), trans(b2);
        for (int i = 0; i < size_fft; ++i) {
            c0[i] = a0[i] * b0[i];
            c1[i] = a1[i] * b1[i];
            c2[i] = a2[i] * b2[i];
        }
        trans(c0, true), trans(c1, true), trans(c2, true);
        static const long long mod0 = MOD0, mod01 = mod0 * MOD1;
        vector<long long> res(N + M - 1);
        for (int i = 0; i < N + M - 1; ++i) {
            int y0 = c0[i].val;
            int y1 = (imod0 * (c1[i] - y0)).val;
            int y2 = (imod01 * (c2[i] - y0) - imod1 * y1).val;
            res[i] = mod01 * y2 + mod0 * y1 + y0;
        }
        return res;
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
int main(){
}
