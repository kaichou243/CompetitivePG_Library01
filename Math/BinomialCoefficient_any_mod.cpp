#include<bits/stdc++.h>
using namespace std;
using ll=long long;
using P=pair<ll,ll>;
// mod function
ll mod(ll a, ll mod) {
    return (a%mod+mod)%mod;
}
ll modpow(ll a, ll n, ll mod) {
    ll res = 1;
    while (n > 0) {
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
  void COMinit(int n){
    for(auto p : pf){
      ll ps=p.first,pfs=pow(p.first,p.second);
      init(n);
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
//example
int main(){
  Comb C(200000,998244353);
  C.COMinit(200000);
  cout<<C.COM(5,2);
}
