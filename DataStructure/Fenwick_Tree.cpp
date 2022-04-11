using ll=long long;
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
