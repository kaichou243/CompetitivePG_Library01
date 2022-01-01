#include<bits/stdc++.h>
using namespace std;
using ll=long long;
using P=pair<ll,ll>;
template<typename T> struct SegmentTree{
  using F=function<T(T,T)>;
  F fT;
  const T et;
  int n;
  vector<T> dat;
  SegmentTree(int N,F fT_,T et_) : fT(fT_),et(et_){
    n=1;
    while(n<N)n*=2;
    dat.assign(2*n-1,et_);
  }
  void add(int k,T x){
    k+=n-1;
    dat[k]+=x;
    while(k){
      k=(k-1)/2;
      dat[k]=fT(dat[k*2+1],dat[k*2+2]);
    }
  }
  void apply(int k,T x){
    k+=n-1;
    dat[k]=fT(dat[k],x);
    while(k){
      k=(k-1)/2;
      dat[k]=fT(dat[k*2+1],dat[k*2+2]);
    }
  }
  void update(int k,T x){
    k+=n-1;
    dat[k]=x;
    while(k){
      k=(k-1)/2;
      dat[k]=fT(dat[k*2+1],dat[k*2+2]);
    }
  }
  T query(int l,int r){
    return query_sub(l,r,0,0,n);
  }
  T query_sub(int l,int r,int k=0,int a=0,int b=-1){
    if(b<0)b=n;
    if(r<=a || b<=l)return et;
    if(l<=a && b<=r)return dat[k];
    return fT(query_sub(l,r,k*2+1,a,(a+b)/2),query_sub(l,r,k*2+2,(a+b)/2,b));
  }
};
//example
//0-indexed
//PAQ+RSQ
auto f = [](ll x1,ll x2) -> ll { return x1+x2; };
const ll e=0;
int main(){
  int N,Q;
  cin>>N>>Q;
  SegmentTree<ll> st(N,f,e);
  for(int i=0;i<N;i++){
    int a;
    cin>>a;
    st.add(i,a);
  }
  while(Q--){
    int c;
    cin>>c;
    if(c==1){
      int k;
      ll x;
      cin>>k>>x;
      //add a_k x
      st.add(k,x);
    }else if(c==2){
      int l,r;
      cin>>l>>r;
      //return product of [l,r)
      cout<<st.query(l,r)<<endl;
    }
  }
}
