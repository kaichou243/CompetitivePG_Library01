#include<bits/stdc++.h>
using namespace std;
using ll=long long;
using P=pair<ll,ll>;
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
int main(){
}
