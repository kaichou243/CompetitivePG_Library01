/*
注意
このセグ木は使わないことを強くお勧めします。
*/
#include<bits/stdc++.h>
using namespace std;
using ll=long long;
using P=pair<ll,ll>;
struct SegmentTree{
  ll g(int a,int b){
    return a+b; //g(a,b) : aとbで二項演算を行った結果を返す関数
  }
  int n;
  vector<ll> dat;
  SegmentTree(int N){
    n=1;
    while(n<N)n*=2;
    dat.assign(2*n-1,0);
  }
  //加算クエリ
  void add(int k,ll x){
    k+=n-1;
    dat[k]+=x;
    while(k){
      k=(k-1)/2;
      dat[k]=g(dat[k*2+1],dat[k*2+2]);
    }
  }
  //更新クエリ
  void update(int k,ll x){
    k+=n-1;
    dat[k]=x;
    while(k){
      k=(k-1)/2;
      dat[k]=g(dat[k*2+1],dat[k*2+2]);
    }
  }
  ll querry(int l,int r){
    return querry_sub(l,r,0,0,n); //区間[l,r)に対する二項演算の総積を返す
  }
  ll querry_sub(int l,int r,int k=0,int a=0,int b=-1){
    if(b<0)b=n;
    if(r<=a || b<=l)return 0; //使う二項演算に合わせて単位元を変える。要注意!
    if(l<=a && b<=r)return dat[k];
    return g(querry_sub(l,r,k*2+1,a,(a+b)/2),querry_sub(l,r,k*2+2,(a+b)/2,b));
  }
};
//example
//0-indexed
int main(){
  int N,Q;
  cin>>N>>Q;
  SegmentTree st(N);
  for(int i=0;i<N;i++){
    int a;
    cin>>a;
    st.add(i,a);
  }
  while(Q--){
    int c;
    cin>>c;
    if(c==1){
      int k,x;
      cin>>k>>x;
      //加算クエリ a_kにxを加算
      st.add(k,x);
    }else if(c==2){
      int l,r;
      cin>>l>>r;
      //区間[l,r)のクエリを返す この場合は[l,r)の総和
      cout<<st.querry(l,r)<<endl;
    }
  }
}
