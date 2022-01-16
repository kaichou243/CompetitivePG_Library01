#include<bits/stdc++.h>
using namespace std;
using ll=long long;
using P=pair<ll,ll>;
struct UnionFind{
  vector<int> data;
  UnionFind(int N) : data(N,-1){}
  inline int root(int x){
    return data[x]<0?x:data[x]=root(data[x]);
  }
  inline int unite(int x, int y) {
    if ((x=root(x))==(y=root(y))) return false;
    if (data[x]>data[y]) swap(x, y);
    data[x]+=data[y];
    data[y]=x;
    return true;
  }
  inline bool same(int x, int y){
    return root(x)==root(y);
  }
  inline int size(int x){
    return -data[root(x)];
  }
};
//example
//G=(N,M) 無向グラフ
//0<=a_i,b_i,S,T<N
int main(){
  int N,M,S,T;
  cin>>N>>M>>S>>T;
  UnionFind uf(N);
  for(int i=0;i<N;i++){
    int a,b;
    cin>>a>>b;
    uf.unite(a,b);
  }
  //頂点Sと頂点Tの連結判定 (連結:1,not 連結:0)
  cout<<same(S,T)<<endl;
  //頂点Sの連結成分の大きさ
  cout<<size(S)<<endl;
}
