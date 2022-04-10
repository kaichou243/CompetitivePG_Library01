#include<bits/stdc++.h>
using namespace std;
using ll=long long;
using P=pair<ll,ll>;
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
int main(){
}
