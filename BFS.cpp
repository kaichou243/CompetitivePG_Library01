#include<bits/stdc++.h>
using namespace std;
using ll=long long;
using P=pair<ll,ll>;
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
//example
//G=(N,M) 無向グラフ
int main(){
  int N,M,S,T;
  cin>>N>>M>>S>>T;
  vector<vector<int>> G(N);
  for(int i=0;i<M;i++){
    int a,b;
    cin>>a>>b;
    G[a].push_back(b);
    G[b].push_back(a);
  }
  vector<int> dist;
  bfs(G,S,dist);
  cout<<dist[T]<<endl; //Path(S,T)
}
