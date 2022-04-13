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
