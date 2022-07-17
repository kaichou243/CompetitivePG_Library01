/**
 * 重み付きの辺を処理できるUnionFindです。
 * 重みT及び距離関数f: (T, T) -> Tが与えられるとします。
 * この時、2点x, yを定めると、xからyへ向かうパスが存在するか否か、またパスが存在する場合はパス間の距離を求めることができます。
 * fはモノイドで良いです。
 * @tparam T 重みの型
 * @tparam Fn 距離関数の型
 */
template<typename T, typename Fn = T(*)(T, T)>
struct WeightedUnionFind {
  std::vector<int> parent;
  std::vector<T> parentWeight, childWeight;
  Fn f;
  /**
   * @param order 頂点数
   * @param f 距離の和、dist(x, z) = f(dist(x, y), dist(y, z))となる
   */
  WeightedUnionFind(int order, Fn f) : parent(order, -1), parentWeight(order), childWeight(order), f(f) {}
  int root(int k) {
    int lastParent = parent[k];
    if (lastParent < 0) return k;
    int par = root(lastParent);
    parentWeight[k] = f(parentWeight[k], parentWeight[lastParent]);
    childWeight[k] = f(childWeight[lastParent], childWeight[k]);
    return parent[k] = par;
  }
  
  /**
   * 頂点x, yをdist(x, y)=weight, dist(y, x)=revWeightとして定めます。
   * なお、既にx->y間に辺が存在する場合は何もしません。
   * @param x 始点
   * @param y 終点
   * @param weight xからyへの距離
   * @param revWeight yからxへの距離
   * @return 元々xとyが連結ならばtrue、そうでないならばfalse
   */
  bool unite(int x, int y, T weight, T revWeight) {
    root(x);
    root(y);
    weight = f(f(childWeight[x], weight), parentWeight[y]);
    revWeight = f(f(childWeight[y], revWeight), parentWeight[x]);
    x = root(x), y = root(y);
    if(x == y) return false;
    parent[y] += parent[x];
    parent[x] = y;
    parentWeight[x] = weight;
    childWeight[x] = revWeight;
    return true;
  }
  /**
   * xとyが連結であるか判定します。
   * @param x 頂点
   * @param y 頂点
   * @return xとyが連結ならばtrue
   */
  bool same(int x,int y){
    return root(x)==root(y);
  }
  /**
   * 頂点xを含む連結成分の大きさを返します。
   * @param x 頂点
   * @return 頂点xを含む連結成分の大きさ
   */
  int size(int x){
    return -parent[root(x)];
  }
  /**
   * xからyへ向かう距離を返します。
   * @param x 始点
   * @param y 終点
   * @return xからyへ向かう距離
   */
  T dist(int x, int y) {
    root(x);
    root(y);
    return f(parentWeight[x], childWeight[y]);
  }
};
