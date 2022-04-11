template <typename T>
vector<int> compress(vector<T> &x){
  vector<T> vals=x;
  vector<int> res;
  sort(vals.begin(),vals.end());
  auto it=unique(vals.begin(),vals.end());
  vals.erase(it,vals.end());
  for(auto xx : x){
    int ret=lower_bound(vals.begin(),vals.end(),xx) - vals.begin();
    res.push_back(ret);
  }
  return res;
}
