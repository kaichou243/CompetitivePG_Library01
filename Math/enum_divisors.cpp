vector<ll> enum_divisors(ll N){
  vector<ll> res;
  for (ll i = 1; i * i <= N; ++i){
    if (N % i == 0){
       res.push_back(i);
       // 重複しないならば i の相方である N/i も push
       if (N/i != i) res.push_back(N/i);
    }
  }
  // 小さい順に並び替える
  sort(res.begin(), res.end());
  return res;
}
