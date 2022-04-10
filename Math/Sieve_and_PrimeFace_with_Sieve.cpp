//[1,N]
vector<int> Eratosthenes(int N) {
  vector<int> isprime(N+1,1);
  isprime[1]=0;
  for (int p=2;p<=N;++p){
    if (!isprime[p]) continue;
    for (int q=p*2;q<=N;q+=p) {
      isprime[q]=0;
    }
  }
  return isprime;
}
//log X 素因数分解
struct PrimeFact{
  vector<int> era;
  PrimeFact(int maxA) : era(maxA+1,-1){
    for(int p=2;p<=maxA;++p){
      if(era[p]!=-1) continue;
      for(int q=p;q<=maxA;q+=p){
        era[q]=p;
      }
    }
  }
  vector<P> prime_fact(int x){
    vector<P> res;
    while(1<x){
      int p=era[x],ex=0;
      while(x%p==0){
        x/=p;
        ex++;
      }
      res.push_back({p,ex});
    }
    reverse(res.begin(),res.end());
    return res;
  }
};
