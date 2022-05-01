template <typename mint> struct PartitionNumber{
    using poly=FPS<mint>;
    poly f;
    PartitionNumber(int maxN,vector<int> st) : f(maxN+1){
        vector<mint> invs(maxN+1);
        for(int i=1;i<=maxN;i++) invs[i]=modinv(mint(i));
        for(auto& e:st){
            for(int j=1;j<=maxN/e;j++){
                f[e*j]+=invs[j];
            }
        }
        f=exp(f);
    }
    PartitionNumber(int maxN) : f(maxN+1){
        vector<mint> invs(maxN+1);
        for(int i=1;i<=maxN;i++) invs[i]=modinv(mint(i));
        for(int i=1;i<=maxN;i++){
            for(int j=1;j<=maxN/i;j++){
                f[i*j]+=invs[j];
            }
        }
        f=exp(f);
    }
    mint operator()(int n){
        return f[n];
    }
};
