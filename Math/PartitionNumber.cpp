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
        f[0]=mint(1);
        for(int k=1;k<=maxN;k++){
            if(1ll*k*(3*k+1)/2<=maxN) f[k*(3*k+1)/2]+=(k % 2 ? mint(-1) : mint(1));
            if(1ll*k*(3*k-1)/2<=maxN) f[k*(3*k-1)/2]+=(k % 2 ? mint(-1) : mint(1));
        }
        f=inv(f);
    }
    mint operator()(int n){
        return f[n];
    }
};
