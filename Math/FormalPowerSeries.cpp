#include<bits/stdc++.h>
#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
#define FOR(i,n) for(int i = 0; i < (n); i++)
#define sz(c) ((int)(c).size())
#define ten(x) ((int)1e##x)
#define all(v) (v).begin(), (v).end()
using namespace std;
using ll=long long;
using P = pair<ll,ll>;
const long double PI=acos(-1);
const ll INF=1e18;
const int inf=1e9;
template<int MOD> struct Fp{
  ll val;
  constexpr Fp(long long v = 0) noexcept : val(v % MOD) {
    if (val < 0) val += MOD;
  }
  static constexpr int getmod() { return MOD; }
  constexpr Fp operator - () const noexcept {
    return val ? MOD - val : 0;
  }
  constexpr Fp operator + (const Fp& r) const noexcept { return Fp(*this) += r; }
  constexpr Fp operator - (const Fp& r) const noexcept { return Fp(*this) -= r; }
  constexpr Fp operator * (const Fp& r) const noexcept { return Fp(*this) *= r; }
  constexpr Fp operator / (const Fp& r) const noexcept { return Fp(*this) /= r; }
  constexpr Fp& operator += (const Fp& r) noexcept {
    val += r.val;
    if (val >= MOD) val -= MOD;
    return *this;
  }
  constexpr Fp& operator -= (const Fp& r) noexcept {
    val -= r.val;
    if (val < 0) val += MOD;
    return *this;
  }
  constexpr Fp& operator *= (const Fp& r) noexcept {
    val = val * r.val % MOD;
    return *this;
  }
  constexpr Fp& operator /= (const Fp& r) noexcept {
    ll a = r.val, b = MOD, u = 1, v = 0;
    while (b) {
      ll t = a / b;
      a -= t * b, swap(a, b);
      u -= t * v, swap(u, v);
    }
    val = val * u % MOD;
    if (val < 0) val += MOD;
    return *this;
  }
  constexpr bool operator == (const Fp& r) const noexcept {
    return this->val == r.val;
  }
  constexpr bool operator != (const Fp& r) const noexcept {
    return this->val != r.val;
  }
  constexpr bool operator < (const Fp& r) const noexcept {
    return this->val < r.val;
  }
  friend constexpr istream& operator >> (istream& is, Fp<MOD>& x) noexcept {
    is >> x.val;
    x.val %= MOD;
    if (x.val < 0) x.val += MOD;
    return is;
  }
  friend constexpr ostream& operator << (ostream& os, const Fp<MOD>& x) noexcept {
    return os << x.val;
  }
  friend constexpr Fp<MOD> modpow(const Fp<MOD>& a, long long n) noexcept {
    Fp<MOD> res=1,r=a;
    while(n){
      if(n&1) res*=r;
      r*=r;
      n>>=1;
    }
    return res;
  }
  friend constexpr Fp<MOD> modinv(const Fp<MOD>& r) noexcept {
        long long a = r.val, b = MOD, u = 1, v = 0;
        while (b) {
            long long t = a / b;
            a -= t * b, swap(a, b);
            u -= t * v, swap(u, v);
        }
        return Fp<MOD>(u);
  }
  explicit operator bool()const{
		return val;
  }
};
ll modpow(ll a,ll n,ll mod){
  ll res=1;
  a%=mod;
  while (n>0){
    if (n & 1) res*=a;
    a *= a;
    a%=mod;
    n >>= 1;
    res%=mod;
  }
  return res;
}
ll modinv(ll a, ll mod) {
    ll b = mod, u = 1, v = 0;
    while (b) {
        ll t = a/b;
        a -= t * b, swap(a, b);
        u -= t * v, swap(u, v);
    }
    u %= mod;
    if (u < 0) u += mod;
    return u;
}
namespace NTT {
    int calc_primitive_root(int mod) {
        if (mod == 2) return 1;
        if (mod == 167772161) return 3;
        if (mod == 469762049) return 3;
        if (mod == 754974721) return 11;
        if (mod == 998244353) return 3;
        int divs[20] = {};
        divs[0] = 2;
        int cnt = 1;
        long long x = (mod - 1) / 2;
        while (x % 2 == 0) x /= 2;
        for (long long i = 3; i * i <= x; i += 2) {
            if (x % i == 0) {
                divs[cnt++] = i;
                while (x % i == 0) x /= i;
            }
        }
        if (x > 1) divs[cnt++] = x;
        for (int g = 2;; g++) {
            bool ok = true;
            for (int i = 0; i < cnt; i++) {
                if (modpow(g, (mod - 1) / divs[i], mod) == 1) {
                    ok = false;
                    break;
                }
            }
            if (ok) return g;
        }
    }
 
    int get_fft_size(int N, int M) {
        int size_a = 1, size_b = 1;
        while (size_a < N) size_a <<= 1;
        while (size_b < M) size_b <<= 1;
        return max(size_a, size_b) << 1;
    }
    constexpr int bsf_constexpr(unsigned int n) {
      int x = 0;
      while (!(n & (1 << x))) x++;
      return x;
    }
    int bsf(unsigned int n) {
      #ifdef _MSC_VER
      unsigned long index;
      _BitScanForward(&index, n);
      return index;
      #else
      return __builtin_ctz(n);
      #endif
    }
    template <class mint>
    struct fft_info{
      static constexpr int rank2 = bsf_constexpr(mint::getmod() - 1);
      std::array<mint, rank2 + 1> root;   // root[i]^(2^i) == 1
      std::array<mint, rank2 + 1> iroot;  // root[i] * iroot[i] == 1
      std::array<mint, std::max(0, rank2 - 2 + 1)> rate2;
      std::array<mint, std::max(0, rank2 - 2 + 1)> irate2;
 
      std::array<mint, std::max(0, rank2 - 3 + 1)> rate3;
      std::array<mint, std::max(0, rank2 - 3 + 1)> irate3;
      int g;
      fft_info(){
        int MOD=mint::getmod();
        g=calc_primitive_root(MOD);
        root[rank2] = modpow(mint(g),(MOD - 1) >> rank2);
        iroot[rank2] = modinv(root[rank2]);
        for (int i = rank2 - 1; i >= 0; i--) {
            root[i] = root[i + 1] * root[i + 1];
            iroot[i] = iroot[i + 1] * iroot[i + 1];
        }
 
        {
            mint prod = 1, iprod = 1;
            for (int i = 0; i <= rank2 - 2; i++) {
                rate2[i] = root[i + 2] * prod;
                irate2[i] = iroot[i + 2] * iprod;
                prod *= iroot[i + 2];
                iprod *= root[i + 2];
            }
        }
        {
            mint prod = 1, iprod = 1;
            for (int i = 0; i <= rank2 - 3; i++) {
                rate3[i] = root[i + 3] * prod;
                irate3[i] = iroot[i + 3] * iprod;
                prod *= iroot[i + 3];
                iprod *= root[i + 3];
            }
        }
      }
    };
    int ceil_pow2(int n) {
      int x = 0;
      while ((1U << x) < (unsigned int)(n)) x++;
      return x;
    }
    // number-theoretic transform
    template <class mint>
    void trans(std::vector<mint>& a) {
      int n = int(a.size());
      int h = ceil_pow2(n);
      int MOD=a[0].getmod();
      static const fft_info<mint> info;
 
      int len = 0;  // a[i, i+(n>>len), i+2*(n>>len), ..] is transformed
      while (len < h) {
        if (h - len == 1) {
            int p = 1 << (h - len - 1);
            mint rot = 1;
            for (int s = 0; s < (1 << len); s++) {
                int offset = s << (h - len);
                for (int i = 0; i < p; i++) {
                    auto l = a[i + offset];
                    auto r = a[i + offset + p] * rot;
                    a[i + offset] = l + r;
                    a[i + offset + p] = l - r;
                }
                if (s + 1 != (1 << len))
                    rot *= info.rate2[bsf(~(unsigned int)(s))];
            }
            len++;
        } else {
            // 4-base
            int p = 1 << (h - len - 2);
            mint rot = 1, imag = info.root[2];
            for (int s = 0; s < (1 << len); s++) {
                mint rot2 = rot * rot;
                mint rot3 = rot2 * rot;
                int offset = s << (h - len);
                for (int i = 0; i < p; i++) {
                    auto mod2 = 1ULL * MOD * MOD;
                    auto a0 = 1ULL * a[i + offset].val;
                    auto a1 = 1ULL * a[i + offset + p].val * rot.val;
                    auto a2 = 1ULL * a[i + offset + 2 * p].val * rot2.val;
                    auto a3 = 1ULL * a[i + offset + 3 * p].val * rot3.val;
                    auto a1na3imag =
                        1ULL * mint(a1 + mod2 - a3).val * imag.val;
                    auto na2 = mod2 - a2;
                    a[i + offset] = a0 + a2 + a1 + a3;
                    a[i + offset + 1 * p] = a0 + a2 + (2 * mod2 - (a1 + a3));
                    a[i + offset + 2 * p] = a0 + na2 + a1na3imag;
                    a[i + offset + 3 * p] = a0 + na2 + (mod2 - a1na3imag);
                }
                if (s + 1 != (1 << len))
                    rot *= info.rate3[bsf(~(unsigned int)(s))];
            }
            len += 2;
        }
      }
    }
    template <class mint>
    void trans_inv(std::vector<mint>& a) {
      int n = int(a.size());
      int h = ceil_pow2(n);
 
      static const fft_info<mint> info;
      int MOD=a[0].getmod();
      int len = h;  // a[i, i+(n>>len), i+2*(n>>len), ..] is transformed
      while (len) {
        if (len == 1) {
            int p = 1 << (h - len);
            mint irot = 1;
            for (int s = 0; s < (1 << (len - 1)); s++) {
                int offset = s << (h - len + 1);
                for (int i = 0; i < p; i++) {
                    auto l = a[i + offset];
                    auto r = a[i + offset + p];
                    a[i + offset] = l + r;
                    a[i + offset + p] =
                        (unsigned long long)(MOD + l.val - r.val) *
                        irot.val;
                    ;
                }
                if (s + 1 != (1 << (len - 1)))
                    irot *= info.irate2[bsf(~(unsigned int)(s))];
            }
            len--;
        } else {
            // 4-base
            int p = 1 << (h - len);
            mint irot = 1, iimag = info.iroot[2];
            for (int s = 0; s < (1 << (len - 2)); s++) {
                mint irot2 = irot * irot;
                mint irot3 = irot2 * irot;
                int offset = s << (h - len + 2);
                for (int i = 0; i < p; i++) {
                    auto a0 = 1ULL * a[i + offset + 0 * p].val;
                    auto a1 = 1ULL * a[i + offset + 1 * p].val;
                    auto a2 = 1ULL * a[i + offset + 2 * p].val;
                    auto a3 = 1ULL * a[i + offset + 3 * p].val;
 
                    auto a2na3iimag =
                        1ULL *
                        mint((MOD + a2 - a3) * iimag.val).val;
 
                    a[i + offset] = a0 + a1 + a2 + a3;
                    a[i + offset + 1 * p] =
                        (a0 + (MOD - a1) + a2na3iimag) * irot.val;
                    a[i + offset + 2 * p] =
                        (a0 + a1 + (MOD - a2) + (MOD - a3)) *
                        irot2.val;
                    a[i + offset + 3 * p] =
                        (a0 + (MOD - a1) + (MOD - a2na3iimag)) *
                        irot3.val;
                }
                if (s + 1 != (1 << (len - 2)))
                    irot *= info.irate3[bsf(~(unsigned int)(s))];
            }
            len -= 2;
        }
      }
    }
    // for garner
    static constexpr int MOD0 = 754974721;
    static constexpr int MOD1 = 167772161;
    static constexpr int MOD2 = 469762049;
    using mint0 = Fp<MOD0>;
    using mint1 = Fp<MOD1>;
    using mint2 = Fp<MOD2>;
    static const mint1 imod0 = 95869806; // modinv(MOD0, MOD1);
    static const mint2 imod1 = 104391568; // modinv(MOD1, MOD2);
    static const mint2 imod01 = 187290749; // imod1 / MOD0;
 
    // small case (T = mint, long long)
    template<class T> vector<T> naive_mul 
    (const vector<T> &A, const vector<T> &B) {
        if (A.empty() || B.empty()) return {};
        int N = (int)A.size(), M = (int)B.size();
        vector<T> res(N + M - 1);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < M; ++j)
                res[i + j] += A[i] * B[j];
        return res;
    }
 
    // mint
    template<class mint>
    vector<mint> mul(vector<mint> A,vector<mint> B) {
        if (A.empty() || B.empty()) return {};
        int n = int(A.size()), m = int(B.size());
        if (min(n, m) < 30) return naive_mul(A, B);
        int MOD = A[0].getmod();
        int z = 1 << ceil_pow2(n + m - 1);
        if (MOD == 998244353) {
          A.resize(z);
          trans(A);
          B.resize(z);
          trans(B);
          for (int i = 0; i < z; i++) {
            A[i] *= B[i];
          }
          trans_inv(A);
          A.resize(n + m - 1);
          mint iz = modinv(mint(z));
          for (int i = 0; i < n + m - 1; i++) A[i] *= iz;
          return A;
        }
        vector<mint0> a0(z, 0), b0(z, 0);
        vector<mint1> a1(z, 0), b1(z, 0);
        vector<mint2> a2(z, 0), b2(z, 0);
        for (int i = 0; i < n; ++i)
            a0[i] = A[i].val, a1[i] = A[i].val, a2[i] = A[i].val;
        for (int i = 0; i < m; ++i)
            b0[i] = B[i].val, b1[i] = B[i].val, b2[i] = B[i].val;
        trans(a0), trans(a1), trans(a2), trans(b0), trans(b1), trans(b2);
        for (int i = 0; i < z; ++i) {
            a0[i] *= b0[i];
            a1[i] *= b1[i];
            a2[i] *= b2[i];
        }
        trans_inv(a0), trans_inv(a1), trans_inv(a2);
        static const mint mod0 = MOD0, mod01 = mod0 * MOD1;
        mint0 i0=modinv(mint0(z));
        mint1 i1=modinv(mint1(z));
        mint2 i2=modinv(mint2(z));
        vector<mint> res(n + m - 1);
        for (int i = 0; i < n + m - 1; ++i) {
            a0[i]*=i0;
            a1[i]*=i1;
            a2[i]*=i2;
            int y0 = a0[i].val;
            int y1 = (imod0 * (a1[i] - y0)).val;
            int y2 = (imod01 * (a2[i] - y0) - imod1 * y1).val;
            res[i] = mod01 * y2 + mod0 * y1 + y0;
        }
        return res;
    }
};
// Formal Power Series
template <typename mint> struct FPS : vector<mint> {
    using vector<mint>::vector;
 /*
    template<class...Args>
    FPS(Args...args) : vector<mint>(args...){}
  */
    // constructor
    FPS(const vector<mint>& r) : vector<mint>(r) {}
 
    // core operator
    inline FPS pre(int siz) const {
        return FPS(begin(*this), begin(*this) + min((int)this->size(), siz));
    }
    inline FPS rev() const {
        FPS res = *this;
        reverse(begin(res), end(res));
        return res;
    }
    inline FPS& normalize() {
        while (!this->empty() && this->back() == 0) this->pop_back();
        return *this;
    }
 
    // basic operator
    inline FPS operator - () const noexcept {
        FPS res = (*this);
        for (int i = 0; i < (int)res.size(); ++i) res[i] = -res[i];
        return res;
    }
    inline FPS operator + (const mint& v) const { return FPS(*this) += v; }
    inline FPS operator + (const FPS& r) const { return FPS(*this) += r; }
    inline FPS operator - (const mint& v) const { return FPS(*this) -= v; }
    inline FPS operator - (const FPS& r) const { return FPS(*this) -= r; }
    inline FPS operator * (const mint& v) const { return FPS(*this) *= v; }
    inline FPS operator * (const FPS& r) const { return FPS(*this) *= r; }
    inline FPS operator / (const mint& v) const { return FPS(*this) /= v; }
    inline FPS operator << (int x) const { return FPS(*this) <<= x; }
    inline FPS operator >> (int x) const { return FPS(*this) >>= x; }
    inline FPS& operator += (const mint& v) {
        if (this->empty()) this->resize(1);
        (*this)[0] += v;
        return *this;
    }
    inline FPS& operator += (const FPS& r) {
        if (r.size() > this->size()) this->resize(r.size());
        for (int i = 0; i < (int)r.size(); ++i) (*this)[i] += r[i];
        return this->normalize();
    }
    inline FPS& operator -= (const mint& v) {
        if (this->empty()) this->resize(1);
        (*this)[0] -= v;
        return *this;
    }
    inline FPS& operator -= (const FPS& r) {
        if (r.size() > this->size()) this->resize(r.size());
        for (int i = 0; i < (int)r.size(); ++i) (*this)[i] -= r[i];
        return this->normalize();
    }
    inline FPS& operator *= (const mint& v) {
        for (int i = 0; i < (int)this->size(); ++i) (*this)[i] *= v;
        return *this;
    }
    inline FPS& operator *= (const FPS& r) {
        return *this = NTT::mul((*this), r);
    }
    inline FPS& operator /= (const mint& v) {
        assert(v != 0);
        mint iv = modinv(v);
        for (int i = 0; i < (int)this->size(); ++i) (*this)[i] *= iv;
        return *this;
    }
    inline FPS& operator <<= (int x) {
        FPS res(x, 0);
        res.insert(res.end(), begin(*this), end(*this));
        return *this = res;
    }
    inline FPS& operator >>= (int x) {
        FPS res;
        res.insert(res.end(), begin(*this) + x, end(*this));
        return *this = res;
    }
    inline mint eval(const mint& v){
        mint res = 0;
        for (int i = (int)this->size()-1; i >= 0; --i) {
            res *= v;
            res += (*this)[i];
        }
        return res;
    }
    inline friend FPS gcd(const FPS& f, const FPS& g) {
        if (g.empty()) return f;
        return gcd(g, f % g);
    }

    // advanced operation
    // df/dx
    inline friend FPS diff(const FPS& f) {
        int n = (int)f.size();
        FPS res(n-1);
        for (int i = 1; i < n; ++i) res[i-1] = f[i] * i;
        return res;
    }

    // \int f dx
    inline friend FPS integral(const FPS& f) {
        int n = (int)f.size();
        FPS res(n+1, 0);
        for (int i = 0; i < n; ++i) res[i+1] = f[i] / (i+1);
        return res;
    }

    // inv(f), f[0] must not be 0
    inline friend FPS inv(const FPS& f, int deg) {
        assert(f[0] != 0);
        if (deg < 0) deg = (int)f.size();
        FPS res({mint(1) / f[0]});
        for (int i = 1; i < deg; i <<= 1) {
            res = (res + res - res * res * f.pre(i << 1)).pre(i << 1);
        }
        res.resize(deg);
        return res;
    }
    inline friend FPS inv(const FPS& f) {
        return inv(f, f.size());
    }

    // division, r must be normalized (r.back() must not be 0)
    inline FPS& operator /= (const FPS& r) {
        const int n=(*this).size(),m=r.size();
        if(n<m){
            (*this).clear();
            return *this;
        }
        assert(r.back() != 0);
        this->normalize();
        if (this->size() < r.size()) {
            this->clear();
            return *this;
        }
        int need = (int)this->size() - (int)r.size() + 1;
        *this = ((*this).rev().pre(need) * inv(r.rev(), need)).pre(need).rev();
        return *this;
    }
    inline FPS& operator %= (const FPS &r) {
        const int n=(*this).size(),m=r.size();
        if(n<m) return (*this);
        assert(r.back() != 0);
        this->normalize();
        FPS q = (*this) / r;
        return *this -= q * r;
    }
    inline FPS operator / (const FPS& r) const { return FPS(*this) /= r; }
    inline FPS operator % (const FPS& r) const { return FPS(*this) %= r; }

    // log(f) = \int f'/f dx, f[0] must be 1
    inline friend FPS log(const FPS& f, int deg) {
        assert(f[0] == 1);
        FPS res = integral((diff(f) * inv(f, deg)).pre(deg-1));
        return res;
    }
    inline friend FPS log(const FPS& f) {
        return log(f, f.size());
    }

    // exp(f), f[0] must be 0
    inline friend FPS exp(const FPS& f, int deg) {
        assert(f[0] == 0);
        FPS res(1, 1);
        for (int i = 1; i < deg; i <<= 1) {
            res = res * (f.pre(i<<1) - log(res, i<<1) + 1).pre(i<<1);
        }
        res.resize(deg);
        return res;
    }
    inline friend FPS exp(const FPS& f) {
        return exp(f, f.size());
    }

    // pow(f) = exp(e * log f)
    inline friend FPS pow(const FPS& f, long long e, int deg) {
        long long i = 0;
        while (i < (int)f.size() && f[i] == 0) ++i;
        if (i == (int)f.size()) return FPS(deg, 0);
        if (i * e >= deg) return FPS(deg, 0);
        mint k = f[i];
        FPS res = exp(log((f >> i) / k, deg) * e, deg) * modpow(k, e) << (e * i);
        res.resize(deg);
        return res;
    }
    inline friend FPS pow(const FPS& f, long long e) {
        return pow(f, e, f.size());
    }

    // sqrt(f), f[0] must be 1
    inline friend FPS sqrt_base(const FPS& f, int deg) {
        assert(f[0] == 1);
        mint inv2 = mint(1) / 2;
        FPS res(1, 1);
        for (int i = 1; i < deg; i <<= 1) {
            res = (res + f.pre(i << 1) * inv(res, i << 1)).pre(i << 1);
            for (mint& x : res) x *= inv2;
        }
        res.resize(deg);
        return res;
    }
    inline friend FPS sqrt_base(const FPS& f) {
        return sqrt_base(f, f.size());
    }
    FPS taylor_shift(mint c) const {
      int n = (int) this->size();
      vector<mint> fact(n), rfact(n);
      fact[0] = rfact[0] = mint(1);
      for(int i = 1; i < n; i++) fact[i] = fact[i - 1] * mint(i);
      rfact[n - 1] = mint(1) / fact[n - 1];
      for(int i = n - 1; i > 1; i--) rfact[i - 1] = rfact[i] * mint(i);
      FPS p(*this);
      for(int i = 0; i < n; i++) p[i] *= fact[i];
      p = p.rev();
      FPS bs(n, mint(1));
      for(int i = 1; i < n; i++) bs[i] = bs[i - 1] * c * rfact[i] * fact[i - 1];
      p = (p * bs).pre(n);
      p = p.rev();
      for(int i = 0; i < n; i++) p[i] *= rfact[i];
      return p;
    }
};
template <typename mint> FPS<mint> product(vector<FPS<mint>> a){
  int siz=1;
  while(siz<int(a.size())){
    siz<<=1;
  }
  vector<FPS<mint>> res(siz*2-1,{1});
  for(int i=0;i<int(a.size());++i){
    res[i+siz-1]=a[i];
  }
  for(int i=siz-2;i>=0;--i){
    res[i]=res[2*i+1]*res[2*i+2];
  }
  return res[0];
}
template<typename mint> FPS<mint> inv_sum(vector<FPS<mint>> f){
  int siz=1;
  while(siz<int(f.size())){
    siz<<=1;
  }
  vector<FPS<mint>> mol(siz*2-1),dem(siz*2-1,{1});
  for(size_t i=0;i<f.size();++i){
    mol[i+siz-1]={1};
    dem[i+siz-1]=f[i];
  }
  for(int i=siz-2;i>=0;--i){
    dem[i]=dem[2*i+1]*dem[2*i+2];
    mol[i]=mol[2*i+1]*dem[2*i+2]+mol[2*i+2]*dem[2*i+1];
  }
  mol[0]*=inv(dem[0]);
  return mol[0];
}
using mint=Fp<998244353>;
int main(){
}
