// AtCoder Typical Contest 001 C 高速フーリエ変換(Fast Fourier Transform)
#include<bits/stdc++.h>
using namespace std;
using ll=long long;
const double PI=acos(-1);
namespace FFT{
  using Real = double;
  struct Comp{
    Real re, im;
    Comp(){};
    Comp(Real re, Real im) : re(re), im(im) {}
    Real slen() const { return re*re+im*im; }
    Real real() { return re; }
    inline Comp operator+(const Comp &c) { return Comp(re+c.re, im+c.im); }
    inline Comp operator-(const Comp &c) { return Comp(re-c.re, im-c.im); }
    inline Comp operator*(const Comp &c) { return Comp(re*c.re-im*c.im, re*c.im+im*c.re); }
    inline Comp operator/(const Real &c) { return Comp(re/c, im/c); } 
    inline Comp conj() const { return Comp(re, -im); }
  };
  void fft(vector<Comp>& a, bool inverse){
    int n = a.size();
    static bool prepared = false;
    static vector<Comp> z(30);
    if(!prepared){
      prepared = true;
      for(int i=0; i<30; i++){
        Real ang = 2*PI/(1 << (i+1));
        z[i] = Comp(cos(ang), sin(ang));
      }
    }
    for (size_t i = 0, j = 1; j < n; ++j) {
      for (size_t k = n >> 1; k > (i ^= k); k >>= 1);
      if (i > j) swap(a[i], a[j]);
    }
    for (int i = 0, t = 1; t < n; ++i, t <<= 1) {
      Comp bw = z[i];
      if(inverse) bw = bw.conj();
      for (int i = 0; i < n; i += t * 2) {
        Comp w(1,0);
        for (int j = 0; j < t; ++j) {
          int l = i + j, r = i + j + t;
          Comp c = a[l], d = a[r] * w;
          a[l] = c + d;
          a[r] = c - d;
          w = w * bw;
        }
      }
    }
    if(inverse)
      for(int i=0; i<n; i++) a[i] = a[i]/(Real)n;
  }
  template<typename T>
  vector<T> mul(vector<T> &a, vector<T> &b, bool issquare){
    int deg = a.size() + b.size();
    int n = 1;
    while(n < deg) n <<= 1;
    vector<Comp> fa,fb;
    fa.resize(n); fb.resize(n);
    for(int i=0;i<n;i++){
      if(i < a.size()) fa[i] = Comp(a[i], 0);
      else fa[i] = Comp(0,0);
      if(i < b.size()) fb[i] = Comp(b[i], 0);
      else fb[i] = Comp(0,0);
    }
    for(int i=0;i<b.size();i++) fb[i] = Comp(b[i], 0);
    fft(fa, false); 
    if(issquare) fb = fa;
    else fft(fb, false);
    for(int i=0;i<n;i++) fa[i] = fa[i] * fb[i];
    fft(fa, true);
    vector<int> res(n);
    for(int i=0;i<n;i++) res[i] = round(fa[i].re);
    return res;
  }
};
int main(){
  int N;
  cin>>N;
  vector<int> A(N),B(N);
  for(int i=0;i<N;i++) cin>>A[i]>>B[i];
  auto conv=FFT::mul(A,B,false);
  cout<<0<<endl;
  for(int i=0;i<2*N-1;i++) cout<<conv[i]<<endl;
}
