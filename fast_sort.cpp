#include <execution>
#include "bits/stdc++.h"

using namespace std;

using ll = long long;

#define sqr(x) (x)*(x)
#define tip unsigned int

#pragma GCC target("avx2")
#pragma GCC optimization("O3")
#pragma GCC optimization("unroll-loops")

template <typename T>
T Binpow (T a, ll n) {T r = 1;while (n) {if (n & 1) r *= a; a *= a; n >>= 1;} return r;}

template<typename T>
vector<T> create_copy(std::vector<T> const &vec) {
    std::vector<T> v(vec);
    return v;
}

template <typename T>
void P(vector<T> &a) {
    for (auto &x : a) cout << x << ' ';
}

namespace Ivan {
    template<typename T>
    T findMax(vector<T>& arr, T l, T r) {
        T max = arr[l];
        for (auto i = l + 1; i < r; i++)
            if (arr[i] > max) max = arr[i];
        return max;
    }

    template<typename T, typename S>
    void countSort(vector<T>& arr, S base, T exp, T l, T r) {
        vector<T> output(r - l);
        vector<T> count(base, (T)0);

        for (int i = l; i < r; i++) {
            count[(arr[i] / exp) % base]++;
        }
        for (int i = 1; i < base; i++) {
            count[i] += count[i - 1];
        }
        for (int i = r - 1; i >= l; i--) {
            output[count[(arr[i] / exp) % base] - 1] = arr[i];
            count[(arr[i] / exp) % base]--;
        }
        for (int i = l; i < r; i++) {
            arr[i] = output[i];
        }
    }

    template <typename T, typename S>
    void RadixSort(vector<T>& a, S base, T l, T r) {
        auto max = findMax(a, l, r);
        for (T exp = l + 1; max / exp > 0; exp *= base) {
            countSort(a, base, exp, l, r);
        }
    }

    template <typename T, typename S>
    void RadixSort(vector<T>& a, S base) {
        return RadixSort(a, base, (T)0, (T)a.size());
    }

    template <typename T>
    void fast_sort(std::vector<T>& a) {
        RadixSort(a, (uint32_t)Binpow(2, 14));
    }
}

namespace Slava{
    const int osn = 1024;
    const int size_num = 10;
    const int steps = 4;
    unsigned get_last_n_bits( unsigned u, int n )
    {
        return u & ~(~0U << n);
    }
    vector<int>all [osn];
    void radix_sort(std::vector<uint32_t>& tmp, int l, int r, int num_buck){
        if (num_buck==-1){
            return;
        }
        if (r-l<=1){
            return;
        }
        for (int i=l; i<r; i++){
            int diff = get_last_n_bits(tmp[i] >> (size_num*num_buck),size_num);
            all[diff].push_back(tmp[i]);
        }
        int pos = 0;
        int last = l;
        for (int i=l; i < r; i++){
            if (all[pos].empty()){
                while (pos < osn && all[pos].empty()){
                    pos++;
                }
                last = i;
            }
            tmp[i]=all[pos].back();
            all[pos].pop_back();
        }
        int last_pos = l;
        int last_type = get_last_n_bits(tmp[l] >> (size_num*num_buck),size_num);
        for (int i=l; i < r; i++){
            int curtype = get_last_n_bits(tmp[i] >> (size_num*num_buck),size_num);
            if (last_type!=curtype){
                radix_sort(tmp,last_pos,i,num_buck - 1);
                last_pos = i;
                last_type = curtype;
            }
        }
        radix_sort(tmp,last_pos,r,num_buck - 1);
    }

    void fast_sort(std::vector<uint32_t>& a){
        radix_sort(a,0,a.size(),steps - 1);
    }
}



int main() {
    tip max = Binpow(2, 32) - 1;
    tip n = 1e7;
    random_device rd;
    uniform_int_distribution<tip> dist(1, max);
    vector<tip> a;
    for (int i = 0; i < n; i++)
        a.push_back(dist(rd));
    auto t = clock()/(double)1000;
    //Ivan::fast_sort(a);
    Slava::fast_sort(a);
    cout << "time: " << clock()/(double)1000 - t << '\n';
    system("pause");
}
