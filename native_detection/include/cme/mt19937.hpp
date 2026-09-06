// MATLAB-compatible Mersenne Twister ("twister" generator of rng()).
//
// rng(seed) with the twister generator seeds mt19937ar with init_genrand(seed)
// (seed 0 is mapped to the MT reference default 5489).  rand() returns
// genrand_res53 doubles, randi(n) == ceil(n*rand()).  Both facts were verified
// against MATLAB R2025a (see PORTING_NOTES.md §3.2).
//
// Reference algorithm: Matsumoto & Nishimura, mt19937ar.c (BSD licence).
#pragma once
#include <cstdint>
#include <cmath>
#include <vector>

namespace cme {

class MatlabTwister {
public:
    explicit MatlabTwister(uint32_t seed = 0) { seed_matlab(seed); }

    // rng(seed): MATLAB maps seed 0 to 5489.
    void seed_matlab(uint32_t seed) {
        init_genrand(seed == 0 ? 5489u : seed);
    }

    void init_genrand(uint32_t s) {
        mt_[0] = s;
        for (mti_ = 1; mti_ < N; ++mti_) {
            mt_[mti_] = (1812433253u * (mt_[mti_ - 1] ^ (mt_[mti_ - 1] >> 30)) + static_cast<uint32_t>(mti_));
        }
    }

    uint32_t genrand_int32() {
        static const uint32_t mag01[2] = {0x0u, MATRIX_A};
        uint32_t y;
        if (mti_ >= N) {
            int kk;
            for (kk = 0; kk < N - M; ++kk) {
                y = (mt_[kk] & UPPER_MASK) | (mt_[kk + 1] & LOWER_MASK);
                mt_[kk] = mt_[kk + M] ^ (y >> 1) ^ mag01[y & 0x1u];
            }
            for (; kk < N - 1; ++kk) {
                y = (mt_[kk] & UPPER_MASK) | (mt_[kk + 1] & LOWER_MASK);
                mt_[kk] = mt_[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1u];
            }
            y = (mt_[N - 1] & UPPER_MASK) | (mt_[0] & LOWER_MASK);
            mt_[N - 1] = mt_[M - 1] ^ (y >> 1) ^ mag01[y & 0x1u];
            mti_ = 0;
        }
        y = mt_[mti_++];
        y ^= (y >> 11);
        y ^= (y << 7) & 0x9d2c5680u;
        y ^= (y << 15) & 0xefc60000u;
        y ^= (y >> 18);
        return y;
    }

    // MATLAB rand: genrand_res53, uniform on [0,1) with 53-bit resolution.
    double rand() {
        uint32_t a = genrand_int32() >> 5;
        uint32_t b = genrand_int32() >> 6;
        return (a * 67108864.0 + b) * (1.0 / 9007199254740992.0);
    }

    // MATLAB randi(n): 1-based integer in [1, n].  Equals ceil(n*rand()).
    long randi(long n) {
        double u = rand();
        long v = static_cast<long>(std::ceil(n * u));
        if (v < 1) v = 1;       // only if u == 0 exactly
        if (v > n) v = n;
        return v;
    }

    // randperm(n,k): partial Fisher-Yates.  For i = 1..k draw j = randi(n-i+1),
    // emit a[j], then move the last remaining element into slot j.
    // Verified against MATLAB for randperm(100,2) = [42 72] (rng(1)); the
    // small-n branch of MATLAB (e.g. randperm(10,3)) is NOT reproduced - this
    // function is only reachable through an unreachable degenerate branch of
    // gmcluster's k-means++ initialisation (see PORTING_NOTES.md).
    std::vector<long> randperm_k(long n, long k) {
        std::vector<long> a(static_cast<size_t>(n));
        for (long i = 0; i < n; ++i) a[static_cast<size_t>(i)] = i + 1;
        std::vector<long> out;
        out.reserve(static_cast<size_t>(k));
        for (long i = 0; i < k; ++i) {
            const long m = n - i;             // remaining count
            const long j = randi(m) - 1;      // 0-based slot in [0, m)
            out.push_back(a[static_cast<size_t>(j)]);
            a[static_cast<size_t>(j)] = a[static_cast<size_t>(m - 1)];
        }
        return out;
    }

private:
    static constexpr int N = 624;
    static constexpr int M = 397;
    static constexpr uint32_t MATRIX_A = 0x9908b0dfu;
    static constexpr uint32_t UPPER_MASK = 0x80000000u;
    static constexpr uint32_t LOWER_MASK = 0x7fffffffu;
    uint32_t mt_[N];
    int mti_ = N + 1;
};

} // namespace cme
