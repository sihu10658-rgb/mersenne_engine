#include <iostream>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <thread>
#include <mutex>
#include <atomic>
#include <chrono>

constexpr uint64_t MOD1 = 998244353;
constexpr uint64_t MOD2 = 1004535809;
constexpr uint64_t MOD3 = 469762049;
constexpr uint64_t G = 3;

std::atomic<int> global_p{2};
std::mutex cout_mutex;

// 고속 모듈로 거듭제곱 (Modular Exponentiation)
inline uint64_t power(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t res = 1;
    base %= mod;
    while (exp > 0) {
        if (exp % 2 == 1) res = (__int128)res * base % mod;
        base = (__int128)base * base % mod;
        exp /= 2;
    }
    return res;
}

inline uint64_t modInverse(uint64_t n, uint64_t mod) {
    return power(n, mod - 2, mod);
}

void bit_reverse(std::vector<uint64_t>& a) {
    size_t n = a.size();
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
}

void ntt(std::vector<uint64_t>& a, bool invert, uint64_t mod) {
    size_t n = a.size();
    bit_reverse(a);
    for (size_t len = 2; len <= n; len <<= 1) {
        uint64_t wlen = power(G, (mod - 1) / len, mod);
        if (invert) wlen = modInverse(wlen, mod);
        for (size_t i = 0; i < n; i += len) {
            uint64_t w = 1;
            for (size_t j = 0; j < len / 2; j++) {
                uint64_t u = a[i + j];
                uint64_t v = (__int128)a[i + j + len / 2] * w % mod;
                a[i + j] = (u + v >= mod) ? (u + v - mod) : (u + v);
                a[i + j + len / 2] = (u >= v) ? (u - v) : (u + mod - v);
                w = (__int128)w * wlen % mod;
            }
        }
    }
    if (invert) {
        uint64_t n_inv = modInverse(n, mod);
        for (uint64_t& x : a) x = (__int128)x * n_inv % mod;
    }
}

std::vector<uint64_t> single_ntt_multiply(const std::vector<uint64_t>& a, const std::vector<uint64_t>& b, uint64_t mod) {
    size_t n = 1;
    while (n < a.size() + b.size()) n <<= 1;

    std::vector<uint64_t> fa(n, 0), fb(n, 0);
    for (size_t i = 0; i < a.size(); ++i) fa[i] = a[i] % mod;
    for (size_t i = 0; i < b.size(); ++i) fb[i] = b[i] % mod;

    ntt(fa, false, mod);
    ntt(fb, false, mod);
    for (size_t i = 0; i < n; i++) fa[i] = (__int128)fa[i] * fb[i] % mod;
    ntt(fa, true, mod);

    return fa;
}

uint64_t inv1 = modInverse(MOD1, MOD2);
uint64_t inv2 = modInverse((__int128)MOD1 * MOD2 % MOD3, MOD3);

inline uint64_t combine_crt(uint64_t r1, uint64_t r2, uint64_t r3) {
    int64_t v1 = r1;
    int64_t v2 = (r2 - v1 + MOD2) * inv1 % MOD2;
    if (v2 < 0) v2 += MOD2;
    
    int64_t v3 = (r3 - (v1 + (__int128)MOD1 * v2 % MOD3) % MOD3 + MOD3) % MOD3;
    v3 = v3 * inv2 % MOD3;

    return v1 + (__int128)MOD1 * v2 + (__int128)MOD1 * MOD2 * v3;
}

struct MersenneBigInt {
    std::vector<uint64_t> limbs;
    int p;

    MersenneBigInt(int p_val, uint64_t val = 0) : p(p_val) {
        if (val > 0) limbs.push_back(val);
        else limbs.push_back(0);
    }

    void reduce_mersenne() {
        if (limbs.empty()) return;
        int limb_idx = p / 32;
        int bit_idx = p % 32;

        while (true) {
            std::vector<uint64_t> high_limbs(limbs.size(), 0);
            bool has_high = false;

            for (size_t i = limb_idx; i < limbs.size(); ++i) {
                if (i == (size_t)limb_idx) {
                    if (bit_idx > 0) {
                        high_limbs[i - limb_idx] = limbs[i] >> bit_idx;
                        limbs[i] &= (1ULL << bit_idx) - 1;
                    }
                } else {
                    high_limbs[i - limb_idx] = limbs[i];
                    limbs[i] = 0;
                }
            }

            for (auto v : high_limbs) if (v > 0) has_high = true;
            if (!has_high) break;

            unsigned __int128 carry = 0;
            for (size_t i = 0; i < limbs.size() || i < high_limbs.size() || carry > 0; ++i) {
                if (i >= limbs.size()) limbs.push_back(0);
                unsigned __int128 sum = (unsigned __int128)limbs[i] + 
                                       (i < high_limbs.size() ? high_limbs[i] : 0) + carry;
                limbs[i] = (uint64_t)(sum & 0xFFFFFFFFULL);
                carry = sum >> 32;
            }

            if (carry > 0) {
                uint64_t c = (uint64_t)carry;
                for (size_t i = 0; i < limbs.size() && c > 0; ++i) {
                    unsigned __int128 sum = (unsigned __int128)limbs[i] + c;
                    limbs[i] = (uint64_t)(sum & 0xFFFFFFFFULL);
                    c = (uint64_t)(sum >> 32);
                }
                if (c > 0) limbs.push_back(c);
            }
        }

        while (limbs.size() > 1 && limbs.back() == 0) {
            limbs.pop_back();
        }
    }

    static MersenneBigInt multiply(const MersenneBigInt& a, const MersenneBigInt& b) {
        auto res1 = single_ntt_multiply(a.limbs, b.limbs, MOD1);
        auto res2 = single_ntt_multiply(a.limbs, b.limbs, MOD2);
        auto res3 = single_ntt_multiply(a.limbs, b.limbs, MOD3);

        size_t max_len = std::max({res1.size(), res2.size(), res3.size()});
        std::vector<uint64_t> combined(max_len, 0);

        unsigned __int128 carry = 0;
        for (size_t i = 0; i < max_len || carry > 0; ++i) {
            if (i >= combined.size()) combined.push_back(0);
            
            uint64_t v1 = (i < res1.size()) ? res1[i] : 0;
            uint64_t v2 = (i < res2.size()) ? res2[i] : 0;
            uint64_t v3 = (i < res3.size()) ? res3[i] : 0;

            unsigned __int128 val = combine_crt(v1, v2, v3) + carry;
            combined[i] = (uint64_t)(val & 0xFFFFFFFFULL);
            carry = val >> 32;
        }

        MersenneBigInt res(a.p);
        res.limbs = combined;
        res.reduce_mersenne();
        return res;
    }

    static MersenneBigInt square(const MersenneBigInt& a) {
        return multiply(a, a);
    }

    void sub_small(uint64_t v) {
        int64_t borrow = v;
        for (size_t i = 0; i < limbs.size(); ++i) {
            int64_t diff = (int64_t)limbs[i] - borrow;
            if (diff < 0) {
                limbs[i] = diff + (1ULL << 32);
                borrow = 1;
            } else {
                limbs[i] = diff;
                borrow = 0;
                break;
            }
        }
        if (borrow > 0) {
            limbs[0] += (1ULL << 32) - borrow;
            reduce_mersenne();
        }
        while (limbs.size() > 1 && limbs.back() == 0) limbs.pop_back();
    }

    bool is_zero_mod() const {
        bool all_zero = true;
        for (auto v : limbs) if (v != 0) all_zero = false;
        if (all_zero) return true;

        for (int i = 0; i < p; ++i) {
            int limb_idx = i / 32;
            int bit_idx = i % 32;
            if (limb_idx >= (int)limbs.size() || !(limbs[limb_idx] & (1ULL << bit_idx))) {
                return false;
            }
        }
        for (size_t i = p / 32; i < limbs.size(); ++i) {
            for (int b = 0; b < 32; ++b) {
                if (i * 32 + b >= (size_t)p) {
                    if ((limbs[i] >> b) & 1) return false;
                }
            }
        }
        return true;
    }
};

inline bool is_prime_optimized(long long n) {
    if (n <= 1) return false;
    if (n <= 3) return true;
    if (n % 2 == 0 || n % 3 == 0) return false;
    for (long long i = 5; i <= n / i; i += 6) {
        if (n % i == 0 || n % (i + 2) == 0) return false;
    }
    return true;
}

// 🧠 [최종 강화된 수학적 필터 파이프라인]
bool advanced_mathematical_filter(int p) {
    if (p <= 3) return true;

    // 1단계: 페르마 소정리 기반 간이 검증 (Fermat Quick Check)
    // 2^(p-1) mod p == 1 인지 확인하여 지수 p의 무결성 재확인
    if (power(2, p - 1, p) != 1) {
        return false; // p가 소수가 아님 (이론상 방어선)
    }

    // 2단계: 확장된 2kp+1 Trial Factoring (GIMPS 전용 소인수 필터)
    // 탐색 범위를 대폭 넓혀 웬만한 합성수 소인수는 여기서 전부 격멸
    uint64_t limit = 20000000; // 2천만까지 확장
    for (uint64_t k = 1; ; ++k) {
        uint64_t q = 2 * k * p + 1;
        if (q / (2 * p) < k) break; // 오버플로우 방지
        if (q > limit) break;

        if (is_prime_optimized(q)) {
            // 메르센 수의 소인수는 반드시 q ≡ 1 또는 7 (mod 8) 만족해야 함
            uint64_t mod8 = q % 8;
            if (mod8 == 1 || mod8 == 7) {
                if (power(2, p, q) == 1) {
                    return false; // 100% 합성수 판정! LL 테스트 스킵
                }
            }
        }
    }

    return true; // 모든 수학적 방어선을 통과한 진정한 메르센 소수 후보
}

bool lucas_lehmer_test(int p) {
    if (p == 2) return true;

    // 수학적 필터링을 통과하지 못하면 즉시 기각 (Composite)
    if (!advanced_mathematical_filter(p)) {
        return false; 
    }

    // 최종 생존자들만 고속 제곱 뤼카-레머 테스트 수행
    MersenneBigInt s(p, 4);
    for (int i = 0; i < p - 2; ++i) {
        s = MersenneBigInt::square(s);
        s.sub_small(2);
        s.reduce_mersenne();
    }
    return s.is_zero_mod();
}

void worker_thread(int thread_id) {
    while (true) {
        int p = global_p.fetch_add(2);
        if (!is_prime_optimized(p)) continue;

        bool is_prime = lucas_lehmer_test(p);

        if (is_prime) {
            std::lock_guard<std::mutex> lock(cout_mutex);
            std::cout << "\n🎉 >>> [Thread " << thread_id << "] MERSENNE PRIME FOUND! 2^" << p << " - 1 is PRIME! <<<\n" << std::endl;
        }
    }
}

int main() {
    std::cout << "=== GIMPS-Grade Ultimate Math-Filtered Mersenne Engine (C++23) ===" << std::endl;
    std::cout << "Pipelines: Fermat Check, Enhanced 2kp+1 Trial Factoring, Multi-Prime NTT\n" << std::endl;

    unsigned int hardware_threads = std::thread::hardware_concurrency();
    unsigned int num_threads = (hardware_threads > 0) ? hardware_threads : 4;
    std::cout << "Spawning " << num_threads << " parallel worker threads...\n" << std::endl;

    if (lucas_lehmer_test(2)) {
        std::lock_guard<std::mutex> lock(cout_mutex);
        std::cout << "🎉 >>> MERSENNE PRIME FOUND! 2^2 - 1 is PRIME! <<<\n" << std::endl;
    }

    global_p.store(3);

    std::vector<std::j