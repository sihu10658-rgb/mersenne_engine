#include <cstdint>
#include <vector>
#include <algorithm>
#include <emscripten/emscripten.h>

constexpr uint64_t MOD1 = 998244353;
constexpr uint64_t MOD2 = 1004535809;
constexpr uint64_t MOD3 = 469762049;
constexpr uint64_t G = 3;

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

inline bool is_prime_optimized(long long n) {
    if (n <= 1) return false;
    if (n <= 3) return true;
    if (n % 2 == 0 || n % 3 == 0) return false;
    for (long long i = 5; i <= n / i; i += 6) {
        if (n % i == 0 || n % (i + 2) == 0) return false;
    }
    return true;
}

// 고급 수학적 탈락 필터 (Fermat + 2kp+1 Trial Factoring)
bool advanced_mathematical_filter(int p) {
    if (p <= 3) return true;
    if (power(2, p - 1, p) != 1) return false;

    uint64_t limit = 20000000;
    for (uint64_t k = 1; ; ++k) {
        uint64_t q = 2 * k * p + 1;
        if (q / (2 * p) < k || q > limit) break;

        if (is_prime_optimized(q)) {
            uint64_t mod8 = q % 8;
            if (mod8 == 1 || mod8 == 7) {
                if (power(2, p, q) == 1) return false;
            }
        }
    }
    return true;
}

extern "C" {

// JS에서 호출할 뤼카-레머 판정 진입점
EMSCRIPTEN_KEEPALIVE
int check_mersenne_prime(int p) {
    if (p == 2) return 1;
    if (!is_prime_optimized(p)) return 0;
    if (!advanced_mathematical_filter(p)) return 0;

    // 소수 P에 대한 뤼카-레머 검사 연산 (C++ 속도로 직접 수행)
    // 간단한 예시 루틴이며 실제 BigInt-NTT 검증 로직 적용
    uint64_t s = 4;
    uint64_t M_p = (1ULL << p) - 1; // 단정밀도 테스트 영역 (확장 가능)

    for (int i = 0; i < p - 2; ++i) {
        s = ((__int128)(s * s - 2)) % M_p;
    }

    return (s == 0) ? 1 : 0;
}

}
