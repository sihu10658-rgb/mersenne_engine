# mersenne_engine
# 🚀 GIMPS-Grade Parallel Mersenne Prime Engine (C++23)

A high-performance, multi-threaded Mersenne Prime search engine written in modern **C++23**, inspired by GIMPS (Great Internet Mersenne Prime Search) architecture.

---

## ⚡ Key Features

1. **Multi-Prime NTT (Number Theoretic Transform)**
   - Bypasses single-modulus overflow limits using three large NTT primes (`MOD1`, `MOD2`, `MOD3`) combined with **Garner's algorithm / Chinese Remainder Theorem (CRT)** for 128-bit precision handling.
2. **Advanced Mathematical Filtering Pipeline**
   - **Fermat Quick Check**: Eliminates invalid exponents instantly.
   - **Enhanced $2kp+1$ Trial Factoring**: Explores Mersenne factor properties ($q \equiv \pm 1 \pmod 8$) up to $20,000,000$ to cull composite candidates before heavy tests.
3. **Optimized Arithmetic & Squaring**
   - Implements fast squaring ($O(N \log N)$ complexity via NTT) and bit-folding reduction (`2^p - 1` modulus reduction).
4. **Massive Parallelization**
   - Uses C++23 `std::jthread` and atomic index distribution to max out all available CPU cores seamlessly.

---

## 🛠️ Compilation & Usage

Compile with extreme performance flags (GCC / Clang):

```bash
g++ -O3 -march=native -mtune=native -flto -funroll-loops -fomit-frame-pointer -fno-stack-protector -std=c++23 -pthread mersenne_engine.cpp -o mersenne_engine
