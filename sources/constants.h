#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <array>
#include <cmath>

constexpr int min_to_m(int m) {
    return m % 2 == 0 ? 1 : -1;
}

constexpr int fact(int n) {
    return n <= 1 ? 1 : (n * fact(n - 1));
}

constexpr int dfact(int n) {
    return n <= 1 ? 1 : ((2 * n - 1) * dfact(n - 1));
}

constexpr int binom(int n, int k) {
    return fact(n) / fact(k) / fact(n - k);
}

constexpr int omega(int i, int j) {
    return std::sqrt((2.0 * i + 1) * fact(i - j) / (2.0 * fact(i + j)));
}

namespace Const_arrays {
constexpr int belt_s  = 500;
constexpr int fact_s  = 15;
constexpr int dfact_s = 15;
constexpr int binom_s = 30;
constexpr int omega_s = 30;

const std::array<int, belt_s> belt = []() {
    std::array<int, belt_s> belt;
    belt[0] = 1;
    for (int i = 1; i < belt_s; ++i)
        belt[i] = -belt[i - 1];
    return belt;
}();

const std::array<double, fact_s> fact = []() {
    std::array<double, fact_s> fact;
    fact[0] = 1.0;
    fact[1] = 1.0;
    for (int i = 2; i < fact_s; ++i)
        fact[i] = fact[i - 1] * i;
    return fact;
}();

const std::array<double, dfact_s> dfact = []() {
    std::array<double, dfact_s> dfact;
    dfact[0] = 1.0;
    dfact[1] = 1.0;
    for (int i = 1; i < dfact_s; ++i)
        dfact[i + 1] = dfact[i] * (2 * i + 1);
    return dfact;
}();

const std::array<std::array<double, binom_s>, binom_s> binom = []() {
    std::array<std::array<double, binom_s>, binom_s> binom;
    for (int i = 0; i < binom_s; ++i) {
        for (int j = 0; j <= i; ++j) {
            binom[i][j] = fact[i];
            binom[i][j] /= fact[j] * fact[i - j];
            binom[j][i] = binom[i][j];
        }
    }
    return binom;
}();

const std::array<std::array<double, omega_s>, omega_s> omega = []() {
    std::array<std::array<double, omega_s>, omega_s> omega;
    for (int i = 0; i < omega_s; i++) {
        for (int j = 0; j <= i; j++)
            omega[i][j] = std::sqrt((2.0 * i + 1) * fact[i - j] / (2.0 * fact[i + j]));
    }
    return omega;
}();

};  // namespace Const_arrays

#endif  // CONSTANTS_H
