#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <array>
#include <cmath>

namespace Const_arrays {
constexpr int belt_s    = 500;
constexpr int fact_s    = 15;
constexpr int dfact_s   = 15;
constexpr int binom_s   = 30;
constexpr int omega_s   = 30;

const std::array<int, belt_s> belt = []() {
    std::array<int, belt_s> belt;
    belt[0] = 1;
    for (unsigned i = 1; i < belt_s; ++i)
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
    for (unsigned i = 1; i < dfact_s; ++i)
        dfact[i + 1] = dfact[i] * (2 * i + 1);
    return dfact;
}();

const std::array<std::array<double, binom_s>, binom_s> binom = []() {
    std::array<std::array<double, binom_s>, binom_s> binom;
    for (unsigned i = 0; i < binom_s; ++i) {
        for (unsigned j = 0; j <= i; ++j) {
            binom[i][j] = fact[i];
            binom[i][j] /= fact[j] * fact[i - j];
            binom[j][i] = binom[i][j];
        }
    }
    return binom;
}();

const std::array<std::array<double, omega_s>, omega_s> omega = []() {
    std::array<std::array<double, omega_s>, omega_s> omega;
    for (unsigned i = 0; i < omega_s; i++) {
        for (unsigned j = 0; j <= i; j++)
            omega[i][j] = std::sqrt((2.0 * i + 1) * fact[i - j] / (2.0 * fact[i + j]));
    }
    return omega;
}();

};  // namespace Const_arrays

#endif  // CONSTANTS_H
