#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <array>
#include <cmath>


constexpr double au_to_meters      = 5.2917721092e-11;
constexpr double barn_to_sqrmeters = 1.0e-28;

constexpr double au_to_barns = au_to_meters * au_to_meters / barn_to_sqrmeters;

constexpr int min_to_m(int m) {
    return m % 2 == 0 ? 1 : -1;
}

constexpr int omega(int i, int j) {
    return std::sqrt((2.0 * i + 1) * fact(i - j) / (2.0 * fact(i + j)));
}

namespace Const_arrays {
constexpr int omega_s = 30;



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
