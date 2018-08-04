#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>
#include <vector>

namespace Const_arrays
{
extern const std::vector<int> belt;
extern const std::vector<double> fact;
extern const std::vector<double> dfact;
extern const std::vector<std::vector<double>> binom;
extern const std::vector<std::vector<double>> omega;
extern const std::vector<unsigned> crt_siz;

constexpr std::size_t belt_s = 500;
constexpr std::size_t fact_s = 15;
constexpr std::size_t dfact_s = 15;
constexpr std::size_t binom_s = 30;
constexpr std::size_t omega_s = 30;
constexpr std::size_t crt_siz_s = 9;
}; // namespace Const_arrays

#endif // CONSTANTS_H
