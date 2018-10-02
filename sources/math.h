#pragma once

#include <boost/math/special_functions/factorials.hpp>

using boost::math::factorial;

namespace math {

constexpr int min_to_m(int m) {
    return m % 2 == 0 ? 1 : -1;
}

double omega(int i, int j) {
    return std::sqrt((2.0 * i + 1) * factorial<double>(i - j) /
                     (2.0 * factorial<double>(i + j)));
}

}
