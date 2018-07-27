#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>
#include <vector>

struct Const_arrays {
	static const std::vector<int> belt;
	static const std::vector<double> fact;
	static const std::vector<double> dfact;
	static const std::vector<std::vector<double> > binom;
	static const std::vector<std::vector<double> > omega;
	static const std::vector<unsigned> crt_siz;

	static constexpr size_t belt_s = 500;
	static constexpr size_t fact_s = 15;
	static constexpr size_t dfact_s = 15;
	static constexpr size_t binom_s = 30;
	static constexpr size_t omega_s = 30;
	static constexpr size_t crt_siz_s = 11;
};

#endif // CONSTANTS_H
