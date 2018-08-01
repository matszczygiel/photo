#ifndef TWO_ELECTRON_INTEGRALS_H
#define TWO_ELECTRON_INTEGRALS_H

#include <complex>
#include <iostream>
#include <vector>
#include <array>
#include <utility>
#include <cassert>

#include <eigen3/unsupported/Eigen/CXX11/Tensor>

template <class type>
class Tensor_2E
{
  public:
	using index_array = std::array<Eigen::Tensor::Index, 4>> ;
	using index = Eigen::Tensor::Index;
	using conjugation_flag = bool;

	enum position
	{
		positive,
		negaive
	};

	void resize(const size_t &s)
	{
		assert(s >= 0);
		size = s;
		pos_cubes.resize(size);
		neg_cubes.resize(size);
		for (int i = 0; i < size; ++i)
		{
			pos_cubes[i].resize(size - i, size - i, size - i);
			neg_cubes[i].resize(i, i, i, i);
		}
	}

	void set_zero()
	{
		for (int i = 0; i < size; ++i)
		{
			pos_cubes[i].setZero();
			neg_cubes[i].setZero();
		}
	}

	inline const auto get_size() const
	{
		return size;
	}

	inline auto coef(const index &i, const index &j, const index &k, const index &l) const
	{
		assert(i < size && j < size && k < size && l < size);
		auto unq_perm = get_unq_combination(i, j, k, l);
		conjugation_flag f = std::get<1>(unq_perm);
		index_array a = std::get<0>(unq_perm);
		if (a[0] < a[3])
		{
			auto val = neg_cubes[a[3]]({a[0], a[1], a[2]});
			return f ? conj(val) : val;
		}
		else
		{
			auto val = pos_cubes[a[3]]({a[0], a[1], a[2]});
			return f ? conj(val) : val;
		}
	}

	inline void assign(const index &i, const index &j, const index &k, const index &l, const std::complex<double> &val)
	{
		assert(i < size && j < size && k < size && l < size);
		auto unq_perm = get_unq_combination(i, j, k, l);
		conjugation_flag f = std::get<1>(unq_perm);
		index_array a = std::get<0>(unq_perm);
		if (a[0] < a[3])
			f ? neg_cubes[a[3]]({a[0], a[1], a[2]}) = conj(val) : neg_cubes[a[3]]({a[0], a[1], a[2]}) = val;
		else
			f ? pos_cubes[a[3]]({a[0], a[1], a[2]}) = conj(val) : npos_cubes[a[3]]({a[0], a[1], a[2]}) = val;
	}

	void print()
	{
		for (size_t i = 0; i < size; ++i)
			for (size_t j = 0; j < size; ++j)
				for (size_t k = 0; k < size; ++k)
					for (size_t l = 0; l < size; ++l)
					{
						std::cout << "  " << i + 1 << " " << j + 1 << " " << k + 1 << " " << l + 1 << "\n";
						std::cout << coef(i, j, k, l) << "\n";
					}
		std::cout << "\n";
	}

  protected:
	bool is_unique_combination(const index &i, const index &j, const index &k, const index &l) const
	{
		if (i >= l && j >= l && k >= l)
			return true;
		else if (i < l && j < l && k < l)
			return true;
		else
			return false;
	}

	std::pair<index_array, conjugation_flag>
	get_unq_combination(const index &i, const index &j, const index &k, const index &l) const
	{
		if (is_unique_combination(i, j, k, l))
			return std::make_pair({i, j, k, l}, false);
		else if (is_unique_combination(k, l, i, j))
			return std::make_pair({k, l, i, j}, false);
		else if (is_unique_combination(j, i, l, k))
			return std::make_pair({j, i, l, k}, true);
		else
			return std::make_pair({l, k, j, i}, true);
	}

  private:
	std::vector<Eigen::Tensor<type>, 3>> pos_cubes;
	std::vector<Eigen::Tensor<type>, 3>> neg_cubes;

	size_t size;
};

using Tensor_2Ecd = Tensor_2E<std::complex<double>>;

#endif // TWO_ELECTRON_INTEGRALS_H
