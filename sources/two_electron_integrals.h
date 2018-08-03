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
	using index = typename Eigen::Tensor<type, 4>::Index;
	using index_array = std::array<index, 4>;
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
			neg_cubes[i].resize(i, i, i);
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

	inline const size_t get_size() const
	{
		return size;
	}

	inline type coef(const index &i, const index &j, const index &k, const index &l) const
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
			f ? pos_cubes[a[3]]({a[0], a[1], a[2]}) = conj(val) : pos_cubes[a[3]]({a[0], a[1], a[2]}) = val;
	}

	inline Eigen::Matrix<type, Eigen::Dynamic, Eigen::Dynamic>
	contract(const Eigen::Matrix<type, Eigen::Dynamic, 1> &vec1,
			 const Eigen::Matrix<type, Eigen::Dynamic, 1> &vec2,
			 const index &i1, const index &i2)
	{
		assert(vec1.size() == size && vec2.size() == size);
		assert(i1 == 0 || i1 == 2);
		assert(i2 == 1 || i2 == 3);

		Eigen::Matrix<type, Eigen::Dynamic, Eigen::Dynamic> res(size, size);
		res.setZero();

		switch (i1)
		{
		case 0:
			switch (i2)
			{
			case 1:
				for (unsigned i = 0; i < size; i++)
					for (unsigned j = 0; j < size; j++)
						for (unsigned k = 0; k < size; k++)
							for (unsigned l = 0; l < size; l++)
								res(i, j) += conj(vec1(k)) * coef(k, l, i, j) * vec2(l);
				break;
			case 3:
				for (unsigned i = 0; i < size; i++)
					for (unsigned j = 0; j < size; j++)
						for (unsigned k = 0; k < size; k++)
							for (unsigned l = 0; l < size; l++)
								res(i, j) += conj(vec1(k)) * coef(k, j, i, l) * vec2(l);
				break;
			}
			break;
		case 2:
			switch (i2)
			{
			case 1:
				for (unsigned i = 0; i < size; i++)
					for (unsigned j = 0; j < size; j++)
						for (unsigned k = 0; k < size; k++)
							for (unsigned l = 0; l < size; l++)
								res(i, j) += conj(vec1(k)) * coef(i, l, k, j) * vec2(l);
				break;
			case 3:
				for (unsigned i = 0; i < size; i++)
					for (unsigned j = 0; j < size; j++)
						for (unsigned k = 0; k < size; k++)
							for (unsigned l = 0; l < size; l++)
								res(i, j) += conj(vec1(k)) * coef(i, j, k, l) * vec2(l);
				break;
			}
			break;
		}
		return res;
	}

	void print(std::ostream &os)
	{
		for (size_t i = 0; i < size; ++i)
			for (size_t j = 0; j < size; ++j)
				for (size_t k = 0; k < size; ++k)
					for (size_t l = 0; l < size; ++l)
					{
						os << "  " << i + 1 << " " << j + 1 << " " << k + 1 << " " << l + 1 << "\n";
						os << coef(i, j, k, l) << "\n";
					}
		os << "\n";
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
			return std::make_pair(index_array({i, j, k, l}), false);
		else if (is_unique_combination(k, l, i, j))
			return std::make_pair(index_array({k, l, i, j}), false);
		else if (is_unique_combination(j, i, l, k))
			return std::make_pair(index_array({j, i, l, k}), true);
		else
			return std::make_pair(index_array({l, k, j, i}), true);
	}

  private:
	std::vector<Eigen::Tensor<type, 3>> pos_cubes;
	std::vector<Eigen::Tensor<type, 3>> neg_cubes;

	size_t size;
};

using Tensor_2Ecd = Tensor_2E<std::complex<double>>;
using Tensor_2Ed = Tensor_2E<double>;

#endif // TWO_ELECTRON_INTEGRALS_H
