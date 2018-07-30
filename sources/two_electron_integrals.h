#ifndef TWO_ELECTRON_INTEGRALS_H
#define TWO_ELECTRON_INTEGRALS_H

#include <complex>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <array>
#include <utility>
#include <initializer_list>
#include <stdexcept>

#include <eigen3/unsupported/Eigen/CXX11/Tensor>

using Two_electron_ints_C = Eigen::Tensor<std::complex<double>, 4>;
using Two_electron_ints_D = Eigen::Tensor<double, 4>;

void read_two_el_from_binary(Two_electron_ints_C& integrals, const std::string file_name, const unsigned basis_size, bool print);

void read_two_el_from_binary_slice(Two_electron_ints_D& integrals, const std::string file_name, const unsigned basis_size, bool print);

void print_two_el_ints(Two_electron_ints_D &integrals);

void print_two_el_ints(Two_electron_ints_C &integrals);

void two_el_int_compare(Two_electron_ints_C &ints_1, Two_electron_ints_D &ints_2, const int b_l);

template<class type>
class Tensor_2E {
public:
	using index_array = std::array<Eigen::Tensor::Index, 4>>;
	using index = Eigen::Tensor::Index;
	using conjugation_flag = bool;

	enum position {
		positive, negaive
	};

	void resize(const size_t& s) {
		assert(s >= 0);
		size = s;
		pos_cubes.resize(size);
		neg_cubes.resize(size);
		for (int i = 0; i < size; ++i) {
			pos_cubes[i].resize(size - i, size - i, size - i);
			neg_cubes[i].resize(i, i, i, i);
		}
	}

	void set_zero() {
		for (int i = 0; i < size; ++i) {
			pos_cubes[i].setZero();
			neg_cubes[i].setZero();
		}
	}

	inline const auto get_size() const {
		return size;
	}

	inline auto& coef(const index& i, const index& j, const index& k, const index& l) {
		auto unq_perm = get_unq_combination(i, j, k, l);
		conjugation_flag f = std::get<1>(unq_perm);
		index_array a = std::get<0>(unq_perm);
		if (a[0] < a[3]) {
			auto val = neg_cubes[a[3]]( { a[0], a[1], a[2] });
			return f ? conj(val) : val;
		} else {
			auto val = pos_cubes[a[3]]( { a[0], a[1], a[2] });
			return f ? conj(val) : val;
		}
	}

	auto& operator()(const index& i, const index& j, const index& k, const index& l) {
		return coef(i, j, k, l);
	}

	void read_from_binary(const std::string& file_name, const unsigned& basis_size) {
		std::ifstream file(file_name, std::ios::in | std::ios::binary | std::ios::ate);
		if (!file.is_open())
			throw std::runtime_error("Cannot open the two electron binary file.");

		resize(basis_size);
		set_zero();

		std::streampos file_size = file.tellg();
		unsigned long number_two_el = file_size * sizeof(char);
		number_two_el /= 2 * sizeof(type) + 4 * sizeof(unsigned short int);
		file.seekg(0, std::ios::beg);

		type re_data, im_data;
		unsigned short int indi, indj, indk, indl;

		for (unsigned long i = 0; i < number_two_el; ++i) {
			file.read(reinterpret_cast<char*>(&indi), sizeof(unsigned short int));
			file.read(reinterpret_cast<char*>(&indj), sizeof(unsigned short int));
			file.read(reinterpret_cast<char*>(&indk), sizeof(unsigned short int));
			file.read(reinterpret_cast<char*>(&indl), sizeof(unsigned short int));

			file.read(reinterpret_cast<char*>(&re_data), sizeof(type));
			file.read(reinterpret_cast<char*>(&im_data), sizeof(type));

			if (!is_unique_combination(--indi, --indj, --indk, --indl))
				break;

			operator()(indi, indj, indk, indl) = std::complex<type>(re_data, im_data);
		}
		file.close();
	}

	void print() {
		for (size_t i = 0; i < size; ++i)
			for (size_t j = 0; j < size; ++j)
				for (size_t k = 0; k < size; ++k)
					for (size_t l = 0; l < size; ++l) {
						std::cout << "  " << i + 1 << " " << j + 1 << " " << k + 1 << " " << l + 1 << "\n";
						std::cout << operator()(i, j, k, l) << "\n";
					}
		std::cout << "\n";
	}

protected:
	bool is_unique_combination(const index& i, const index& j, const index& k, const index& l) const {
		if (i >= l && j >= l && k >= l)
			return true;
		else if (i < l && j < l && k < l)
			return true;
		else
			return false;
	}

	std::pair<index_array, conjugation_flag> get_unq_combination(const index& i, const index& j, const index& k, const index& l) const {
		if (is_unique_combination(i, j, k, l))
			return std::make_pair( { i, j, k, l }, false);
		else if (is_unique_combination(k, l, i, j))
			return std::make_pair( { k, l, i, j }, false);
		else if (is_unique_combination(j, i, l, k))
			return std::make_pair( { j, i, l, k }, true);
		else
			return std::make_pair( { l, k, j, i }, true);
	}
private:
	std::vector<Eigen::Tensor<std::complex<type>, 3>> pos_cubes;
	std::vector<Eigen::Tensor<std::complex<type>, 3>> neg_cubes;

	size_t size;
};

#endif // TWO_ELECTRON_INTEGRALS_H
