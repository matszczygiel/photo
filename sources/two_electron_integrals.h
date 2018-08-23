#ifndef TWO_ELECTRON_INTEGRALS_H
#define TWO_ELECTRON_INTEGRALS_H

#include <algorithm>
#include <array>
#include <cassert>
#include <complex>
#include <iostream>
#include <utility>

#include <eigen3/Eigen/Dense>

template <class type>
class Tensor_2E {
   private:
    type**** pos_cubes = nullptr;
    type**** neg_cubes = nullptr;

    unsigned size = 0;

   public:
    typedef int index;
    typedef std::array<index, 4> index_array;
    typedef bool conjugation_flag;

    enum position { positive,
                    negaive };

    Tensor_2E() = default;

    inline const unsigned get_size() const { return size; }

    void free() {
        if (pos_cubes == nullptr && neg_cubes == nullptr) {
            size = 0;
            return;
        }

        for (unsigned i = 0; i < size; ++i) {
            for (unsigned it1 = 0; it1 < size - i; ++it1) {
                for (unsigned it2 = 0; it2 < size - i; ++it2)
                    delete[] pos_cubes[i][it1][it2];

                delete[] pos_cubes[i][it1];
            }

            for (unsigned it1 = 0; it1 < i; ++it1) {
                for (unsigned it2 = 0; it2 < i; ++it2)
                    delete[] neg_cubes[i][it1][it2];

                delete[] neg_cubes[i][it1];
            }

            delete[] pos_cubes[i];
            delete[] neg_cubes[i];
        }
        delete[] pos_cubes;
        delete[] neg_cubes;
        size = 0;
    }

    ~Tensor_2E() { free(); }

    void resize(const int& s) {
        assert(s >= 0);

        free();
        size = s;

        pos_cubes = new type***[size];
        neg_cubes = new type***[size];

        for (unsigned i = 0; i < size; ++i) {
            pos_cubes[i] = new type**[size - i];

            for (unsigned it1 = 0; it1 < size - i; ++it1) {
                pos_cubes[i][it1] = new type*[size - i];
                for (unsigned it2 = 0; it2 < size - i; ++it2)
                    pos_cubes[i][it1][it2] = new type[size - i];
            }

            neg_cubes[i] = new type**[i];

            for (unsigned it1 = 0; it1 < i; ++it1) {
                neg_cubes[i][it1] = new type*[i];
                for (unsigned it2 = 0; it2 < i; ++it2)
                    neg_cubes[i][it1][it2] = new type[i];
            }
        }
    }

    void zero() {
        for (unsigned i = 0; i < size; ++i) {
            for (unsigned it1 = 0; it1 < size - i; ++it1)
                for (unsigned it2 = 0; it2 < size - i; ++it2)
                    std::fill(&pos_cubes[i][it1][it2][0],
                              &pos_cubes[i][it1][it2][size - i], 0);

            for (unsigned it1 = 0; it1 < i; ++it1)
                for (unsigned it2 = 0; it2 < i; ++it2)
                    std::fill(&neg_cubes[i][it1][it2][0], &neg_cubes[i][it1][it2][i], 0);
        }
    }

    inline type coef(const index& i,
                     const index& j,
                     const index& k,
                     const index& l) const {
        assert(i < size && j < size && k < size && l < size);
        auto unq_perm      = get_unq_combination(i, j, k, l);
        conjugation_flag f = std::get<1>(unq_perm);
        index_array a      = std::get<0>(unq_perm);
        if (a[0] < a[3]) {
            auto val = neg_cubes[a[3]][a[0]][a[1]][a[2]];
            return f ? conj(val) : val;
        } else {
            auto val = pos_cubes[a[3]][a[0] - a[3]][a[1] - a[3]][a[2] - a[3]];
            return f ? conj(val) : val;
        }
    }

    inline void assign(const index& i,
                       const index& j,
                       const index& k,
                       const index& l,
                       const type& val) {
        assert(i < size && j < size && k < size && l < size);
        auto unq_perm      = get_unq_combination(i, j, k, l);
        conjugation_flag f = std::get<1>(unq_perm);
        index_array a      = std::get<0>(unq_perm);
        assert((a[0] < a[3] && a[1] < a[3] && a[2] < a[3]) ||
               (a[0] >= a[3] && a[1] >= a[3] && a[2] >= a[3]));
        if (a[0] < a[3])
            f ? neg_cubes[a[3]][a[0]][a[1]][a[2]] = conj(val)
              : neg_cubes[a[3]][a[0]][a[1]][a[2]] = val;
        else
            f ? pos_cubes[a[3]][a[0] - a[3]][a[1] - a[3]][a[2] - a[3]] = conj(val)
              : pos_cubes[a[3]][a[0] - a[3]][a[1] - a[3]][a[2] - a[3]] = val;
    }

    inline Eigen::Matrix<type, Eigen::Dynamic, Eigen::Dynamic> contract(
        const Eigen::Matrix<type, Eigen::Dynamic, 1>& vec1,
        const Eigen::Matrix<type, Eigen::Dynamic, 1>& vec2,
        const index& i1,
        const index& i2) {
        assert(vec1.size() == size && vec2.size() == size);
        assert(i1 == 0 || i1 == 2);
        assert(i2 == 1 || i2 == 3);

        Eigen::Matrix<type, Eigen::Dynamic, Eigen::Dynamic> res(size, size);
        res.setZero();

        switch (i1) {
            case 0:
                switch (i2) {
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
                switch (i2) {
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

    void print(std::ostream& os) {
        for (unsigned i = 0; i < size; ++i)
            for (unsigned j = 0; j < size; ++j)
                for (unsigned k = 0; k < size; ++k)
                    for (unsigned l = 0; l < size; ++l) {
                        os << "  " << i << " " << j << " " << k << " " << l << "\n";
                        os << coef(i, j, k, l) << "\n";
                    }
        os << "\n";
    }

   protected:
    bool is_unique_combination(const index& i,
                               const index& j,
                               const index& k,
                               const index& l) const {
        if (i >= l && j >= l && k >= l)
            return true;
        else if (i < l && j < l && k < l)
            return true;
        else
            return false;
    }

    std::pair<index_array, conjugation_flag> get_unq_combination(
        const index& i,
        const index& j,
        const index& k,
        const index& l) const {
        if (is_unique_combination(i, j, k, l))
            return std::make_pair(index_array({i, j, k, l}), false);
        else if (is_unique_combination(k, l, i, j))
            return std::make_pair(index_array({k, l, i, j}), false);
        else if (is_unique_combination(j, i, l, k))
            return std::make_pair(index_array({j, i, l, k}), true);
        else
            return std::make_pair(index_array({l, k, j, i}), true);
    }
};

typedef Tensor_2E<std::complex<double>> Tensor_2Ecd;

#endif  // TWO_ELECTRON_INTEGRALS_H
