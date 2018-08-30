#ifndef TWO_ELECTRON_INTEGRALS_H
#define TWO_ELECTRON_INTEGRALS_H

#include <algorithm>
#include <array>
#include <cassert>
#include <complex>
#include <iostream>
#include <utility>

#include <eigen3/Eigen/Dense>

template <typename T>
class Tensor_2E {
   private:
    std::complex<T>**** pos_cubes = nullptr;
    std::complex<T>**** neg_cubes = nullptr;

    int size = 0;

   public:
    typedef int index;
    typedef std::array<index, 4> index_array;
    typedef bool conjugation_flag;

    enum position { positive,
                    negaive };

    void free() {
        if (pos_cubes == nullptr && neg_cubes == nullptr) {
            size = 0;
            return;
        }

        for (int i = 0; i < size; ++i) {
            for (int it1 = 0; it1 < size - i; ++it1) {
                for (int it2 = 0; it2 < size - i; ++it2)
                    delete[] pos_cubes[i][it1][it2];

                delete[] pos_cubes[i][it1];
            }

            for (int it1 = 0; it1 < i; ++it1) {
                for (int it2 = 0; it2 < i; ++it2)
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

    void resize(const int& s) {
        assert(s >= 0);

        free();
        size = s;

        pos_cubes = new std::complex<T>***[size];
        neg_cubes = new std::complex<T>***[size];

        for (int i = 0; i < size; ++i) {
            pos_cubes[i] = new std::complex<T>**[size - i];

            for (int it1 = 0; it1 < size - i; ++it1) {
                pos_cubes[i][it1] = new std::complex<T>*[size - i];
                for (int it2 = 0; it2 < size - i; ++it2)
                    pos_cubes[i][it1][it2] = new std::complex<T>[size - i];
            }

            neg_cubes[i] = new std::complex<T>**[i];

            for (int it1 = 0; it1 < i; ++it1) {
                neg_cubes[i][it1] = new std::complex<T>*[i];
                for (int it2 = 0; it2 < i; ++it2)
                    neg_cubes[i][it1][it2] = new std::complex<T>[i];
            }
        }
    }

    Tensor_2E() = default;

    Tensor_2E(Tensor_2E& other) {
        resize(other.size);
        for (int i = 0; i < size; ++i) {
            for (int it1 = 0; it1 < size - i; ++it1)
                for (int it2 = 0; it2 < size - i; ++it2)
                    std::copy(&other.pos_cubes[i][it1][it2][0],
                              &other.pos_cubes[i][it1][it2][size - i],
                              pos_cubes[i][it1][it2]);

            for (int it1 = 0; it1 < i; ++it1)
                for (int it2 = 0; it2 < i; ++it2)
                    std::copy(&other.neg_cubes[i][it1][it2][0],
                              &other.neg_cubes[i][it1][it2][i],
                              neg_cubes[i][it1][it2]);
        }
    }

    Tensor_2E& operator=(Tensor_2E& other) {
        if (this != &other) {
            free();
            resize(other.size);

            for (int i = 0; i < size; ++i) {
                for (int it1 = 0; it1 < size - i; ++it1)
                    for (int it2 = 0; it2 < size - i; ++it2)
                        std::copy(&other.pos_cubes[i][it1][it2][0],
                                  &other.pos_cubes[i][it1][it2][size - i],
                                  pos_cubes[i][it1][it2]);

                for (int it1 = 0; it1 < i; ++it1)
                    for (int it2 = 0; it2 < i; ++it2)
                        std::copy(&other.neg_cubes[i][it1][it2][0],
                                  &other.neg_cubes[i][it1][it2][i],
                                  neg_cubes[i][it1][it2]);
            }
        }
        return *this;
    }

    Tensor_2E& operator=(Tensor_2E&& other) {
        if (this != &other) {
            free();
            pos_cubes       = std::move(other.pos_cubes);
            neg_cubes       = std::move(other.neg_cubes);
            size            = std::move(other.size);
            other.neg_cubes = nullptr;
            other.pos_cubes = nullptr;
            other.size      = 0;
        }
        return *this;
    }

    Tensor_2E(Tensor_2E&& other)
        : pos_cubes(std::move(other.pos_cubes)),
          neg_cubes(std::move(other.neg_cubes)),
          size(std::move(other.size)) {
        other.neg_cubes = nullptr;
        other.pos_cubes = nullptr;
        other.size      = 0;
    }

    ~Tensor_2E() { free(); }

    inline const int get_size() const { return size; }

    void zero() {
        for (int i = 0; i < size; ++i) {
            for (int it1 = 0; it1 < size - i; ++it1)
                for (int it2 = 0; it2 < size - i; ++it2)
                    std::fill(&pos_cubes[i][it1][it2][0],
                              &pos_cubes[i][it1][it2][size - i], 0);

            for (int it1 = 0; it1 < i; ++it1)
                for (int it2 = 0; it2 < i; ++it2)
                    std::fill(&neg_cubes[i][it1][it2][0],
                              &neg_cubes[i][it1][it2][i], 0);
        }
    }

    inline std::complex<T> coef(const index& i,
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
                       const std::complex<T>& val) {
        assert(i < size && j < size && k < size && l < size);
        auto unq_perm      = get_unq_combination(i, j, k, l);
        conjugation_flag f = std::get<1>(unq_perm);
        index_array a      = std::get<0>(unq_perm);
        assert((a[0] < a[3] && a[1] < a[3] && a[2] < a[3]) ||
               (a[0] >= a[3] && a[1] >= a[3] && a[2] >= a[3]));
        if (a[0] < a[3])
            neg_cubes[a[3]][a[0]][a[1]][a[2]] = f ? conj(val) : val;
        else
            pos_cubes[a[3]][a[0] - a[3]][a[1] - a[3]][a[2] - a[3]] = f ? conj(val) : val;
    }

    inline Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> contract(
        const Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1>& vec1,
        const Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1>& vec2,
        const index& i1,
        const index& i2) {
        assert(vec1.size() == size && vec2.size() == size);
        assert(i1 == 0 || i1 == 2);
        assert(i2 == 1 || i2 == 3);

        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic> res(size, size);
        res.setZero();

        switch (i1) {
            case 0:
                switch (i2) {
                    case 1:
                        for (int i = 0; i < size; i++)
                            for (int j = 0; j < size; j++)
                                for (int k = 0; k < size; k++)
                                    for (int l = 0; l < size; l++)
                                        res(i, j) += conj(vec1(k)) * coef(k, l, i, j) * vec2(l);
                        break;
                    case 3:
                        for (int i = 0; i < size; i++)
                            for (int j = 0; j < size; j++)
                                for (int k = 0; k < size; k++)
                                    for (int l = 0; l < size; l++)
                                        res(i, j) += conj(vec1(k)) * coef(k, j, i, l) * vec2(l);
                        break;
                }
                break;
            case 2:
                switch (i2) {
                    case 1:
                        for (int i = 0; i < size; i++)
                            for (int j = 0; j < size; j++)
                                for (int k = 0; k < size; k++)
                                    for (int l = 0; l < size; l++)
                                        res(i, j) += conj(vec1(k)) * coef(i, l, k, j) * vec2(l);
                        break;
                    case 3:
                        for (int i = 0; i < size; i++)
                            for (int j = 0; j < size; j++)
                                for (int k = 0; k < size; k++)
                                    for (int l = 0; l < size; l++)
                                        res(i, j) += conj(vec1(k)) * coef(i, j, k, l) * vec2(l);
                        break;
                }
                break;
        }
        return res;
    }

    void print(std::ostream& os) {
        for (int i = 0; i < size; ++i)
            for (int j = 0; j < size; ++j)
                for (int k = 0; k < size; ++k)
                    for (int l = 0; l < size; ++l) {
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
            return std::move(std::make_pair(index_array({i, j, k, l}), false));
        else if (is_unique_combination(k, l, i, j))
            return std::move(std::make_pair(index_array({k, l, i, j}), false));
        else if (is_unique_combination(j, i, l, k))
            return std::move(std::make_pair(index_array({j, i, l, k}), true));
        else
            return std::move(std::make_pair(index_array({l, k, j, i}), true));
    }
};

typedef Tensor_2E<double> Tensor_2Ecd;

#endif  // TWO_ELECTRON_INTEGRALS_H
