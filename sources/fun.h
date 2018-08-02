#ifndef FUN_H
#define FUN_H

#include <complex>

#include "constants.h"
#include "harmonics.h"
#include "gamess.h"

#include <eigen3/Eigen/Dense>

template <class type>
inline type scalar_prod(const Eigen::Matrix<type, Eigen::Dynamic, 1> &vec1,
                        const Eigen::Matrix<type, Eigen::Dynamic, Eigen::Dynamic> &kernel,
                        const Eigen::Matrix<type, Eigen::Dynamic, 1> &vec2)
{
    return vec1.dot(kernel * vec2));
}

double Norm(double a, int k, int l, int m);

double Energy(const double k, const double ionization_pot);

double Sigma(const double k, const double factor, const double ionization_pot);

double photonEeV(const double k, const double ionization_pot);

double k_length(const double energy_eV, const double ionization_pot);

void Prepare_cont_coef(const std::string norms_file, std::complex<double> *&cont_coef, const int l_max, const int basis_length_k, const int basis_non_k_length, bool &error);

#endif
