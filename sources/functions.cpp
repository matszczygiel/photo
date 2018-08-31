#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

#include "constants.h"
#include "functions.h"
#include "harmonics.h"

using namespace std;

double Norm(double a, int k, int l, int m) {
    double norm_sqrt = pow((2 * a), -k - l - m - 1.5) * Const_arrays::dfact[k] * Const_arrays::dfact[l] * Const_arrays::dfact[m] * pow(0.5, k + l + m) * pow(M_PI, 1.5);
    return sqrt(norm_sqrt);
}

double Energy(const double k, const double ionization_pot) {
    double energy_ion = -2.8616175470;
    double res        = energy_ion + k * k * 0.5 + ionization_pot;
    return res;
}

double dsigma(const double& photon, const Eigen::Vector3d& polarization, const Eigen::Vector3cd& dipole) {
    double c        = 137.035999139;
    double pre_fct  = 4. * M_PI * M_PI * photon / c;
    double post_fct = std::norm(polarization.dot(dipole));
    return pre_fct * post_fct;
}

double photonEeV(const double k, const double ionization_pot) {
    double res = (k * k * 0.5 + ionization_pot) * 27.211385;
    return res;
}

double k_length(const double energy_eV, const double ionization_pot) {
    return sqrt(2.0 * (energy_eV / 27.211385 - ionization_pot));
}

Eigen::VectorXcd fetch_coulomb_wf(const int& lmax, const Eigen::Vector3d& kvec, const Eigen::VectorXd& norms) {
    std::vector<std::vector<std::vector<double>>> Dfact(lmax + 1);

    for (int l = 0; l <= lmax; l++)
        Dfact[l].resize(l + 1);
    for (int l = 0; l <= lmax; l++)
        for (int p = 0; p <= l; p++)
            Dfact[l][p].resize(l - p + 1);

    using std::pow;
    using std::sqrt;

    for (int l = 0; l <= lmax; l++)
        for (int p = 0; p <= l; p++)
            for (int q = 0; q <= l - p; q++) {
                Dfact[l][p][q] = 0.0;
                for (int m = -l; m <= l; m++)
                    Dfact[l][p][q] += Harmonics::NoNormCalcClmR(l, m, p, q, l - p - q) * Harmonics::Y(l, m, kvec);
                Dfact[l][p][q] *= pow(M_PI, 0.25) * sqrt(Const_arrays::dfact[p] * Const_arrays::dfact[q] * Const_arrays::dfact[l - p - q]) / pow(2.0, 0.25 + l);
            };

    std::vector<std::complex<double>> vec_cont = Gamess::order_set<std::complex<double>, double>(Dfact);

    Eigen::VectorXcd vec = Eigen::Map<Eigen::VectorXcd>(vec_cont.data(), vec_cont.size());
    for (int i = 0; i < vec.size(); i++)
        vec(i) /= norms.tail(vec.size())(i);

    return std::move(vec);
}
