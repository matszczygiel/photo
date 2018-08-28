#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cstring>
#include <string>

#include "functions.h"
#include "constants.h"
#include "harmonics.h"

using namespace std;


double Norm(double a, int k, int l, int m)
{
    double norm_sqrt = pow( (2 * a), -k -l -m -1.5) * Const_arrays::dfact[k] * Const_arrays::dfact[l] * Const_arrays::dfact[m] * pow(0.5, k+l+m) * pow(M_PI, 1.5);
    return sqrt(norm_sqrt);
}

double Energy(const double k, const double ionization_pot) {
    double energy_ion = -2.8616175470;
    double res = energy_ion + k*k*0.5 + ionization_pot;
    return res;
}

double Sigma(const double k, const double factor, const double ionization_pot) {
    double sigma = 2 * 1.24183899890039 * factor * photonEeV(k, ionization_pot) * k;
    return sigma;
}

double photonEeV(const double k, const double ionization_pot) {
    double res = (k*k*0.5 + ionization_pot) * 27.211385;
    return res;
}

double k_length(const double energy_eV, const double ionization_pot)
{
    return sqrt(2.0 * (energy_eV/27.211385 - ionization_pot));
}








