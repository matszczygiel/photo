#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cstring>
#include <string>

#include "fun.h"
#include "rec.h"
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



void Prepare_cont_coef(const string norms_file, complex<double> *& cont_coef,
                       const int l_max, const int basis_length_k, const int basis_non_k_length,  bool & error )
{
    //store in memory in Gamess order

    double * norms;
    Get_norms(norms_file, norms, error , false);

    double *** D_factor = new double **[l_max+1];
    for(int l =0; l<l_max + 1; l++) D_factor[l] = new double *[l+1];
    for(int l =0; l<l_max + 1; l++) for(int p =0; p<l+1; p++) D_factor[l][p] = new double [l-p+1];

    for(int l =0; l<l_max + 1; l++) for(int p =0; p<l+1; p++) for(int q=0; q<l-p+1; q++)
    {
        
        D_factor[l][p][q] = 0.0;
        for(int m =-l; m <l+1; m++) D_factor[l][p][q]+= Harmonics::NoNormCalcClmR(l,m,p,q,l-p-q) * Harmonics::NoNormCalcClmR(l,m,l,0,0);
        D_factor[l][p][q] *= pow(M_PI, 0.25) * sqrt( Const_arrays::dfact[p] * Const_arrays::dfact[q] * Const_arrays::dfact[l - p - q])  / pow(2.0 , 0.25 + l);
    };

    cont_coef = new complex<double> [basis_length_k];

    GamessOrderWholeSet(D_factor, cont_coef, l_max);

    //    print_matrix_rowmajor(" D factor", basis_length_k, 1, cont_coef);

    for(int l =0; l<l_max + 1; l++) for(int p =0; p<l+1; p++) delete [] D_factor[l][p];
    for(int l =0; l<l_max + 1; l++) delete [] D_factor[l];
    delete [] D_factor;

    for(int i =0; i< basis_length_k; i++) cont_coef[i] /= norms[basis_non_k_length + i];

    delete[] norms;
}






