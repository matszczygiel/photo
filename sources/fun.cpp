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

void GamessOrder( complex<double> const* const* shl_crt, complex<double>* shl_crt_dum,
                  const int l) {
    /* transform to the Gamess shell indexing */
    int p,q,p_gam,pos;

    for(p=0; p<l+1; p++) for(q=0; q<l+1-p; q++){
        pos = (l+1)*p - p*(p-1)/2 + q;
        p_gam = pos_change_gamess( l, pos);
        shl_crt_dum[p_gam] = shl_crt[p][q];
    };
};

void GamessOrder( double const* const* shl_crt, double* shl_crt_dum,
                  const int l) {
    /* transform to the Gamess shell indexing */
    int p,q,p_gam,pos;

    for(p=0; p<l+1; p++) for(q=0; q<l+1-p; q++){
        pos = (l+1)*p - p*(p-1)/2 + q;
        p_gam = pos_change_gamess( l, pos);
        shl_crt_dum[p_gam] = shl_crt[p][q];
    };
};

void GamessOrderWholeSet( double const* const* const* shl_crt, complex<double>* shl_crt_dum,
                          const int l_max) {

    int count=0;
    for(int l=0; l<l_max+1; l++)
    {
        double *mem = new double [crt_siz[l]];
        GamessOrder(shl_crt[l], mem,l);

        for(int i=0; i<crt_siz[l]; i++)
        {
            shl_crt_dum[count] = mem[i];
            count++;
        }
        delete[] mem;
    }
    //    cout << count <<endl;

}



void Get_2E(const string file_name, complex<double> ****& two_el_matrix, const int basis_length, bool & error, bool print)
{
    streampos size2E;
    int number_two_el;

    two_el_matrix  = new complex<double>*** [basis_length];
    for(int i =0; i<basis_length; i++) two_el_matrix[i] = new complex<double>** [basis_length];
    for(int i =0; i<basis_length; i++) for(int j =0; j<basis_length; j++)
        two_el_matrix[i][j] = new complex<double>* [basis_length];
    for(int i =0; i<basis_length; i++) for(int j =0; j<basis_length; j++)
        for(int k =0; k<basis_length; k++) two_el_matrix[i][j][k] = new complex<double> [basis_length];

    for(int i =0; i<basis_length; i++) for(int j =0; j<basis_length; j++)
        for(int k =0; k<basis_length; k++)
            memset( two_el_matrix[i][j][k], 0.0, basis_length * sizeof(complex<double>));

    ifstream file2E (file_name, ios::in|ios::binary|ios::ate);
    if (file2E.is_open()){

        size2E = file2E.tellg();
        cout << " The size of 2E file is: " << size2E << endl;
        number_two_el = size2E *sizeof(char)/(2*sizeof(double)+4*sizeof(unsigned short int));

        file2E.seekg (0, ios::beg);

        double re_data, im_data;
        unsigned short int indi, indj, indk, indl;

        for(int i = 0; i < number_two_el; i++){

            file2E.read(reinterpret_cast<char*>(&indi),sizeof(unsigned short int));
            file2E.read(reinterpret_cast<char*>(&indj),sizeof(unsigned short int));
            file2E.read(reinterpret_cast<char*>(&indk),sizeof(unsigned short int));
            file2E.read(reinterpret_cast<char*>(&indl),sizeof(unsigned short int));

            file2E.read(reinterpret_cast<char*>(&re_data),sizeof(double));
            file2E.read(reinterpret_cast<char*>(&im_data),sizeof(double));

            two_el_matrix[indi-1][indj-1][indk-1][indl-1] = two_el_matrix[indk-1][indl-1][indi-1][indj-1] = re_data + im_data * J;
            two_el_matrix[indj-1][indi-1][indl-1][indk-1] = two_el_matrix[indl-1][indk-1][indj-1][indi-1] = re_data - im_data * J;

            if (print)
            {
                cout << "  " << indi << " " << indj << " " << indk << " " << indl <<endl;
                printf("  %e + i%e \n", re_data, im_data);
            }
        }
        file2E.close();
        cout << " The entire 2E file content is in memory. Number of loaded 2E integrals: " << number_two_el<<endl;
    }
    else{
        cout << " Unable to open 2E file"<< endl;
        error = true;
    }
}

void Get_norms(const string file_name, double *& norms, bool & error, bool print){

    streampos size;
    ifstream file (file_name, ios::in|ios::binary|ios::ate);
    if (file.is_open()){

        size = file.tellg();

        cout << " Size of norms file is: " << size << endl;
        int double_size = size*sizeof(char)/sizeof(double);

        norms = new double [double_size];

        file.seekg (0, ios::beg);
        file.read ( reinterpret_cast<char*>(norms), size);
        file.close();
        cout << " The entire norms file content is in memory. Number of loaded norms is: " << double_size <<endl;

        if(print)
        {
            for(int i =0; i<double_size; i++) cout << norms[i]<<endl;
        }
    }
    else {
        cout << " Unable to open norms file"<< endl;
        error = true;
    }
}


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
    //compute Y_lm (k), assume that kx =k ,ky=0, kz=0
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

void Get_HF_data(const string file_name, double *& data, const int basis_non_k_length_sqrt, bool & error, bool print){

    data = new double [basis_non_k_length_sqrt];

    ifstream in;
    string line;

    in.open(file_name);
    if (in.is_open())
    {
        for(int i =0; i<basis_non_k_length_sqrt; i++) in >> data[i];
        if(print){
            cout << " HF data:\n";
            for(int i =0; i<basis_non_k_length_sqrt; i++) cout << data[i] << "\n";
        }
    }
    else error = true;
}


Eigen::MatrixXd get_CI_coefficients(const std::string file_name, const int basis_non_k_length) {
    std::ifstream file(file_name);
    int b_sqrt = basis_non_k_length * basis_non_k_length;
    Eigen::MatrixXd mat ( basis_non_k_length, basis_non_k_length);
    mat.setZero();

    if(file.is_open()) {
        std::stringstream ss;
        std::string line;

        int i, j;
        double v_ij;
        while ( getline ( file, line )) {
            ss.clear();
            ss << line;
            line.clear();
            ss >> i >> j >> v_ij;
            mat(i, j) = v_ij;
        }
    }
    else {
        std::cout << " Unable to open CI file !\n";
    }
    file.close();
    return mat;
}


Eigen::VectorXd get_HF_energies(const std::string file_name, const int basis_non_k_length) {
    std::ifstream file(file_name);
    Eigen::VectorXd vec ( basis_non_k_length);
    vec.setZero();

    if(file.is_open()) {
        std::stringstream ss;
        std::string line;

        for(unsigned i = 0; i < basis_non_k_length; ++i) {

        if (!getline ( file, line )) break;
            ss.clear();
            ss << line;
            line.clear();
            ss >> vec(i);
        }
    }
    else {
        std::cout << " Unable to open HF energies file !\n";
    }
    file.close();
    return vec;
}




