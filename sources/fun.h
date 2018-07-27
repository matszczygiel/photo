#ifndef FUN_H
#define FUN_H

#include "rec.h"
#include <iostream>
#include <cstring>
#include <complex>

#include "constants.h"
#include "harmonics.h"

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SparseCore>


inline int pos_change_gamess( const int L, const int pos ) {
    /* change ordering within the shell to the Gamess convention */
    switch( L ) {
    case 0: // s-type
        return 0;
    case 1: // p-type
        switch( pos ) {
        case 0: return 2;
        case 1: return 1;
        case 2: return 0;
        default:
            std::cout << " Error in pos_change!" << std::endl;
            return EXIT_FAILURE;
        }
    case 2: // d-type
        switch( pos ) {
        case 0: return 2;
        case 1: return 5;
        case 2: return 1;
        case 3: return 4;
        case 4: return 3;
        case 5: return 0;
        default:
            std::cout << " Error in pos_change!" << std::endl;
            return EXIT_FAILURE;
        }
    case 3: // f-type
        switch( pos ) {
        case 0: return 2;
        case 1: return 8;
        case 2: return 6;
        case 3: return 1;
        case 4: return 7;
        case 5: return 9;
        case 6: return 5;
        case 7: return 4;
        case 8: return 3;
        case 9: return 0;
        default:
            std::cout << " Error in pos_change!" << std::endl;
            return EXIT_FAILURE;
        }
    case 4: // g-type
        switch( pos ) {
        case 0 : return 2 ;
        case 1 : return 8 ;
        case 2 : return 11;
        case 3 : return 6 ;
        case 4 : return 1 ;
        case 5 : return 7 ;
        case 6 : return 14;
        case 7 : return 13;
        case 8 : return 5 ;
        case 9 : return 10;
        case 10: return 12;
        case 11: return 9 ;
        case 12: return 4 ;
        case 13: return 3 ;
        case 14: return 0 ;
        default:
            std::cout << " Error in pos_change!" << std::endl;
            return EXIT_FAILURE;
        }
    case 5: // h-type
        switch( pos ) {
        case 0 : return 2 ;
        case 1 : return 8 ;
        case 2 : return 14;
        case 3 : return 12;
        case 4 : return 6 ;
        case 5 : return 1 ;
        case 6 : return 7 ;
        case 7 : return 17;
        case 8 : return 20;
        case 9 : return 16;
        case 10: return 5 ;
        case 11: return 13;
        case 12: return 19;
        case 13: return 18;
        case 14: return 11;
        case 15: return 10;
        case 16: return 15;
        case 17: return 9 ;
        case 18: return 4 ;
        case 19: return 3 ;
        case 20: return 0 ;
        default:
            std::cout << " Error in pos_change!" << std::endl;
            return EXIT_FAILURE;
        }
    case 6: // i-type
        switch( pos ) {
        case 0 : return 2 ;
        case 1 : return 8 ;
        case 2 : return 14;
        case 3 : return 20;
        case 4 : return 12;
        case 5 : return 6 ;
        case 6 : return 1 ;
        case 7 : return 7 ;
        case 8 : return 17;
        case 9 : return 26;
        case 10: return 24;
        case 11: return 16;
        case 12: return 5;
        case 13: return 13;
        case 14: return 25;
        case 15: return 27;
        case 16: return 23;
        case 17: return 11;
        case 18: return 19;
        case 19: return 22;
        case 20: return 21;
        case 21: return 18;
        case 22: return 10;
        case 23: return 15;
        case 24: return 9 ;
        case 25: return 4 ;
        case 26: return 3 ;
        case 27: return 0 ;
        default:
            std::cout << " Error in pos_change!" << std::endl;
            return EXIT_FAILURE;
        }
    default: // higher L not supported in Gamess
        return pos;
    };
};

void GamessOrder(const std::complex<double> * const *shl_crt, std::complex<double> *shl_crt_dum,
                  const int l);

void GamessOrder(double const* const* shl_crt, double *shl_crt_dum,
                  const int l);

void GamessOrderWholeSet(double const* const* const* shl_crt, std::complex<double> *shl_crt_dum,
                  const int l_max);

void Get_1E(const std::string file_name, std::complex<double> *&one_el, const int basis_length_sqrt, bool &error, bool print);

void Get_2E(const std::string file_name, std::complex<double> ****&two_el_matrix, const int basis_length, bool & error, bool print);

double Norm( double a, int k, int l, int m);

void Get_norms(const std::string file_name, double *& norms, bool & error, bool print);

double Energy(const double k, const double ionization_pot);

double Sigma(const double k, const double factor, const double ionization_pot);

double photonEeV(const double k, const double ionization_pot);

double k_length(const double energy_eV, const double ionization_pot);

void Prepare_cont_coef(const std::string norms_file, std::complex<double> *& cont_coef, const int l_max, const int basis_length_k, const int basis_non_k_length,  bool & error );

void Get_HF_data(const std::string file_name, double *& data, const int basis_non_k_length_sqrt, bool & error, bool print);

Eigen::MatrixXd get_CI_coefficients(const std::string file_name, const int basis_non_k_length);

Eigen::VectorXd get_HF_energies(const std::string file_name, const int basis_non_k_length);

#endif
