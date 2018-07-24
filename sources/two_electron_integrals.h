#ifndef TWO_ELECTRON_INTEGRALS_H
#define TWO_ELECTRON_INTEGRALS_H

#include <complex>
#include <iostream>
#include <string>
#include <fstream>

#include <eigen3/unsupported/Eigen/CXX11/Tensor>


using Two_electron_ints_C = Eigen::Tensor<std::complex<double>, 4>;
using Two_electron_ints_D = Eigen::Tensor<double, 4>;

void read_two_el_from_binary( Two_electron_ints_C& integrals,
                              const std::string file_name,
                              const unsigned basis_size,
                              bool print );


void read_two_el_from_binary_slice( Two_electron_ints_D& integrals,
                                    const std::string file_name,
                                    const unsigned basis_size,
                                    bool print );

void print_two_el_ints(Two_electron_ints_D &integrals);


void print_two_el_ints(Two_electron_ints_C &integrals);

void two_el_int_compare(Two_electron_ints_C &ints_1, Two_electron_ints_D &ints_2, const int b_l);


#endif // TWO_ELECTRON_INTEGRALS_H
