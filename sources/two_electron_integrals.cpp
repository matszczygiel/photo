#include "two_electron_integrals.h"

using namespace std;
using Cdouble = std::complex<double>;

void read_two_el_from_binary( Two_electron_ints_C& integrals,
                              const string file_name,
                              const unsigned basis_size,
                              bool print ) {
    ifstream file ( file_name, ios::in | ios::binary | ios::ate );
    if ( file.is_open() ) {
        integrals.resize(basis_size, basis_size, basis_size, basis_size);
        integrals.setZero();

        streampos file_size = file.tellg();
        cout << " The size of 2E file is: " << file_size << " bytes. \n";

        unsigned long number_two_el = file_size * sizeof(char);
        number_two_el /= 2 * sizeof(double) + 4 * sizeof(unsigned short int);

        file.seekg (0, ios::beg);

        double              re_data, im_data;
        unsigned short int  indi, indj, indk, indl;
        Cdouble             data, data_bar;

        for(unsigned long i = 0; i < number_two_el; ++i) {
            file.read( reinterpret_cast<char*> (&indi), sizeof(unsigned short int) );
            file.read( reinterpret_cast<char*> (&indj), sizeof(unsigned short int) );
            file.read( reinterpret_cast<char*> (&indk), sizeof(unsigned short int) );
            file.read( reinterpret_cast<char*> (&indl), sizeof(unsigned short int) );

            file.read( reinterpret_cast<char*> (&re_data), sizeof(double) );
            file.read( reinterpret_cast<char*> (&im_data), sizeof(double) );

            data        = Cdouble( re_data, im_data );
            data_bar    = Cdouble( re_data, -im_data );

            if (print) {
                cout << "  " << indi << " " << indj << " " << indk << " " << indl <<"\n";
                cout << data << "\n";
            }

            indi--;
            indj--;
            indk--;
            indl--;

            integrals( indi, indj, indk, indl ) = integrals( indk, indl, indi, indj ) = data;
            integrals( indj, indi, indl, indk ) = integrals( indl, indk, indj, indi ) = data_bar;

        }
        file.close();

        cout << endl;
        cout << " The entire 2E file content is in the memory.\n";
        cout << " Number of loaded 2E integrals: " << number_two_el << "\n\n";
    }
    else {
        cout << " Unable to open 2E file \n";
        return;
    }
}


void read_two_el_from_binary_slice(Two_electron_ints_D &integrals,
                                   const std::string file_name,
                                   const unsigned basis_size,
                                   bool print ) {
    ifstream file ( file_name, ios::in | ios::binary | ios::ate );
    if ( file.is_open() ) {
        integrals.resize(basis_size, basis_size, basis_size, basis_size);
        integrals.setZero();

        streampos file_size = file.tellg();
        cout << " The size of 2E file is: " << file_size << " bytes. \n";

        unsigned long number_two_el = file_size * sizeof(char);
        number_two_el /= 2 * sizeof(double) + 4 * sizeof(unsigned short int);

        file.seekg (0, ios::beg);

        double              re_data, im_data;
        unsigned short int  indi, indj, indk, indl;

        for(unsigned long i = 0; i < number_two_el; ++i) {
            file.read( reinterpret_cast<char*> (&indi), sizeof(unsigned short int) );
            file.read( reinterpret_cast<char*> (&indj), sizeof(unsigned short int) );
            file.read( reinterpret_cast<char*> (&indk), sizeof(unsigned short int) );
            file.read( reinterpret_cast<char*> (&indl), sizeof(unsigned short int) );

            file.read( reinterpret_cast<char*> (&re_data), sizeof(double) );
            file.read( reinterpret_cast<char*> (&im_data), sizeof(double) );

            if(indi <= basis_size && indj <= basis_size && indk <= basis_size && indl <= basis_size) {
                if (print) {
                    cout << "  " << indi << " " << indj << " " << indk << " " << indl <<"\n";
                    cout << re_data << "\n";
                }

                indi--;
                indj--;
                indk--;
                indl--;

                integrals( indi, indj, indk, indl ) = integrals( indk, indl, indi, indj ) = re_data;
                integrals( indj, indi, indl, indk ) = integrals( indl, indk, indj, indi ) = re_data;
                integrals( indj, indi, indk, indl ) = integrals( indk, indl, indj, indi ) = re_data;
                integrals( indi, indj, indl, indk ) = integrals( indl, indk, indi, indj ) = re_data;
            }
        }
        file.close();
    }
    else {
        cout << " Unable to open 2E file \n";
        return;
    }
}


void print_two_el_ints(Two_electron_ints_D &integrals) {
    const auto& dims = integrals.dimensions();

    for(unsigned i = 0; i < dims[0]; ++i) for(unsigned j = 0; j < dims[1]; ++j)
        for(unsigned k = 0; k < dims[2]; ++k) for(unsigned l = 0; l < dims[3]; ++l) {
            cout << "  " << i+1 << " " << j+1 << " " << k+1 << " " << l+1 <<"\n";
            cout << integrals(i, j, k ,l) << "\n";
        }
    cout << "\n";
}


void print_two_el_ints(Two_electron_ints_C &integrals) {
    const auto& dims = integrals.dimensions();

    for(unsigned i = 0; i < dims[0]; ++i) for(unsigned j = 0; j < dims[1]; ++j)
        for(unsigned k = 0; k < dims[2]; ++k) for(unsigned l = 0; l < dims[3]; ++l) {
            cout << "  " << i+1 << " " << j+1 << " " << k+1 << " " << l+1 <<"\n";
            cout << integrals(i, j, k ,l) << "\n";
        }
    cout << "\n";
}


void two_el_int_compare(Two_electron_ints_C &ints_1, Two_electron_ints_D &ints_2, const int b_l) {
    for(unsigned i = 0; i < b_l; ++i) for(unsigned j = 0; j < b_l; ++j)
        for(unsigned k = 0; k < b_l; ++k) for(unsigned l = 0; l < b_l; ++l)
            if( real(ints_1(i, j, k ,l)) != ints_2(i, j, k, l) ) {
                cout << " Diffrence on: " << i << " " << j << " " <<  k << " " << l << "\n";
                cout << ints_2(i, j, k, l) << "\n";
                cout << ints_1(i, j, k ,l) << "\n";
            }
    cout << "\n";
}



