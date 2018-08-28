#include "disk_reader.h"

#include <algorithm>
#include <stdexcept>
#include <fstream>
#include <complex>


Disk_reader::Disk_reader(const Input_data &data) {
    auto bnkl = std::stoi(data.first("NUMBER_GTO"));
    auto bkl = std::stoi(data.first("NUMBER_PWGTO"));
    
    basis_lnk = bnkl;    
    basis_l = bnkl + bkl;
}


Eigen::MatrixXcd Disk_reader::load_matrix1E_bin(const std::string &path, const int &position) const {
    std::ifstream file1E(path, std::ios::in | std::ios::binary | std::ios::ate);
    if (!file1E.is_open())
        throw std::runtime_error("Unable to open 1E file.");

    auto bl_sqrt = basis_l * basis_l;

    std::streampos size1E = file1E.tellg();
    int complex_size      = size1E * sizeof(char) / sizeof(double);
    complex_size /= 2;

    if (complex_size != bl_sqrt * matrices1E_number)
        throw std::runtime_error("The size of 1E file does not match the basis. Have you used the correct 1E file?");

    std::vector<double> real(bl_sqrt);
    std::vector<double> imag(bl_sqrt);

    int chunk_size = size1E / matrices1E_number;
    file1E.seekg(position * chunk_size, std::ios::beg);
    chunk_size /= 2;
    file1E.read(reinterpret_cast<char *>(real.data()), chunk_size);
    file1E.read(reinterpret_cast<char *>(imag.data()), chunk_size);
    file1E.close();

    Eigen::MatrixXcd ints(basis_l, basis_l);

    std::transform(real.begin(), real.end(), imag.begin(), ints.data(),
                   [](double &dr, double &di) {
                       return std::complex<double>(dr, di);
                   });
    return ints.transpose();
}

Eigen::MatrixXcd Disk_reader::load_S(const std::string &path) const {
    return load_matrix1E_bin(path, 0);
}

Eigen::MatrixXcd Disk_reader::load_H(const std::string &path) const {
    return load_matrix1E_bin(path, 3);
}

Eigen::MatrixXcd Disk_reader::load_Dipx(const std::string &path) const {
    return load_matrix1E_bin(path, 4);
}

Eigen::MatrixXcd Disk_reader::load_Dipy(const std::string &path) const {
    return load_matrix1E_bin(path, 5);
}

Eigen::MatrixXcd Disk_reader::load_Dipz(const std::string &path) const {
    return load_matrix1E_bin(path, 6);
}

Eigen::MatrixXcd Disk_reader::load_Gradx(const std::string &path) const {
    return load_matrix1E_bin(path, 13);
}

Eigen::MatrixXcd Disk_reader::load_Grady(const std::string &path) const {
    return load_matrix1E_bin(path, 14);
}

Eigen::MatrixXcd Disk_reader::load_Gradz(const std::string &path) const {
    return load_matrix1E_bin(path, 15);
}

Tensor_2Ecd Disk_reader::load_Rints(const std::string &path) const {
    std::ifstream file2E(path, std::ios::in | std::ios::binary | std::ios::ate);
    if (!file2E.is_open())
        throw std::runtime_error("Cannot open the 2E file.");

    Tensor_2Ecd ints;
    ints.resize(basis_l);
    ints.zero();

    std::streampos file_size    = file2E.tellg();
    unsigned long number_two_el = file_size * sizeof(char);
    number_two_el /= 2 * sizeof(double) + 4 * sizeof(unsigned short int);
    file2E.seekg(0, std::ios::beg);

    double re, im;
    unsigned short int indi, indj, indk, indl;

    for (unsigned long i = 0; i < number_two_el; ++i) {
        file2E.read(reinterpret_cast<char *>(&indi), sizeof(unsigned short int));
        file2E.read(reinterpret_cast<char *>(&indj), sizeof(unsigned short int));
        file2E.read(reinterpret_cast<char *>(&indk), sizeof(unsigned short int));
        file2E.read(reinterpret_cast<char *>(&indl), sizeof(unsigned short int));

        file2E.read(reinterpret_cast<char *>(&re), sizeof(double));
        file2E.read(reinterpret_cast<char *>(&im), sizeof(double));

        ints.assign(--indi, --indj, --indk, --indl, std::complex<double>(re, im));
    }
    file2E.close();
    return ints;
}

Eigen::MatrixXd Disk_reader::load_HFv(const std::string &path) const {

    Eigen::MatrixXd mat(basis_lnk, basis_lnk);

    std::ifstream file(path);
    if (!file.is_open())
        throw std::runtime_error("Cannot open the HFv file.");

    auto size = basis_lnk * basis_lnk;
    for (int i = 0; i < size; ++i)
        file >> mat.data()[i];

    file.close();
    return mat;
}

Eigen::VectorXd Disk_reader::load_HFe(const std::string &path) const {
    Eigen::VectorXd vec(basis_lnk);

    std::ifstream file(path);
    if (!file.is_open())
        throw std::runtime_error("Cannot open the HFv file.");

    for (int i = 0; i < basis_lnk; ++i)
        file >> vec(i);

    file.close();

    return vec;
}

Eigen::MatrixXd Disk_reader::load_CI(const std::string &path) const {
    Eigen::MatrixXd mat(basis_lnk, basis_lnk);
    mat.setZero();

    std::ifstream file(path);
    if (!file.is_open())
        throw std::runtime_error("Cannot open the CI file.");

    std::stringstream ss;
    std::string line;

    int i, j;
    double v_ij;
    while (std::getline(file, line)) {
        ss.clear();
        ss << line;
        line.clear();
        ss >> i >> j >> v_ij;
        mat(i, j) = v_ij;
    }
    file.close();

    return mat;
}

Eigen::VectorXd Disk_reader::load_norms(const std::string &path) const {
    Eigen::VectorXd vec(basis_l);

    std::ifstream file(path, std::ios::in | std::ios::binary | std::ios::ate);
    if (!file.is_open())
        throw std::runtime_error("Cannot open the norms file.");

    std::streampos size = file.tellg();
    int double_size     = size * sizeof(char) / sizeof(double);
    if (double_size != basis_l)
        throw std::runtime_error("Size of the norms file is not consistent with basis. Have you used the correct norms file?");

    file.seekg(0, std::ios::beg);
    file.read(reinterpret_cast<char *>(vec.data()), size);
    file.close();

    return vec;
}
