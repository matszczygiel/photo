#include "disk_reader.h"

void Disk_reader::initialize(const Job_control &jc) {
    job    = std::make_shared<Job_control>(jc);
    status = initialized;
    read_file_basis();
}

void Disk_reader::read_file_basis() {
    assert(status == initialized);

    std::ifstream bfile(job->get_file_basis());
    if (!bfile.is_open())
        throw std::runtime_error("Cannot open the basis file.");

    for (std::string line; std::getline(bfile, line);)
        if (line == "$BASIS")
            break;

    std::vector<Basis> basis;
    Basis b;
    int basis_cont_l = 0;
    basis_l          = 0;
    lmax             = 0;
    while (b.read(bfile)) {
        basis.push_back(b);
        basis_l += b.functions_number();
        int max_shl = shell2int(b.get_max_shell());
        if (max_shl > lmax)
            lmax = max_shl;
        if (b.get_label() == "CONT") {
            kvec(0) = b.get_kvec()[0];
            kvec(1) = b.get_kvec()[1];
            kvec(2) = b.get_kvec()[2];
            basis_cont_l += b.functions_number();
        }
    }

    bfile.close();

    basis_lnk = basis_l - basis_cont_l;
    status    = ready;
}

Eigen::MatrixXcd Disk_reader::load_matrix1E_bin(const int &position) const {
    assert(status == ready);

    std::ifstream file1E(job->get_file_1E(), std::ios::in | std::ios::binary | std::ios::ate);
    if (!file1E.is_open())
        throw std::runtime_error("Unable to open 1E file.");

    auto bl_sqrt = get_basis_length_sqrt();

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

Eigen::MatrixXcd Disk_reader::load_S() const {
    return load_matrix1E_bin(0);
}

Eigen::MatrixXcd Disk_reader::load_H() const {
    return load_matrix1E_bin(3);
}

Eigen::MatrixXcd Disk_reader::load_Dipx() const {
    return load_matrix1E_bin(4);
}

Eigen::MatrixXcd Disk_reader::load_Dipy() const {
    return load_matrix1E_bin(5);
}

Eigen::MatrixXcd Disk_reader::load_Dipz() const {
    return load_matrix1E_bin(6);
}

Eigen::MatrixXcd Disk_reader::load_Gradx() const {
    return load_matrix1E_bin(13);
}

Eigen::MatrixXcd Disk_reader::load_Grady() const {
    return load_matrix1E_bin(14);
}

Eigen::MatrixXcd Disk_reader::load_Gradz() const {
    return load_matrix1E_bin(15);
}

Tensor_2Ecd Disk_reader::load_Rints() const {
    assert(status == ready);

    std::ifstream file2E(job->get_file_2E(), std::ios::in | std::ios::binary | std::ios::ate);
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

        --indi;
        --indj;
        --indk;
        --indl;
        ints.assign(indi, indj, indk, indl, std::complex<double>(re, im));
    }
    file2E.close();
    return ints;
}



Eigen::MatrixXd Disk_reader::load_HFv() const {
    assert(status == ready);

    Eigen::MatrixXd mat(basis_lnk, basis_lnk);

    std::ifstream file(job->get_file_HFv());
    if (!file.is_open())
        throw std::runtime_error("Cannot open the HFv file.");

    for (int i = 0; i < get_basis_length_nk_sqrt(); ++i)
        file >> mat.data()[i];

    file.close();

    return mat;
}

Eigen::MatrixXd Disk_reader::load_CI() const {
    assert(status == ready);

    Eigen::MatrixXd mat(basis_lnk, basis_lnk);
    mat.setZero();

    std::ifstream file(job->get_file_CI());
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

Eigen::VectorXcd Disk_reader::load_norms() const {
    assert(status == ready);

    Eigen::VectorXcd vec(basis_l);

    std::ifstream file(job->get_file_norm(), std::ios::in | std::ios::binary | std::ios::ate);
    if (!file.is_open())
        throw std::runtime_error("Cannot open the norms file.");

    std::streampos size = file.tellg();
    int double_size     = size * sizeof(char) / sizeof(double);
    std::cout << size << std::endl;
    std::cout << double_size << std::endl;
    if (double_size != basis_l)
        throw std::runtime_error("Size of the norms file is not consistent with basis. Have you used the correct norms file?");

    file.seekg(0, std::ios::beg);
    file.read(reinterpret_cast<char *>(vec.data()), size);
    file.close();

    return vec;
}
