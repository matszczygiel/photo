#include "disk_reader.h"

void Disk_reader::initialize(const Job_control &jc)
{
    job = std::make_shared<Job_control>(jc);
    read_file_basis();
    status = initialized;
}

void Disk_reader::read_file_basis()
{
    assert(status == initialized);

    std::ifstream bfile(job->get_file_basis());
    if (!bfile.is_open())
        throw std::runtime_error("Cannot open the basis file.");

    basis_l = 0;
    basis_lnk = 0;
    bool cont = false;

    for (std::string line; std::getline(bfile, line);)
        if (line == "$BASIS")
            break;

    for (std::string line; std::getline(bfile, line);)
    {
        if (line == "$END")
            break;

        boost::tokenizer<> tok(line);
        *tok.begin() == job->get_continuum_id() ? cont = true : cont = false;

        while (!line.empty())
        {
            std::getline(bfile, line);
            tok.assign(line);
            std::string moment = *tok.begin();
            int end = std::stoi(*(++tok.begin()));

            if (moment == "S")
            {
                basis_l += Const_arrays::crt_siz.at(0);
                if (!cont)
                    basis_lnk += Const_arrays::crt_siz.at(0);
            }
            else if (moment == "P")
            {
                basis_l += Const_arrays::crt_siz.at(1);
                if (!cont)
                    basis_lnk += Const_arrays::crt_siz.at(1);
            }
            else if (moment == "D")
            {
                basis_l += Const_arrays::crt_siz.at(2);
                if (!cont)
                    basis_lnk += Const_arrays::crt_siz.at(2);
            }
            else if (moment == "F")
            {
                basis_l += Const_arrays::crt_siz.at(3);
                if (!cont)
                    basis_lnk += Const_arrays::crt_siz.at(3);
            }
            else if (moment == "G")
            {
                basis_l += Const_arrays::crt_siz.at(4);
                if (!cont)
                    basis_lnk += Const_arrays::crt_siz.at(4);
            }
            else if (moment == "H")
            {
                basis_l += Const_arrays::crt_siz.at(5);
                if (!cont)
                    basis_lnk += Const_arrays::crt_siz.at(5);
            }
            else if (moment == "I")
            {
                basis_l += Const_arrays::crt_siz.at(6);
                if (!cont)
                    basis_lnk += Const_arrays::crt_siz.at(6);
            }
            else if (moment == "K")
            {
                basis_l += Const_arrays::crt_siz.at(7);
                if (!cont)
                    basis_lnk += Const_arrays::crt_siz.at(7);
            }
            else if (moment == "L")
            {
                basis_l += Const_arrays::crt_siz.at(8);
                if (!cont)
                    basis_lnk += Const_arrays::crt_siz.at(8);
            }
            else
                throw std::runtime_error("The basis file contains unrecognized shell.");

            for (int i = 0; i < end; i++)
                std::getline(bfile, line);
        }
    }
    bfile.close();
    status = ready;
}

Eigen::MatrixXcd Disk_reader::load_matrix1E(const int &position) const
{
    assert(status == ready);

    std::ifstream file1E(job->get_file_1E(), std::ios::in | std::ios::binary);
    if (!file1E.is_open())
        throw std::runtime_error("Unable to open 1E file.");

    auto bl_sqrt = get_basis_length_sqrt();

    std::streampos size1E = file1E.tellg();
    int complex_size = size1E * sizeof(char) / sizeof(double);
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

Eigen::MatrixXcd Disk_reader::load_S() const
{
    return load_matrix1E(0);
}

Eigen::MatrixXcd Disk_reader::load_H() const
{
    return load_matrix1E(3);
}

Eigen::MatrixXcd Disk_reader::load_Dipx() const
{
    return load_matrix1E(4);
}

Eigen::MatrixXcd Disk_reader::load_Dipy() const
{
    return load_matrix1E(5);
}

Eigen::MatrixXcd Disk_reader::load_Dipz() const
{
    return load_matrix1E(6);
}

Eigen::MatrixXcd Disk_reader::load_Gradx() const
{
    return load_matrix1E(13);
}

Eigen::MatrixXcd Disk_reader::load_Grady() const
{
    return load_matrix1E(14);
}

Eigen::MatrixXcd Disk_reader::load_Gradz() const
{
    return load_matrix1E(15);
}

Tensor_2Ecd &Disk_reader::load_Rints() const
{
    assert(status == ready);

    std::ifstream file2E(job->get_file_2E(), std::ios::in | std::ios::binary);
    if (!file2E.is_open())
        throw std::runtime_error("Cannot open the 2E file.");

    Tensor_2Ecd ints;
    ints.resize(basis_l);
    ints.set_zero();

    std::streampos file_size = file2E.tellg();
    unsigned long number_two_el = file_size * sizeof(char);
    number_two_el /= 2 * sizeof(double) + 4 * sizeof(unsigned short int);
    file2E.seekg(0, std::ios::beg);

    double re, im;
    unsigned short int indi, indj, indk, indl;

    for (unsigned long i = 0; i < number_two_el; ++i)
    {
        file2E.read(reinterpret_cast<char *>(&indi), sizeof(unsigned short int));
        file2E.read(reinterpret_cast<char *>(&indj), sizeof(unsigned short int));
        file2E.read(reinterpret_cast<char *>(&indk), sizeof(unsigned short int));
        file2E.read(reinterpret_cast<char *>(&indl), sizeof(unsigned short int));

        file2E.read(reinterpret_cast<char *>(&re), sizeof(double));
        file2E.read(reinterpret_cast<char *>(&im), sizeof(double));

        ints.assign(indi, indj, indk, indl, std::complex<double>(re, im));
    }
    file2E.close();
    return ints;
}

Eigen::MatrixXcd Disk_reader::load_Gaugex() const
{
    switch (job->get_gauge())
    {
    case Job_control::gauge_t::dipole:
        return load_Dipx();

    case Job_control::gauge_t::velocity:
        return load_Gradx();

    default:
        return load_Dipx();
    }
}

Eigen::MatrixXcd Disk_reader::load_Gaugey() const
{
    switch (job->get_gauge())
    {
    case Job_control::gauge_t::dipole:
        return load_Dipy();

    case Job_control::gauge_t::velocity:
        return load_Grady();

    default:
        return load_Dipy();
    }
}

Eigen::MatrixXcd Disk_reader::load_Gaugez() const
{
    switch (job->get_gauge())
    {
    case Job_control::gauge_t::dipole:
        return load_Dipz();

    case Job_control::gauge_t::velocity:
        return load_Gradz();

    default:
        return load_Dipz();
    }
}