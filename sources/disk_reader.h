#ifndef DISK_READER_H
#define DISK_READER_H

#include <algorithm>
#include <cassert>
#include <complex>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/tokenizer.hpp>
#include <eigen3/Eigen/Dense>

#include "basis.h"
#include "constants.h"
#include "job_control.h"
#include "two_electron_integrals.h"

class Disk_reader {
   public:
    enum Status_t {
        constructed,
        initialized,
        ready
    };

    Disk_reader() = default;

    void initialize(const Job_control &jc);

    int get_basis_length() const { return basis_l; }
    int get_basis_length_sqrt() const { return basis_l * basis_l; }
    int get_basis_length_k() const { return basis_l - basis_lnk; }
    int get_basis_length_k_sqrt() const { return (basis_l - basis_lnk) * (basis_l - basis_lnk); }
    int get_basis_length_nk() const { return basis_lnk; }
    int get_basis_length_nk_sqrt() const { return basis_lnk * basis_lnk; }
    Eigen::Vector3d get_kvec() const { return kvec; }
    double get_kval() const { return kvec.norm(); }
    int get_lmax() const { return lmax; }

    bool is_ready() const {
        return status == ready ? true : false;
    }

    Eigen::MatrixXcd load_S() const;
    Eigen::MatrixXcd load_H() const;
    Eigen::MatrixXcd load_Dipx() const;
    Eigen::MatrixXcd load_Dipy() const;
    Eigen::MatrixXcd load_Dipz() const;
    Eigen::MatrixXcd load_Gradx() const;
    Eigen::MatrixXcd load_Grady() const;
    Eigen::MatrixXcd load_Gradz() const;

    Eigen::MatrixXd load_HFv(const std::string &path) const;
    Eigen::MatrixXd load_CI() const;
    Eigen::VectorXd load_HFe(const std::string &path) const;

    Eigen::VectorXcd load_norms() const;

    Tensor_2Ecd load_Rints() const;

   protected:
    void read_file_basis();
    Eigen::MatrixXcd load_matrix1E_bin(const int &position) const;

   private:
    static constexpr int matrices1E_number = 20;

    std::shared_ptr<Job_control> job = nullptr;
    int basis_l                      = 0;
    int basis_lnk                    = 0;
    int lmax                         = 0;
    Status_t status                  = constructed;
    Eigen::Vector3d kvec             = {0, 0, 0};
};

#endif