#ifndef DISK_READER_H
#define DISK_READER_H

#include <memory>
#include <fstream>
#include <stdexcept>
#include <string>
#include <cassert>
#include <vector>
#include <complex>
#include <algorithm>

#include <boost/tokenizer.hpp>
#include <eigen3/Eigen/Dense>

#include "job_control.h"
#include "constants.h"
#include "two_electron_integrals.h"

class Disk_reader
{
public:
  enum Status_t
  {
    constructed,
    initialized,
    ready
  };

  Disk_reader() = default;

  void initialize(const Job_control &jc);

  inline int get_basis_length() const { return basis_l; }
  inline int get_basis_length_sqrt() const { return basis_l * basis_l; }
  inline int get_basis_length_k() const { return basis_l - basis_lnk; }
  inline int get_basis_length_k_sqrt() const { return (basis_l - basis_lnk) * (basis_l - basis_lnk); }
  inline int get_basis_length_nk() const { return basis_lnk; }
  inline int get_basis_length_nk_sqrt() const { return basis_lnk * basis_lnk; }

  inline bool is_ready() const { return status == ready ? true : false; }

  Eigen::MatrixXcd load_S() const;
  Eigen::MatrixXcd load_H() const;
  Eigen::MatrixXcd load_Dipx() const;
  Eigen::MatrixXcd load_Dipy() const;
  Eigen::MatrixXcd load_Dipz() const;
  Eigen::MatrixXcd load_Gradx() const;
  Eigen::MatrixXcd load_Grady() const;
  Eigen::MatrixXcd load_Gradz() const;

  Eigen::MatrixXcd load_Gaugex() const;
  Eigen::MatrixXcd load_Gaugey() const;
  Eigen::MatrixXcd load_Gaugez() const;

  Eigen::MatrixXcd load_HF() const;

  Tensor_2Ecd &load_Rints() const;

protected:
  void read_file_basis();
  Eigen::MatrixXcd load_matrix1E(const int &position) const;

private:
  std::shared_ptr<Job_control> job = nullptr;
  int basis_l = 0;
  int basis_lnk = 0;
  Status_t status = constructed;
  const int matrices1E_number = 20;
};

#endif