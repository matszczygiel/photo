#ifndef DATA_HOLDER_H
#define DATA_HOLDER_H

#include <cassert>

#include <eigen3/Eigen/Dense>

#include "two_electron_integrals.h"
#include "disk_reader.h"

class Data_holder
{
public:
  void free();

  void load(const Disk_reader &r);

  inline bool is_full() const { return full; }

  Eigen::MatrixXcd H;
  Eigen::MatrixXcd S;
  Eigen::MatrixXcd Gaugex;
  Eigen::MatrixXcd Gaugey;
  Eigen::MatrixXcd Gaugez;
  Eigen::MatrixXcd HFv;
  Eigen::MatrixXcd CI;

  Tensor_2Ecd Rints;

private:
  bool full = false;
};

#endif