#ifndef DISK_READER_H
#define DISK_READER_H

#include <memory>
#include <fstream>
#include <stdexcept>
#include <string>
//#include <iomanip>

#include "job_control.h"
#include "constants.h"

class Disk_reader
{
  public:
    Disk_reader() = default;

    void initialize(const Job_control &jc);

    inline int get_basis_length() const { return basis_l; }
    inline int get_basis_length_sqrt() const { return basis_l * basis_l; }
    inline int get_basis_length_k() const { return basis_l - basis_lnk; }
    inline int get_basis_length_k_sqrt() const { return (basis_l - basis_lnk) * (basis_l - basis_lnk); }
    inline int get_basis_length_nk() const { return basis_lnk; }
    inline int get_basis_length_nk_sqrt() const { return basis_lnk * basis_lnk; }

  protected:
    void read_file_basis();

  private:
    std::shared_ptr<Job_control> job;
    int basis_l;
    int basis_lnk;
};

#endif