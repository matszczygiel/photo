#ifndef JOB_CONTROL_H
#define JOB_CONTROL_H

#include <boost/algorithm/string.hpp>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "input_data.h"

enum selection_mth_t {
    energy,
    norm
};

enum gauge_t {
    velocity,
    dipole
};

class Job_control {
   public:
    Job_control();

    void read(Input_data &in);

    inline std::string get_file_basis() const { return path_in + file_basis; }
    inline std::string get_file_1E() const { return path_in + file_1E; }
    inline std::string get_file_2E() const { return path_in + file_2E; }
    inline std::string get_file_norm() const { return path_in + file_norm; }
    inline std::string get_file_HFv_I() const { return path_in + file_HFv_i; }
    inline std::string get_file_HFv_I() const { return path_in + file_HFv_i; }
    inline std::string get_file_HFe_F() const { return path_in + file_HFe_f; }
    inline std::string get_file_CI_I() const { return path_in + file_CI_i; }
    inline std::string get_continuum_id() const { return continuum_id; }
    inline bool is_computeI() const { return compute_bound; }
    inline bool is_computeC() const { return compute_cont; }
    inline double get_energy_i() const { return energy_i; }
    inline double get_photon() const { return photonEeV; };
    inline double get_nucl_rep_en() const { return 1. / R; };

    friend class PhotoSCF;

   private:
    bool compute_bound;
    bool compute_cont;

    double energy_i;
    double photonEv;
    double R;

    std::string path_in;
    std::string file_basis;
    std::string file_1E;
    std::string file_2E;
    std::string file_norm;

    std::string file_HFv_i;
    std::string file_CI_i;

    std::string file_HFv_f;
    std::string file_HFe_f;

    std::string continuum_id;

    selection_mth_t selection_m;
    gauge_t gauge;
    bool force_orth;

    bool write;
    std::string path_out;
    std::string file_out;
};

#endif  // JOB_CONTROL_H
