#ifndef PHOTO_SCF_H
#define PHOTO_SCF_H

#include "basis.h"
#include "disk_reader.h"
#include "job_control.h"

#include "eigen3/Eigen/Dense"

class PhotoSCF {
   public:
    enum status {
        not_started,
        ready,
        running,
        self_consistent,
        finished,
        iterations_limit
    };

    PhotoSCF(const Job_control &job);

    void run(const Eigen::VectorXcd &vec_ion,
             const Eigen::VectorXcd &vec_cont);
    status one_step();
    void free_ints();

    void set_energy(const double &en) { energy = en; };

    status get_info() const { return info; }
    int get_iter_count() const { return iter_count; }
    int get_max_iter_count() const { return max_iter_count; }
    int get_self_sc_cout() const { return self_sc_cout; }
    double get_treshold() const { return treshold; }
    int get_max_sc_count() const { return max_sc_count; }

    void set_max_iter_count(const int &max_iter_count = 50) { this->max_iter_count = max_iter_count; }
    void set_max_sc_count(const int &max_sc_count = 5) { this->max_sc_count = max_sc_count; }
    void set_ready() { info = ready; }
    void set_treshold(const double &treshold = 0.00001) { this->treshold = treshold; }

   private:
    bool debug = true;

    Disk_reader reader;

    bool evalI, evalC;
    bool force_orth;
    selection_mth_t selection;

    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> es;

    bool starting_vec_loaded = false;
    Eigen::VectorXcd vecI, vecC;
    Eigen::VectorXcd vecIr, vecCr;
    Eigen::VectorXcd vecIrs, vecCrs;

    bool matrices_loaded = false;
    Tensor_2Ecd Rints;
    Eigen::MatrixXcd H, S;
    Eigen::MatrixXd Hnk, Snk;
    Eigen::MatrixXcd Str;
    Eigen::MatrixXcd U;

    bool energy_loaded = false;
    double energy;

    int bl;
    int bkl;
    int bnkl;

   private:
    int self_sc_cout   = 0;
    int max_sc_count   = 5;
    int iter_count     = 0;
    int max_iter_count = 50;
    double treshold    = 0.00001;
    status info        = not_started;
};

#endif