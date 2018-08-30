#ifndef PHOTO_SCF_H
#define PHOTO_SCF_H

#include "disk_reader.h"
#include "input_data.h"

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

    enum selection_mth_t {
        by_energy,
        by_norm
    };

    PhotoSCF(const Input_data &data, const std::string & k);

    void run(const Eigen::VectorXcd &vec_ion,
             const Eigen::VectorXcd &vec_cont);
    status one_step();
    void free_ints();


    int get_iter_count() const { return iter_count; }

    void set_max_iter_count(const int &max_iter_count = 50) { this->max_iter_count = max_iter_count; }
    void set_max_sc_count(const int &max_sc_count = 5) { this->max_sc_count = max_sc_count; }
    void set_treshold(const double &treshold = 0.00001) { this->treshold = treshold; }

   private:
    const Disk_reader reader;

    std::string path_2eints;
    std::string path_1eints;

    bool evalI, evalC;
    selection_mth_t selection;

    Eigen::VectorXcd vecI, vecC;
    Eigen::VectorXcd vecIr, vecCr;
    Eigen::VectorXcd vecIrs, vecCrs;

    Tensor_2Ecd Rints;
    Eigen::MatrixXcd H, S;
    Eigen::MatrixXd Hnk, Snk;
    Eigen::MatrixXcd Str;
    Eigen::MatrixXcd U;

    double energy;

    int bl;
    int bkl;
    int bnkl;

   private:
    int self_sc_cout   = 0;
    int max_sc_count   = 5;
    int iter_count     = 0;
    int max_iter_count = 10;
    double treshold    = 0.00001;
    status info        = not_started;
};

#endif
