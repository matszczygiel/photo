#ifndef PHOTO_SCF_H
#define PHOTO_SCF_H

#include "basis.h"
#include "disk_reader.h"
#include "job_control.h"

class PhotoSCF {
   public:

    PhotoSCF(const Job_control& job);

    void run();



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
    double En;
};
}
;

#endif