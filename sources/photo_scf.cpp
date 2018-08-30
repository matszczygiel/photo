#include <cassert>
#include <future>
#include <iostream>
#include <stdexcept>
#include <thread>

#include "functions.h"
#include "photo_scf.h"

using namespace std;
using namespace Eigen;

PhotoSCF::PhotoSCF(const Input_data &data, const std::string &k)
    : reader(data) {
    char token;
    std::string arg = "COMPUTE_BOUND";

    token = std::tolower(*data.first(arg).begin());
    if (token == 'y')
        evalI = true;
    else if (token == 'n')
        evalI = false;
    else
        throw std::runtime_error("Invalid argument for" + arg + ".");

    arg = "COMPUTE_CONTINUUM";

    token = std::tolower(*data.first(arg).begin());
    if (token == 'y')
        evalC = true;
    else if (token == 'n')
        evalC = false;
    else
        throw std::runtime_error("Invalid argument for" + arg + ".");

    arg = "SELECTION_METHOD";

    if (data.first(arg) == "energy")
        selection = by_energy;
    else if (data.first(arg) == "norm")
        selection = by_norm;
    else
        throw std::runtime_error("Invalid argument for" + arg + ".");

    path_2eints = data.first("PATH_IN") + data.first("FILES_2E") + k + data.second("FILES_2E");
    path_1eints = data.first("PATH_IN") + data.first("FILES_1E") + k + data.second("FILES_1E");

    auto en_i = std::stof(data.first("ENERGY_I_STATE"));
    en_i -= 1. / std::stof(data.first("R"));

    double photon = std::stof(data.first("PHOTON_EN")) / 27.211385;

    energy = en_i + photon;

    bnkl = std::stoi(data.first("NUMBER_GTO"));
    bkl  = std::stoi(data.first("NUMBER_PWGTO"));
    bl   = bnkl + bkl;
}

void PhotoSCF::run(const Eigen::VectorXcd &vec_ion,
                   const Eigen::VectorXcd &vec_cont) {
    vecI = VectorXcd::Zero(bl);
    vecC = VectorXcd::Zero(bl);

    vecI.head(bnkl) = vec_ion;
    vecC.tail(bkl)  = vec_cont;

    if (evalC || evalI)
        Rints = reader.load_Rints(path_2eints);

    H = reader.load_H(path_1eints);
    S = reader.load_S(path_1eints);

    Hnk = H.topLeftCorner(bnkl, bnkl).real();
    Snk = S.topLeftCorner(bnkl, bnkl).real();

    U                         = MatrixXcd::Zero(bl, bnkl + 1);
    U.col(0)                  = vecC;
    U.block(0, 1, bnkl, bnkl) = MatrixXd::Identity(bnkl, bnkl);

    vecCrs    = VectorXcd::Zero(bnkl + 1);
    vecCrs(0) = 1;

    Str = U.adjoint() * S * U;

    vecCr = vecCrs;

    vecIrs = vec_ion;
    vecIr  = vecIrs;

    info = ready;

    if (evalI || evalC) {
        std::cout << "\n"
                  << " Starting iteration.\n";

        info = running;
        while (info == running) {
            iter_count++;
            std::cout << "\n";
            std::cout << "========= Iteration " << iter_count << " ============="
                      << "\n\n";

            info = one_step();

            if (info == self_consistent) {
                self_sc_cout++;
                info = running;
            } else
                self_sc_cout = 0;

            if (self_sc_cout == max_sc_count) {
                info = finished;
                break;
            }
            if (iter_count == max_iter_count) {
                info = iterations_limit;
                std::cout << " Iterations limit reached. \n";
                break;
            }
            if (self_sc_cout != 0)
                std::cout << "\n"
                          << " Self consistency counter:" << self_sc_cout << "\n\n";
        }
        std::cout << "\n";
        std::cout << " ======= End of iteration "
                  << " ==========="
                  << "\n\n";
    }
}

void PhotoSCF::free_ints() {
    Rints.resize(0);
    H.resize(0, 0);
    S.resize(0, 0);
    Hnk.resize(0, 0);
    Snk.resize(0, 0);
    Str.resize(0, 0);
    U.resize(0, 0);
}

PhotoSCF::status PhotoSCF::one_step() {
    // R matrices preparation
    auto kRk_future = std::async(std::launch::async, &Tensor_2Ecd::contract, &Rints, vecC, vecC, 0, 3);
    auto kkR_future = std::async(std::launch::async, &Tensor_2Ecd::contract, &Rints, vecC, vecC, 0, 1);
    auto ppR_future = std::async(std::launch::async, &Tensor_2Ecd::contract, &Rints, vecI, vecI, 0, 1);
    auto pRp_future = std::async(std::launch::async, &Tensor_2Ecd::contract, &Rints, vecI, vecI, 0, 3);

    MatrixXcd kRk = kRk_future.get().topLeftCorner(bnkl, bnkl);
    MatrixXcd kkR = kkR_future.get().topLeftCorner(bnkl, bnkl);
    MatrixXcd ppR = ppR_future.get();
    MatrixXcd pRp = pRp_future.get();

    /*
    auto kRk = Rints.contract(vecC, vecC, 0, 3);
    kRk      = kRk.topLeftCorner(bnkl, bnkl);
    auto kkR = Rints.contract(vecC, vecC, 0, 1);
    kkR      = kkR.topLeftCorner(bnkl, bnkl);
    auto ppR = Rints.contract(vecI, vecI, 0, 1);
    auto pRp = Rints.contract(vecI, vecI, 0, 3);
*/
    // H matrix elements
    auto Hpp = real(scalar_prod(vecI, H, vecI));
    auto Hkk = real(scalar_prod(vecC, H, vecC));
    auto Hpk = scalar_prod(vecI, H, vecC);

    // Print H matrix elements

    cout << " Energy of the system is: " << energy << "\n";
    cout << "\n"
         << " H matrix elements: "
         << "\n";
    cout << " H_pp = " << Hpp << "\n"
         << " H_kk = " << Hkk << "\n"
         << " H_pk = " << Hpk << "\n\n";

    auto normC = real(scalar_prod(vecC, S, vecC));
    auto normI = real(scalar_prod(vecI, S, vecI));

    cout << "\n"
         << " Norm squared of the continuum state is: " << normC << "\n";
    cout << "\n"
         << " Norm squared of the bound state is: " << normI << "\n";

    /// Compose the set of equation of form Ax=ESx

    //A matrices prep

    MatrixXcd AmatI =
        (normC * Hnk) +
        (H.topRows(bnkl) * vecC) * (vecC.adjoint() * S.leftCols(bnkl)) +
        (S.topRows(bnkl) * vecC) * (vecC.adjoint() * H.leftCols(bnkl)) +
        (Hkk * Snk) + kRk + kkR;

    MatrixXcd AmatC =
        U.adjoint() *
        ((normI * H) + (H * vecI) * (vecI.adjoint() * S) + (S * vecI) * (vecI.adjoint() * H) + Hpp * S + pRp + ppR) *
        U;

    // assert(AmatC.rows() == (bnkl + 1));

    cout << " A matrices preparation done.\n\n";

    //S matrices prep

    MatrixXcd SmatI =
        (normC * Snk) + (S.topRows(bnkl) * vecC) * (vecC.adjoint() * S.leftCols(bnkl));

    MatrixXcd SmatC =
        U.adjoint() *
        ((normI * S) + (S * vecI) * (vecI.adjoint() * S)) *
        U;

    cout << " S matrices preparation done. \n\n";
    cout << " Solving eigenvalue problem. \n\n";

    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> es;

    VectorXcd vecIr_dum, vecCr_dum;

    if (evalC) {
        es.compute(AmatC, SmatC);
        VectorXd temp;
        int itC;
        int sizeC = AmatC.cols();

        cout << " C eigenvalues:\n";
        cout << es.eigenvalues() << "\n";
        /*
        for (int i = 0; i < bl; ++i) {
            cout << "**********\n";
            cout << es.eigenvalues()[i];
            cout << "\n";
            cout << es.eigenvectors().col(i);
            cout << "\n";
        }
*/
        switch (selection) {
            case selection_mth_t::by_energy:
                temp = es.eigenvalues() - energy * VectorXd::Ones(es.eigenvalues().size());
                temp = temp.cwiseAbs2();
                temp.minCoeff(&itC);
                cout << " Iterator to the matching energy of C: " << itC
                     << " \n and it's enegry: " << es.eigenvalues()[itC] << "\n\n";

                break;

            case selection_mth_t::by_norm:
                temp.resize(sizeC);
                for (int i = 0; i < sizeC; ++i) {
                    vecCr_dum = es.eigenvectors().col(i) / es.eigenvectors()(0, i);
                    temp(i)   = real(
                        vecCr_dum.tail(sizeC - 1).dot(Str.bottomRightCorner(sizeC - 1, sizeC - 1) * vecCr_dum.tail(sizeC - 1)));
                }
                temp.minCoeff(&itC);
                cout << " Iterator to the matching norm of C: " << itC;
                cout << "\n and it's energy: " << es.eigenvalues()[itC] << "\n\n";
                break;
        }
        vecCr_dum = es.eigenvectors().col(itC) / es.eigenvectors()(0, itC);
    }

    if (evalI) {
        es.compute(AmatI, SmatI);
        VectorXd temp;
        int itI;

        cout << " I eigenvalues:\n";
        cout << es.eigenvalues();
        cout << "\n";

        switch (selection) {
            case selection_mth_t::by_energy:
                temp = es.eigenvalues() - energy * VectorXd::Ones(es.eigenvalues().size());
                temp = temp.cwiseAbs2();
                temp.minCoeff(&itI);
                cout << " Iterator to the matching energy of I: " << itI
                     << " \n and it's enegry: " << es.eigenvalues()[itI] << "\n\n";

                break;

            case selection_mth_t::by_norm:
                temp.resize(bnkl);
                for (int i = 0; i < bnkl; ++i) {
                    vecIr_dum = es.eigenvectors().col(i);
                    vecIr_dum /= sqrt(real(vecIr_dum.dot(Snk * vecIr_dum)));
                    temp(i) = std::norm(vecIrs.dot(Snk * vecIr_dum));
                }
                temp.maxCoeff(&itI);
                cout << " Iterator to the matching norm of I: " << itI;
                cout << "\n and it's energy: " << es.eigenvalues()[itI] << "\n\n";

                break;
        }

        vecIr_dum = es.eigenvectors().col(itI);
        vecIr_dum /= sqrt(real(vecIr_dum.dot(Snk * vecIr_dum)));
    }

    //Check for self consistency

    double conv_paramC = 0.0;
    double conv_paramI = 0.0;

    if (evalC) {
        conv_paramC = (vecCr_dum.cwiseAbs() - vecCr.cwiseAbs()).norm();
        cout << "\n Norm of difference for eq1: " << conv_paramC << "\n";
        vecCr = vecCr_dum;
    }
    if (evalI) {
        conv_paramI = (vecIr_dum.cwiseAbs() - vecIr.cwiseAbs()).norm();
        cout << "\n Norm of difference for eq2: " << conv_paramI << "\n";
        vecIr = vecIr_dum;
    }

    vecC = U * vecCr;
    vecI << vecIr, VectorXcd::Zero(bl - bnkl);

    if (!(evalI && evalC))
        return finished;

    if (conv_paramC < treshold && conv_paramI < treshold)
        return self_consistent;
    else
        return running;
}