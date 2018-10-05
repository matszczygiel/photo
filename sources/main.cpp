/* A photoionization cros section evaluation program. The calculation is made by fixing asymptotic of the continuum wave function.
 * The final state is build as antisymetrized product of bound and continuum state. The mean field approximation is used to correct input orbitals.
 * The input is made by specyfing the asymptotic of continuum wave function.
 * Orthogonality of the continuum wave function to the ground state orbital is optional.
 * Two selection methods of orbitals are avaliable.
 *
 * The program uses the Eigen linear algebra libary.
 *
 * M. S. Szczygiel 2018
 */

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#include "constants.h"
#include "disk_reader.h"
#include "functions.h"
#include "gamess.h"
#include "harmonics.h"
#include "photo_scf.h"
#include "two_electron_integrals.h"

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[]) {
    auto start = chrono::system_clock::now();

    if (!(argc == 3 || argc == 2)) {
        cout << " Proper usage: ./photo <input name> <settings>\n";
        return EXIT_SUCCESS;
    }

    string input   = argv[1];
    string setting = "n";
    if (argc == 3)
        setting = argv[2];

    ifstream ifile(input);
    const Input_data data(ifile);
    ifile.close();

    const Disk_reader reader(data);

    string orbs_f = data.first("PATH_IN") + data.first("FILE_HF_F_VEC");
    string enrg_f = data.first("PATH_IN") + data.first("FILE_HF_F_EN");

    cout.precision(5);

    ////////////////////////////

    auto orbitals_ion = reader.load_HFv(orbs_f);
    auto energies_ion = reader.load_HFe(enrg_f);

    auto en_i = std::stof(data.first("ENERGY_I_STATE"));
    en_i -= 1. / std::stof(data.first("R"));

    double photon = std::stof(data.first("PHOTON_EN")) / 27.211385;

    double en_final = en_i + photon;

    vector<int> indices;
    vector<double> kvals;
    for (int i = 0; i < energies_ion.size(); ++i) {
        if (energies_ion(i) < en_final) {
            indices.push_back(i);
            kvals.push_back(sqrt(2. * (en_final - energies_ion(i))));
        } else
            break;
    }

    cout << "K:  ";
    for (auto &x : kvals) {
        std::stringstream stream;
        stream << std::fixed << std::setprecision(3) << x;
        cout << stream.str() << "  ";
    }
    cout << "\n";

    if (setting == "-dump")
        return 0;

    ///////////////////////////////

    std::stringstream stream;
    stream << std::fixed << std::setprecision(3) << kvals[0];
    std::string k_str = stream.str();

    auto lmax = std::stoi(data.first("MAX_L"));

    double ptheta = std::stof(data.first("POL_THETA"));
    double pphi   = std::stof(data.first("POL_PHI"));
    /////////////
    vector<string> norm_files, ints_files, Rints_files;
    vector<double> theta;

    int job_size = data.size("FILES_NORM") - 1;

    for (int i = 0; i < job_size; ++i) {
        norm_files.push_back(data.first("PATH_IN") + data.first("FILES_NORM") + k_str + data("FILES_NORM", i + 1));
        ints_files.push_back(data.first("PATH_IN") + data.first("FILES_1E") + k_str + data("FILES_1E", i + 1));
        Rints_files.push_back(data.first("PATH_IN") + data.first("FILES_2E") + k_str + data("FILES_2E", i + 1));

        theta.push_back(std::stof(data("K_THETA", i)));
    }

    double phi = std::stof(data.first("K_PHI"));
    vector<vector<double>> sigma(job_size, vector<double>(kvals.size(), 0.));

    for (int k = 0; k < static_cast<int>(kvals.size()); ++k) {
        for (int i = 0; i < job_size; ++i) {
            Vector3d kvec;
            kvec(0) = kvals[k] * sin(theta.at(i)) * cos(phi);
            kvec(1) = kvals[k] * sin(theta.at(i)) * sin(phi);
            kvec(2) = kvals[k] * cos(theta.at(i));

            auto norms = reader.load_norms(norm_files.at(i));

            VectorXcd cont_vec = fetch_coulomb_wf(lmax, kvec, norms);

            //////////////////////

            PhotoSCF sys(data, k_str);
            sys.run(orbitals_ion.col(indices[k]), cont_vec);
            sys.free_ints();

            auto vecI = sys.getI();
            auto vecC = sys.getC();

            auto bnkl = std::stoi(data.first("NUMBER_GTO"));

            //////////////////////////////////

            auto veci = vecI.head(bnkl);

            auto Dx = reader.load_Dipx(ints_files.at(i));
            auto Dy = reader.load_Dipy(ints_files.at(i));
            auto Dz = reader.load_Dipz(ints_files.at(i));

            auto S   = reader.load_S(ints_files.at(i));
            auto Snk = S.topLeftCorner(bnkl, bnkl);

            auto Dxnk = Dx.topLeftCorner(bnkl, bnkl);
            auto Dynk = Dy.topLeftCorner(bnkl, bnkl);
            auto Dznk = Dz.topLeftCorner(bnkl, bnkl);

            string orbs_i_s = data.first("PATH_IN") + data.first("FILE_HF_I_VEC");
            auto orbs_i     = reader.load_HFv(orbs_i_s);
            VectorXd grHF   = orbs_i.col(0);

            complex<double> S_kI  = vecC.dot(S.leftCols(bnkl) * grHF);
            complex<double> S_pI  = veci.dot(Snk * grHF);
            complex<double> Dx_pI = veci.dot(Dxnk * grHF);
            complex<double> Dy_pI = veci.dot(Dynk * grHF);
            complex<double> Dz_pI = veci.dot(Dznk * grHF);
            complex<double> Dx_kI = vecC.dot(Dx.leftCols(bnkl) * grHF);
            complex<double> Dy_kI = vecC.dot(Dy.leftCols(bnkl) * grHF);
            complex<double> Dz_kI = vecC.dot(Dz.leftCols(bnkl) * grHF);

            Vector3cd T;
            T(0) = S_kI * Dx_pI + S_pI * Dx_kI;
            T(1) = S_kI * Dy_pI + S_pI * Dy_kI;
            T(2) = S_kI * Dz_pI + S_pI * Dz_kI;
            T *= sqrt(2.);

            //normalize to energy scale
            T *= sqrt(kvals[k]);

            cout << " Dipole moment: \n";
            cout << T << "\n\n";

            Vector3d j;
            j(0) = sin(ptheta) * cos(pphi);
            j(1) = sin(ptheta) * sin(pphi);
            j(2) = cos(ptheta);

 /*           double r_pI_sqrt, r_kI_sqrt;

            r_pI_sqrt = norm(Dx_pI) + norm(Dy_pI) + norm(Dz_pI);
            r_kI_sqrt = norm(Dx_kI) + norm(Dy_kI) + norm(Dz_kI);
*/
            auto S_pk = vecI.dot(S * vecC);
/*
            cout << "\n"
                 << " S_pk : " << S_pk << "\n\n";

            cout << " Dx_kI:     " << Dx_kI << "\n";
            cout << " Dy_kI:     " << Dy_kI << "\n";
            cout << " Dz_kI:     " << Dz_kI << "\n";
            cout << " \n\n";

            cout << " S_kI:      " << S_kI << "\n";
            cout << " S_pI:      " << S_pI << "\n";
            cout << " r_pI_sqrt: " << r_pI_sqrt << "\n";
            cout << " r_kI_sqrt: " << r_kI_sqrt << "\n";
            cout << " \n\n";

    auto Rints = reader.load_Rints(Rints_file);

    MatrixXcd pp_R  = Rints.contract(vecI, vecI, 0, 1);
    MatrixXcd p_R_p = Rints.contract(vecI, vecI, 0, 3);

    Rints.resize(0);

    pk_R_pk     = real(vecC.dot(pp_R * vecC));
    pk_R_kp     = real(vecC.dot(p_R_p * vecC));
    */
            double norm_k_sqrt = real(vecC.dot(S * vecC));

            double norm_psi = 2 * (norm_k_sqrt + norm(S_pk));
            cout << " Norm of the whole state: " << norm_psi << "\n\n";

            //    double E_rep_corr = (pk_R_kp + pk_R_pk) / (0.5 * norm_psi);

            //  cout << " Electron repulsion energy: " << E_rep_corr << "\n\n";



            // keep this for He
            sigma.at(i).at(k) += sigma_tot_spherical_symetry(photon, T);

            //sigma.at(i).at(k) += dsigma(photon, j, T);

            cout << " To state:            " << fixed << indices[k] << "\n";
            cout << " Photon energy [eV]:  " << fixed << setprecision(3) << data.first("PHOTON_EN") << "\n";
            cout << " Cross section :      " << fixed << setprecision(4) << sigma.at(i).at(k) << "\n";
            cout << " \n\n\n";
        }
    }

    //write results

    bool write;
    char token;
    std::string arg = "WRITE";

    token = std::tolower(*(data.first(arg).begin()));
    if (token == 'y')
        write = true;
    else if (token == 'n')
        write = false;
    else
        throw std::runtime_error("Invalid argument for" + arg + ".");

    if (write) {
        string res_path = data.first("PATH_OUT") + data.first("FILE_OUT");

        std::ofstream outfile(res_path, std::ios_base::app);
        outfile << std::fixed;
        outfile << "****** " << data.first("NAME") << " ******\n";
        outfile << "Photon [eV]         " << data.first("PHOTON_EN") << "\n";
        outfile << "k ( phi)            " << std::setprecision(3) << phi << "\n";
        outfile << "j (theta, phi)      " << std::setprecision(3) << ptheta << "\t" << std::setprecision(3) << pphi << "\n";
        outfile << "========\n";
        outfile << "k(theta)";
        for (const auto &x : indices)
            outfile << "\t" << x << "\t";
        outfile << "\ttot sigma\n";
        for (int i = 0; i < job_size; ++i) {
            double sig_tot = 0;
            outfile << std::setprecision(3) << theta.at(i) << "\t\t";
            for (const auto &x : sigma.at(i)) {
                outfile << std::setprecision(5) << x << "\t\t";
                sig_tot += x;
            }
            outfile << std::setprecision(5) << sig_tot << "\n";
        }
        outfile << "\n";
        outfile.close();
    }

    // if ( jc.getIfWrite() ) jc.writeRes(sig);

    // TEST OF CI
    // H transformation
    /*	MatrixXd H = HF.transpose() * Hnk.real() * HF;
        MatrixXd CI = get_CI_coefficients(jc.getFileCI(), bnkl);
        double En_CI = 2. * (CI * H * CI).trace();

        Two_electron_ints_D two_el_real, temp_two;
        read_two_el_from_binary_slice(two_el_real, jc.getFile2E(), bnkl, false);

        TensorMap<Tensor<double, 2>> CI_tens(CI.data(), bnkl, bnkl);
        TensorMap<Tensor<double, 2>> HF_tens(HF.data(), bnkl, bnkl);

        // Two_el transformation
        Eigen::array<Eigen::IndexPair<int>, 1> dims_0, dims_1, dims_2, dims_3;
        dims_0 = {IndexPair<int>(0, 0)};
        dims_1 = {IndexPair<int>(0, 1)};
        dims_2 = {IndexPair<int>(0, 2)};
        dims_3 = {IndexPair<int>(0, 3)};

        temp_two = HF_tens.contract(two_el_real, dims_0);
        two_el_real = temp_two;
        temp_two = HF_tens.contract(two_el_real, dims_1);
        two_el_real = temp_two;
        temp_two = HF_tens.contract(two_el_real, dims_2);
        two_el_real = temp_two;
        temp_two = HF_tens.contract(two_el_real, dims_3);
        two_el_real = temp_two;

        Eigen::array<Eigen::IndexPair<int>, 2> dims_pair_0, dims_pair_1;
        dims_pair_0 = {IndexPair<int>(0, 0), IndexPair<int>(1, 2)};
        dims_pair_1 = {IndexPair<int>(0, 0), IndexPair<int>(1, 1)};
        Tensor<double, 2> temp = CI_tens.contract(two_el_real, dims_pair_0);
        Tensor<double, 0> temp_2 = CI_tens.contract(temp, dims_pair_1);

        En_CI += temp_2(0);

        auto HF_energies = get_HF_energies(jc.getFileHFEnergy(), bnkl);

        cout << " CI energy: " << En_CI << "\n";

        // Dipole moment
        MatrixXcd T_ij_x, T_ij_y, T_ij_z;
        T_ij_x = (vecC.adjoint() * Dx.leftCols(bnkl)).transpose() * (veci.adjoint() * Snk);
        T_ij_x += (veci.adjoint() * Dxnk).transpose() * (vecC.adjoint() * S.leftCols(bnkl));
        T_ij_x = HF.transpose() * T_ij_x * HF;
        T_ij_x += T_ij_x.transpose().eval();

        T_ij_y = (vecC.adjoint() * Dy.leftCols(bnkl)).transpose() * (veci.adjoint() * Snk);
        T_ij_y += (veci.adjoint() * Dynk).transpose() * (vecC.adjoint() * S.leftCols(bnkl));
        T_ij_y = HF.transpose() * T_ij_y * HF;
        T_ij_y += T_ij_y.transpose().eval();

        T_ij_z = (vecC.adjoint() * Dz.leftCols(bnkl)).transpose() * (veci.adjoint() * Snk);
        T_ij_z += (veci.adjoint() * Dznk).transpose() * (vecC.adjoint() * S.leftCols(bnkl));
        T_ij_z = HF.transpose() * T_ij_z * HF;
        T_ij_z += T_ij_z.transpose().eval();

        Vector3cd T;
        T(0) = (CI * T_ij_x).trace();
        T(1) = (CI * T_ij_y).trace();
        T(2) = (CI * T_ij_z).trace();

        factor = T.squaredNorm();
        factor /= 4.;
        sig = Sigma(jc.getK(), factor, jc.getIB());

        cout << " CI sigma: " << sig << "\n";
        */
    // successful exit
    auto end = chrono::system_clock::now();

    chrono::duration<double> elapsed_seconds = end - start;
    cout << " CPU time: " << setprecision(5) << fixed << elapsed_seconds.count() << " s\n";
    cout << "=========================================================================="
         << "\n"
         << "\n";

    return EXIT_SUCCESS;
}
