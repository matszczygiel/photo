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

#include <iostream>
#include <fstream>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <stdexcept>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>

#include "rec.h"
#include "constants.h"
#include "harmonics.h"
#include "job_control.h"
#include "two_electron_integrals.h"
#include "iteration.h"
#include "disk_reader.h"
#include "calculator.h"

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[])
{
	auto start = chrono::system_clock::now();

	if (argc != 2)
	{
		cout << " Proper usage: ./rec <input name> \n";
		return EXIT_SUCCESS;
	}

	string input = argv[1];

	Job_control controler;
	Input_data in_data;
	Disk_reader reader;
	Calculator calc;
	Iteration iter;

	cout << scientific;
	cout.precision(5);

	try
	{
		ifstream ifile(input);
		if (!ifile.is_open())
			throw runtime_error("Invalid input file parsed.");

		in_data.read_input(ifile);
		ifile.close();

		controler.read(in_data);
		reader.initialize(controler);
		calc.initialize(controler, reader);
		iter.initialize(controler, reader);
	}
	catch (exception &e)
	{
		cerr << e.what() << "\n";
		return EXIT_FAILURE;
	}

	double energy = calc.energy();
	auto vec_cont = calc.continuum_vec();
	VectorXd vec_ion;
	VectorXd grHF = reader.load_HFv().col(0);

	if (controler.get_compute_ion_state())
		vec_ion = calc.bound_vec();
	else
		vec_ion = grHF;

	try
	{
		iter.set_energy(energy);
		iter.set_starting_vecs(vec_ion, vec_cont);
		iter.load_matrices();
	}
	catch (exception &e)
	{
		cerr << e.what() << "\n";
		return EXIT_FAILURE;
	}

	iter.iterate();

	iter.free_ints();

	auto vecI = iter.get_vecI();
	auto vecC = iter.get_vecC();

/*


	Corr_coef = iter.Corr_coef;
	Coeff_ion = iter.Coeff_ion;

	VectorXcd Cont = iter.get_cont();
	VectorXcd Ion = iter.get_ion();

	/// cross section evaluation

	MatrixXcd Dx_non_k = Dx_matrix.topLeftCorner(b_nk_l, b_nk_l);
	MatrixXcd Dy_non_k = Dy_matrix.topLeftCorner(b_nk_l, b_nk_l);
	MatrixXcd Dz_non_k = Dz_matrix.topLeftCorner(b_nk_l, b_nk_l);

	complex<double> S_kI, S_pI, S_pk, Dx_pI, Dx_kI, Dy_pI, Dy_kI, Dz_pI, Dz_kI;

	S_kI = Cont.dot(S_matrix.leftCols(b_nk_l) * HF_ground);
	S_pI = Coeff_ion.dot(S_non_k * HF_ground);
	Dx_pI = Coeff_ion.dot(Dx_non_k * HF_ground);
	Dy_pI = Coeff_ion.dot(Dy_non_k * HF_ground);
	Dz_pI = Coeff_ion.dot(Dz_non_k * HF_ground);
	Dx_kI = Cont.dot(Dx_matrix.leftCols(b_nk_l) * HF_ground);
	Dy_kI = Cont.dot(Dy_matrix.leftCols(b_nk_l) * HF_ground);
	Dz_kI = Cont.dot(Dz_matrix.leftCols(b_nk_l) * HF_ground);

	double r_pI_sqrt, r_kI_sqrt;

	r_pI_sqrt = norm(Dx_pI) + norm(Dy_pI) + norm(Dz_pI);
	r_kI_sqrt = norm(Dx_kI) + norm(Dy_kI) + norm(Dz_kI);

	S_pk = Ion.dot(S_matrix * Cont);

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

	double pk_R_pk, pk_R_kp, norm_k_sqrt;

	MatrixXcd pp_R(b_l, b_l);
	MatrixXcd p_R_p(b_l, b_l);

	pp_R.setZero();

	for (unsigned n = 0; n < b_l; n++)
		for (unsigned m = 0; m < b_l; m++)
			for (unsigned i = 0; i < b_l; i++)
				for (unsigned j = 0; j < b_l; j++)
					pp_R(n, m) += conj(Ion(i)) * Ion(j) * two_el_matrix(i, j, n, m);

	p_R_p.setZero();

	for (unsigned n = 0; n < b_l; n++)
		for (unsigned m = 0; m < b_l; m++)
			for (unsigned i = 0; i < b_l; i++)
				for (unsigned j = 0; j < b_l; j++)
					p_R_p(i, j) += conj(Ion(i)) * Ion(j) * two_el_matrix(i, n, m, j);

	pk_R_pk = real(Cont.dot(pp_R * Cont));
	pk_R_kp = real(Cont.dot(p_R_p * Cont));
	norm_k_sqrt = real(Cont.dot(S_matrix * Cont));

	double norm_psi = 2 * (norm_k_sqrt + norm(S_pk));

	cout << " Norm of the whole state: " << norm_psi << "\n\n";

	// energy is of the form: E = E_bound + 0.5 k^2 + E_repulsion_corr

	double E_rep_corr = (pk_R_kp + pk_R_pk) / (0.5 * norm_psi);

	cout << " Electron repulsion energy: " << E_rep_corr << "\n\n";

	double factor;
	factor = norm(S_kI) * r_pI_sqrt + norm(S_pI) * r_kI_sqrt;

	double sig;
	sig = Sigma(jc.getK(), factor, jc.getIB());

	double En = Energy(jc.getK(), jc.getIB());

	cout << " System energy :      " << fixed << setprecision(3) << En << "\n";
	cout << " k :                  " << fixed << setprecision(3) << jc.getK() << "\n";
	cout << " Photon energy [eV]:  " << fixed << setprecision(3) << photonEeV(jc.getK(), jc.getIB()) << "\n";
	cout << " Cross section :      " << fixed << setprecision(4) << sig << "\n";
	cout << " \n\n\n";

	// if ( jc.getIfWrite() ) jc.writeRes(sig);

	two_el_matrix.resize(0, 0, 0, 0);

	// TEST OF CI
	// H transformation
	MatrixXd H = HF_matrix.transpose() * H_non_k.real() * HF_matrix;
	MatrixXd CI_matrix = get_CI_coefficients(jc.getFileCI(), b_nk_l);
	double En_CI = 2. * (CI_matrix * H * CI_matrix).trace();

	Two_electron_ints_D two_el_real, temp_two;
	read_two_el_from_binary_slice(two_el_real, jc.getFile2E(), b_nk_l, false);

	TensorMap<Tensor<double, 2>> CI_tens(CI_matrix.data(), b_nk_l, b_nk_l);
	TensorMap<Tensor<double, 2>> HF_tens(HF_matrix.data(), b_nk_l, b_nk_l);

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

	auto HF_energies = get_HF_energies(jc.getFileHFEnergy(), b_nk_l);

	cout << " CI energy: " << En_CI << "\n";

	// Dipole moment
	MatrixXcd T_ij_x, T_ij_y, T_ij_z;
	T_ij_x = (Cont.adjoint() * Dx_matrix.leftCols(b_nk_l)).transpose() * (Coeff_ion.adjoint() * S_non_k);
	T_ij_x += (Coeff_ion.adjoint() * Dx_non_k).transpose() * (Cont.adjoint() * S_matrix.leftCols(b_nk_l));
	T_ij_x = HF_matrix.transpose() * T_ij_x * HF_matrix;
	T_ij_x += T_ij_x.transpose().eval();

	T_ij_y = (Cont.adjoint() * Dy_matrix.leftCols(b_nk_l)).transpose() * (Coeff_ion.adjoint() * S_non_k);
	T_ij_y += (Coeff_ion.adjoint() * Dy_non_k).transpose() * (Cont.adjoint() * S_matrix.leftCols(b_nk_l));
	T_ij_y = HF_matrix.transpose() * T_ij_y * HF_matrix;
	T_ij_y += T_ij_y.transpose().eval();

	T_ij_z = (Cont.adjoint() * Dz_matrix.leftCols(b_nk_l)).transpose() * (Coeff_ion.adjoint() * S_non_k);
	T_ij_z += (Coeff_ion.adjoint() * Dz_non_k).transpose() * (Cont.adjoint() * S_matrix.leftCols(b_nk_l));
	T_ij_z = HF_matrix.transpose() * T_ij_z * HF_matrix;
	T_ij_z += T_ij_z.transpose().eval();

	Vector3cd T;
	T(0) = (CI_matrix * T_ij_x).trace();
	T(1) = (CI_matrix * T_ij_y).trace();
	T(2) = (CI_matrix * T_ij_z).trace();

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
