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
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>


#include "rec.h"
#include "fun.h"
#include "constants.h"
#include "harmonics.h"
#include "job_control.h"
#include "two_electron_integrals.h"
#include "iteration.h"

using namespace std;
using namespace Eigen;

using Cdouble = std::complex<double>;


JobControl jc;


int main( int argc, char* argv[]) {
	// put a battery in the clock
	auto start = chrono::system_clock::now();

	if ( argc != 2 ) {
		cout << " Proper usage: ./rec <input name> \n";
		return EXIT_SUCCESS;
	}

	string input = argv[1];
	bool error = false ;
	jc.readInput(input);
	jc.calcualteBasisFunctNumber();
	jc.print();

	const unsigned b_l          = jc.getBL();
	const unsigned b_l_sqrt     = jc.getBLsqrt();
	const unsigned b_nk_l       = jc.getBNKL();
	const unsigned b_nk_l_sqrt  = jc.getBNKLsqrt();
	const unsigned b_k_l        = jc.getBKL();

	cout << scientific;
	cout.precision(5);

	///load integrals to memory

	//one-electron integrals
	Cdouble * one_el;
	Get_1E(jc.getFile1E(), one_el, b_l_sqrt, error, false);

	MatrixXcd S_matrix  = Map<MatrixXcd>(one_el, b_l, b_l).transpose();
	MatrixXcd H_matrix  = Map<MatrixXcd>(&one_el[3*b_l_sqrt], b_l, b_l).transpose();
	MatrixXcd Dx_matrix = Map<MatrixXcd>(&one_el[4*b_l_sqrt], b_l, b_l).transpose();
	MatrixXcd Dy_matrix = Map<MatrixXcd>(&one_el[5*b_l_sqrt], b_l, b_l).transpose();
	MatrixXcd Dz_matrix = Map<MatrixXcd>(&one_el[6*b_l_sqrt], b_l, b_l).transpose();

	delete[] one_el;

	//two-electron integrals

	Tensor<Cdouble, 4> two_el_matrix;
	cout.precision(10);
	read_two_el_from_binary(two_el_matrix, jc.getFile2E(), b_l, false);


	/// import HF orbitlas

	double * HF_coef;
	Get_HF_data(jc.getFileHF(), HF_coef, b_nk_l_sqrt, error, false);

	VectorXd HF_ground  = Map<VectorXd>(HF_coef, b_nk_l);
	MatrixXd U_matrix   = Map<MatrixXd>(&HF_coef[b_nk_l], b_nk_l, b_nk_l-1);
	MatrixXd HF_matrix  = Map<MatrixXd>(HF_coef, b_nk_l, b_nk_l);

	delete[] HF_coef;

	/// import the continuum coefficients and transforamtion to more direct form
	//store in the memory as Gamess order

	Cdouble * cont_coef;
	Prepare_cont_coef( jc.getFileNorm(), cont_coef, jc.getLmax(), b_k_l, b_nk_l, error );

	VectorXcd Cont_coef = Map<VectorXcd>(cont_coef, b_k_l);
	delete[] cont_coef;


	///////Prepare data for the iteration

	///Prepare matrices

	VectorXcd Corr_coef;
	MatrixXcd U_matrix_1eq;
	MatrixXcd H_non_k = H_matrix.topLeftCorner(b_nk_l, b_nk_l);
	MatrixXcd S_non_k = S_matrix.topLeftCorner(b_nk_l, b_nk_l);

	if ( jc.getForceOrthogonality() == true ) {
		Corr_coef.resize(b_nk_l);
		Corr_coef << 1.0, VectorXcd::Zero(b_nk_l-1);

		U_matrix_1eq =  MatrixXcd::Zero(b_l, b_nk_l);

		U_matrix_1eq.col(0) <<  -HF_ground.dot(S_matrix.topRightCorner( b_nk_l, b_k_l) * Cont_coef) * HF_ground , Cont_coef;
		U_matrix_1eq.block(0, 1, b_nk_l, b_nk_l-1) << U_matrix;
	}
	else {
		Corr_coef.resize(b_nk_l+1);
		Corr_coef << 1.0, VectorXcd::Zero(b_nk_l);

		U_matrix_1eq = MatrixXcd::Zero(b_l, b_nk_l+1);

		U_matrix_1eq.col(0) << VectorXcd::Zero(b_nk_l), Cont_coef;
		U_matrix_1eq.block(0, 1, b_nk_l, b_nk_l) << MatrixXd::Identity(b_nk_l, b_nk_l);
	}

	VectorXcd Coeff_ion;

	if ( jc.getSolveIon() == false) {
		Coeff_ion= HF_ground;
	}
	else {
		cout << "\n Computing starting ion orbital. \n";

		GeneralizedSelfAdjointEigenSolver < MatrixXcd > es;
		es.compute( H_non_k, S_non_k );

		cout <<" Ion orbital energies: \n\n" << es.eigenvalues() << "\n\n";

		Coeff_ion = es.eigenvectors().col(0);

		cout << " Eigenvector:\n" << Coeff_ion << "\n\n";
		cout << " Done. \n\n";
	}

	///Let the iteration begin
	Iteration iter(two_el_matrix);

	iter.H_matrix = H_matrix;
	iter.S_matrix = S_matrix;
	iter.U_matrix_1eq = U_matrix_1eq;

	iter.Corr_coef = Corr_coef;
	iter.Coeff_ion = Coeff_ion;

	iter.initialize(jc);

	iter.iterate();

	Corr_coef = iter.Corr_coef;
	Coeff_ion = iter.Coeff_ion;

	VectorXcd Cont = iter.get_cont();
	VectorXcd Ion = iter.get_ion();

	/// cross section evaluation

	MatrixXcd Dx_non_k = Dx_matrix.topLeftCorner(b_nk_l, b_nk_l);
	MatrixXcd Dy_non_k = Dy_matrix.topLeftCorner(b_nk_l, b_nk_l);
	MatrixXcd Dz_non_k = Dz_matrix.topLeftCorner(b_nk_l, b_nk_l);

	complex<double> S_kI, S_pI, S_pk, Dx_pI, Dx_kI, Dy_pI, Dy_kI, Dz_pI, Dz_kI;

	S_kI  = Cont.dot( S_matrix.leftCols(b_nk_l) * HF_ground);
	S_pI  = Coeff_ion.dot( S_non_k * HF_ground);
	Dx_pI = Coeff_ion.dot( Dx_non_k * HF_ground);
	Dy_pI = Coeff_ion.dot( Dy_non_k * HF_ground);
	Dz_pI = Coeff_ion.dot( Dz_non_k * HF_ground);
	Dx_kI = Cont.dot( Dx_matrix.leftCols(b_nk_l) * HF_ground );
	Dy_kI = Cont.dot( Dy_matrix.leftCols(b_nk_l) * HF_ground );
	Dz_kI = Cont.dot( Dz_matrix.leftCols(b_nk_l) * HF_ground );

	double r_pI_sqrt, r_kI_sqrt;

	r_pI_sqrt = norm(Dx_pI) + norm(Dy_pI) + norm(Dz_pI);
	r_kI_sqrt = norm(Dx_kI) + norm(Dy_kI) + norm(Dz_kI);

	S_pk = Ion.dot( S_matrix * Cont);

	cout << "\n" << " S_pk : " << S_pk << "\n\n";

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

	for(unsigned n =0; n< b_l; n++) for(unsigned m = 0; m< b_l; m++)
		for(unsigned i =0; i< b_l; i++) for(unsigned j =0; j<b_l; j++)
			pp_R(n,m) += conj(Ion(i)) * Ion(j) * two_el_matrix( i, j, n, m );

	p_R_p.setZero();

	for(unsigned n =0; n< b_l; n++) for(unsigned m = 0; m< b_l; m++)
		for(unsigned i =0; i< b_l; i++) for(unsigned j =0; j<b_l; j++)
			p_R_p(i,j) += conj(Ion(i)) * Ion(j) * two_el_matrix( i, n, m, j );

	pk_R_pk = real(Cont.dot(pp_R * Cont));
	pk_R_kp = real(Cont.dot(p_R_p * Cont));
	norm_k_sqrt = real( Cont.dot(S_matrix * Cont));

	double norm_psi = 2 * ( norm_k_sqrt + norm(S_pk) );

	cout << " Norm of the whole state: " << norm_psi << "\n\n";

	// energy is of the form: E = E_bound + 0.5 k^2 + E_repulsion_corr

	double E_rep_corr = (pk_R_kp + pk_R_pk) / (0.5 * norm_psi);

	cout << " Electron repulsion energy: " << E_rep_corr << "\n\n";

	double factor;
	factor = norm(S_kI) * r_pI_sqrt + norm(S_pI) * r_kI_sqrt;

	double sig;
	sig = Sigma(jc.getK(), factor, jc.getIB());

	double En = Energy(jc.getK(), jc.getIB());

	cout << " System energy :      " << fixed << setprecision(3) << En 				     <<"\n";
	cout << " k :                  " << fixed << setprecision(3) << jc.getK() 			     <<"\n";
	cout << " Photon energy [eV]:  " << fixed << setprecision(3) << photonEeV(jc.getK(), jc.getIB()) <<"\n";
	cout << " Cross section :      " << fixed << setprecision(4) << sig 			     <<"\n";
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
    dims_pair_0 = { IndexPair<int>(0, 0), IndexPair<int>(1, 2) };
    dims_pair_1 = { IndexPair<int>(0, 0), IndexPair<int>(1, 1) };
    Tensor<double, 2> temp = CI_tens.contract(two_el_real, dims_pair_0);
    Tensor<double, 0> temp_2 = CI_tens.contract(temp, dims_pair_1);

    En_CI += temp_2(0);

    auto HF_energies = get_HF_energies(jc.getFileHFEnergy(), b_nk_l);

    cout << " CI energy: " << En_CI << "\n";

	// Dipole moment
	MatrixXcd T_ij_x, T_ij_y, T_ij_z;
	T_ij_x =  ( Cont.adjoint() * Dx_matrix.leftCols(b_nk_l) ).transpose() * ( Coeff_ion.adjoint() * S_non_k );
	T_ij_x += ( Coeff_ion.adjoint() * Dx_non_k ).transpose() * ( Cont.adjoint() * S_matrix.leftCols(b_nk_l) );
	T_ij_x = HF_matrix.transpose() * T_ij_x * HF_matrix;
	T_ij_x += T_ij_x.transpose().eval();

	T_ij_y =  ( Cont.adjoint() * Dy_matrix.leftCols(b_nk_l) ).transpose() * ( Coeff_ion.adjoint() * S_non_k );
	T_ij_y += ( Coeff_ion.adjoint() * Dy_non_k ).transpose() * ( Cont.adjoint() * S_matrix.leftCols(b_nk_l) );
	T_ij_y = HF_matrix.transpose() * T_ij_y * HF_matrix;
	T_ij_y += T_ij_y.transpose().eval();

	T_ij_z =  ( Cont.adjoint() * Dz_matrix.leftCols(b_nk_l) ).transpose() * ( Coeff_ion.adjoint() * S_non_k );
	T_ij_z += ( Coeff_ion.adjoint() * Dz_non_k ).transpose() * ( Cont.adjoint() * S_matrix.leftCols(b_nk_l) );
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

	if ( jc.getIfWrite() ) jc.writeRes(sig);




	// successful exit
	auto end = chrono::system_clock::now();
	chrono::duration <double> elapsed_seconds = end - start;
	cout<<" CPU time: " << setprecision(5) << fixed << elapsed_seconds.count() << " s\n";
	cout << "==========================================================================" << "\n" << "\n";

	return EXIT_SUCCESS;
}

